function [rhoLBM_met,flux_bio,flux,IDpos,r2t_met]=cLBM_aTFA_SPT(metaB,FD,cA,num_sp,indCsat,rD_ind,Csat,r2t_met,metlbm,AvN,Ax_lbm,IDlbm,frac,W,Prob,Q,FBD,aQ1,Nx_lbm,Ny_lbm,Nz_lbm,time_lbm,Vmt,At_lbm,Ncells,num_met,rhoLBM_met,VavLBM,maxMcell,IDpos,n,Vmax,Km,VmaxDiff,vshrinkage,num_metInert,cellDist,frac_medium,vol_lbm,num_eps,dens,dens_int,r,rint,Rmet,P_kT,vff,NNlist,GEMs,varargin)

% if cA==1 || cA==2
%     Biofilm=(frac(:)<1);
% else
    Biofilm=frac(:,1)*0;
    Biofilm(IDpos(1:Ncells,4))=1;
    Biofilm=repmat(Biofilm,num_met,1);
% end

%Spatial distribution of metabolites and metabolic fluxes
if cellDist==0 && n>1 %Cells of the same species under the same environmental concentration behave in the same way.
    posBio=IDpos(1:Ncells,4)+Nx_lbm*Ny_lbm*Nz_lbm*(IDpos(1:Ncells,2)-1);
    [~,ia,ic]=unique(posBio,'stable');
    Ncells1=length(ia);
else %Cells are model individually because they are different from ecah other 
    ia=1:Ncells;
    ic=[];
    Ncells1=Ncells;
    posBio=[];
end

%Create a position vector of the cells in LBM-grid     
lbmgrid=repmat(IDpos(ia,4),1,num_met)+repmat((0:num_met-1)*Nx_lbm*Ny_lbm*Nz_lbm,Ncells1,1);
tupt=1; %Frequency with which the metabolic fluxes are evaluated (i.e. once every "tupt" LBM-iterations, see below). 
%To speed up tupt can be selected greater than 1
tprev=0;

 
for t=1:floor(time_lbm/At_lbm)
    
    %Compute metabolic fluxes using either TFA or NN  at every time step At_lbm
    if rem(t,tupt)==0 || t==1 
        tt=t-tprev;
        tprev=t;
             
        %(Re-)compute exchange fluxes and growth rate
        [Vmt,flux_bio,flux]=MetabolicFluxes(Vmt*0,At_lbm*tt,Ncells,rhoLBM_met,VavLBM,lbmgrid,maxMcell,IDpos,n,Vmax,Km,VmaxDiff,vshrinkage,num_metInert,rhoLBM_met*0,NNlist,num_sp,cellDist,posBio,ia,ic,GEMs);

        %Update the amount of metabolites
        rho=rhoLBM_met+(Vmt*At_lbm*tt);
        rho(:)=(rho(:)>0).*rho(:);

        %Update the cell mass
        IDpos(1:Ncells,3)=IDpos(1:Ncells,3)+flux_bio.*IDpos(1:Ncells,3)*At_lbm*tt/(1000*3600); %[gDW]
    
    else
        rho=rhoLBM_met;
    end   
    
%%  cLBM
%     if t==1
        %For efficiency we assume that the cell size is constant during
        %time_lbm, thus the probability P to find available space is
        %computed once every time_lbm
        if cA==1 || cA==2 %Only these crowding assumptions consider that cells have volume
            %Here, we compute the volume fraction available in each LBM-voxel 
            %using Scaled Particle Theory, which depends on the size and 
            %abundance of the cells, macromolecules, and metabolites
            [dens,dens_int,r,rint,vff]=PavailableSpace_SPT_part1(metaB,CellSizeVar,vol_lbm,v_sp,IDpos,frac_medium,Ncells,Reps1,Rcell_cte,Rmet,Rprot,metlbm,rhoLBM_met,AvN,fepsCell,num_eps,num_met,Ny_lbm,Nx_lbm,Nz_lbm,MMprot);
            [frac,VavLBM]=PavailableSpace_SPT_part2(dens,dens_int,r,rint,Ny_lbm,Nx_lbm,Nz_lbm,num_met,Rmet,metlbm,P_kT,frac*0,vff,vol_lbm);
        end
     
        P=[frac(:);0;1].*[rD_ind(:);1;1];
        F2_7=W.*repmat(Prob(metlbm),1,Q).*P(IDlbm);
        F2_7(1:Nx_lbm*Ny_lbm*Nz_lbm*num_met)=1-sum(F2_7(:,2:end),2);
        FBD1=FBD.*P(IDlbm(:,1));
%     end
 
    
    
    %If a neighbouring voxel is occupied by air, we assume that the amount O2 in such voxel is equal to the oxygen-saturated water (Csat) 
    rho(Nx_lbm*Ny_lbm*Nz_lbm+indCsat)=Csat*P(Nx_lbm*Ny_lbm*Nz_lbm+indCsat);
    
    %Streaming step
    aQ2=repmat(rho(IDlbm(:,1)),1,Q).*F2_7;
    rho1=accumarray(aQ1,aQ2(:));

    %Update the amount of metabolites in the Dirichlet boundaries, i.e. constant concentration       
    rhoLBM_met(:)=rho1(1:Nx_lbm*Ny_lbm*Nz_lbm*num_met).*(FD==0)+FBD1;
    
    
%%  Compute the mean square displacement of the metabolites in the biofilm   
    r2met=accumarray(metlbm,sum(aQ2(:,2:end).*Biofilm,2)*(AvN/1000)*(Ax_lbm^2)); %Mean square displacement of the metabolites [mm^2].
    molec=accumarray(metlbm,sum(aQ2(:,1:end).*Biofilm,2)*AvN/1000);              %Number of molecules of each metabolite [molecules].
    r2t_met=r2t_met+(r2met./((molec(:)==0)+molec(:)));
end
 