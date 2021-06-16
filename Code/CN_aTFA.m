function [rhoLBM_met,flux_bio,flux,IDpos]=CN_aTFA(frac_medium,cA,CellSizeVar,Vcell_cte,vol_lbm,v_sp,fwater,FD,indCsat,rD_ind,Csat,frac,FBD,Nx_lbm,Ny_lbm,Nz_lbm,time_lbm,Vmt,At_lbm,Ncells,num_met,rhoLBM_met,VavLBM,maxMcell,IDpos,Vmax,Km,VmaxDiff,vshrinkage,num_metInert,contt,a1y,c1y,a1x,c1x,a1z,c1z,b1y,b2y,b1x,b2x,b1z,b2z,ap1x,cp1x,ap1y,cp1y,ap1z,cp1z,CNx,CNy,CNz,lambda_x,lambda_y,lambda_z,n,num_sp,cellDist,NNlist,GEMs)

%Spatial distribution of metabolites and metabolic fluxes
if cellDist==0 && n>1 %Cells of the same species under the same environmental concentration behave in the same way.
    posBio=IDpos(1:Ncells,4)+Nx_lbm*Ny_lbm*Nz_lbm*(IDpos(1:Ncells,2)-1);
    [~,ia,ic]=unique(posBio,'stable');
    Ncells1=length(ia);
else %Cells are model individually because they are different from each other 
    ia=1:Ncells;
    ic=[];
    Ncells1=Ncells;
    posBio=[];
end
%Create a position vector of the cells in LBM-grid     
lbmgrid=repmat(IDpos(ia,4),1,num_met)+repmat((0:num_met-1)*Nx_lbm*Ny_lbm*Nz_lbm,Ncells1,1);


for t=1:round(time_lbm/At_lbm)

    %Compute metabolic fluxes using either TFA or NN  at every time step At_lbm
    %(Re-)compute exchange fluxes and growth rate
    [Vmt,flux_bio,flux]=MetabolicFluxes(Vmt*0,At_lbm,Ncells,rhoLBM_met,VavLBM,lbmgrid,maxMcell,IDpos,n,Vmax,Km,VmaxDiff,vshrinkage,num_metInert,rhoLBM_met*0,NNlist,num_sp,cellDist,posBio,ia,ic,GEMs);

    %Update the amount of metabolites
    rho=rhoLBM_met+(Vmt*At_lbm);
    rho(:)=(rho(:)>0).*rho(:);
    
    %Update the cell mass
    IDpos(1:Ncells,3)=IDpos(1:Ncells,3)+flux_bio.*IDpos(1:Ncells,3)*At_lbm/(1000*3600); %[gDW]

    
%%  Crank-Nicholson
%     if t==1
        %For efficiency we assume that the cell size is constant during
        %time_lbm, thus the probability P to find available space is
        %computed once every time_lbm
        if cA==1 || cA==2 %Only these crowding assumptions consider that cells have volume
            %Here, we compute the volume fraction available in each LBM-voxel 
            %using Scaled Particle Theory, which depends on the size and 
            %abundance of the cells, macromolecules, and metabolites
            [frac,VavLBM]=PavailableSpace(CellSizeVar,Vcell_cte,vol_lbm,v_sp,fwater,IDpos,Ny_lbm,Nx_lbm,Nz_lbm,num_met,frac_medium,Ncells,frac*0);

        end
        
        P=[frac(:);0;1].*[rD_ind(:);1;1];

        ay=lambda_y.*P(a1y);
        by=lambda_y.*(P(b1y)+P(b2y));
        cy=lambda_y.*P(c1y);
        apy=lambda_y.*P(ap1y);
        cpy=lambda_y.*P(cp1y);
        
        ax=lambda_x.*P(a1x);
        bx=lambda_x.*(P(b1x)+P(b2x));
        cx=lambda_x.*P(c1x);
        apx=lambda_x.*P(ap1x);
        cpx=lambda_x.*P(cp1x);
        
        if Nz_lbm>1
            az=lambda_z.*P(a1z);
            bz=lambda_z.*(P(b1z)+P(b2z));
            cz=lambda_z.*P(c1z);
            apz=lambda_z.*P(ap1z);
            cpz=lambda_z.*P(cp1z);
        end

        %Matrixes for Crank-Nicholson
        %Axis Y:
        Aty=spdiags([-apy,-ay,1+by,-cy,-cpy],[-Ny_lbm+1,-1:1,Ny_lbm-1],Ny_lbm*Nx_lbm*Nz_lbm*num_met,Ny_lbm*Nx_lbm*Nz_lbm*num_met);
        Ay=spdiags([apy,ay,1-by,cy,cpy],[-Ny_lbm+1,-1:1,Ny_lbm-1],Ny_lbm*Nx_lbm*Nz_lbm*num_met,Ny_lbm*Nx_lbm*Nz_lbm*num_met);
        
        %Axis X:
        Atx=spdiags([-apx,-ax,1+bx,-cx,-cpx],[-Nx_lbm+1,-1:1,Nx_lbm-1],Ny_lbm*Nx_lbm*Nz_lbm*num_met,Ny_lbm*Nx_lbm*Nz_lbm*num_met);
        Ax=spdiags([apx,ax,1-bx,cx,cpx],[-Nx_lbm+1,-1:1,Nx_lbm-1],Ny_lbm*Nx_lbm*Nz_lbm*num_met,Ny_lbm*Nx_lbm*Nz_lbm*num_met);
        
        if Nz_lbm>1
            %Axis Z:
            Atz=spdiags([-apz,-az,1+bz,-cz,-cpz],[-Nz_lbm+1,-1:1,Nz_lbm-1],Ny_lbm*Nx_lbm*Nz_lbm*num_met,Ny_lbm*Nx_lbm*Nz_lbm*num_met);
            Az=spdiags([apz,az,1-bz,cz,cpz],[-Nz_lbm+1,-1:1,Nz_lbm-1],Ny_lbm*Nx_lbm*Nz_lbm*num_met,Ny_lbm*Nx_lbm*Nz_lbm*num_met);
        end
        
        FBD1=FBD.*P(CNy);
%     end
    
    %If a neighbouring voxel is occupied by air, we assume that the amount O2 in such voxel is equal to the oxygen-saturated water (Csat)
    rho(Nx_lbm*Ny_lbm*Nz_lbm+indCsat)=Csat*P(Nx_lbm*Ny_lbm*Nz_lbm+indCsat);
        
    sfactor=min(rho.*(rho>0)+(10^10).*(rho==0)); %Scale factor
    rho=rho./sfactor;
      
    %Semi-implicit Crank-Nicholson
    %Axis Y:
    rho(CNy)=Aty\(Ay*rho(CNy));

    %Axis X:
    rho(CNx)=Atx\(Ax*rho(CNx));

    if Nz_lbm>1
        %Axis Z:
        rho(CNz)=Atz\(Az*rho(CNz));
    end

   
    rhoLBM_met=rho.*sfactor;
    rhoLBM_met(:)=rhoLBM_met(1:Nx_lbm*Ny_lbm*Nz_lbm*num_met)'.*(FD==0)+FBD1;
end
  
    