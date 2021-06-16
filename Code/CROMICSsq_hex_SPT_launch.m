function []=CROMICSsq_hex_SPT_launch(rep,cA,nameDataFile)

load (nameDataFile);


tStart=tic;

if cA==1        %cA1. Rcell~Mcell
    CellSizeVar=1;          %1=cell size is variable;  0=cell size constant
    
elseif cA==2    %cA2. Rcell cte
    CellSizeVar=0;          %1=cell size is variable;  0=cell size constant
    Rcell_cte=ones(num_sp,1)*6.2035e-04;%Cell radius (assumed to be constant) [mm]. 
    
elseif cA==3    %cA3. Dmet=0.25Dmet
    CellSizeVar=1;          %1=cell size is variable;  0=cell size constant
    Dmet=.25*Dmet;
    v_sp(:)=0;              %Neglect crowding 
    Rmet(:)=0;

elseif cA==4    %cA4. Dmet=Dmet
    CellSizeVar=1;          %1=cell size is variable;  0=cell size constant
    v_sp(:)=0;              %Neglect crowding 
    Rmet(:)=0;
end

if CellSizeVar==1
    Rcell_cte=zeros(num_sp,1);            %Cell radius (assumed to be constant) [mm]. 
end

Km=(Km==0)+Km;
fwater=((6.7)/2.8)*(1/((v_sp(1)==0)+v_sp(1)))*(1e3); %Fraction of intracellular water in a cell [dimensionless].
Vcell_cte=4*pi*(Rcell_cte.^3)/3;    %Cell volume (assumed to be constant) [mm]. This value is only used if CellSizeVar=0, i.e. R constant


%Create matrixes to store time-dependet variables 
rhotCA_sp=zeros(Ny_ca,Nx_ca,Nz_ca,num_sp,floor(tsim/(nAt*At))+1);       %Matrix-CA of the density distribution for cell species [gDW] at different time.
rhotLBM_met=zeros(Ny_lbm,Nx_lbm,Nz_lbm,num_met,floor(tsim/(nAt*At))+1); %Matrix-LBM of the density distribution for metabolite [mmol] at different time.
phentCA=zeros(Ny_ca,Nx_ca,Nz_ca,num_sp,floor(tsim/(nAt*At))+1);         %Matrix-CA of the phenotype of each cell [dimensionless].
postCA_sp=zeros(Ny_ca,Nx_ca,Nz_ca,floor(tsim/(nAt*At))+1);              %Matrix-CA of the position of each cell [dimensionless].
tvec=zeros(floor(tsim/(nAt*At))+1,1);                                   %Time vector [h].
r2tt_met=zeros(num_met,floor(tsim/(nAt*At))+1);                         %Matrix of the mean squared displacement of the metabolites [mm^2] at different time.
r2tt_sp=zeros(num_sp,floor(tsim/(nAt*At))+1);                           %Matrix of the mean squared displacement of the cells [mm^2] at different time.
VbiotCA_sp=zeros(Ny_ca,Nx_ca,Nz_ca,num_sp,floor(tsim/(nAt*At))+1);      %Matrix-CA of growth rate of the cells [h^-1] at different time.

%ID vectors of the possible directions for LBM
Q=(dim*2)+1;                                    %Posible directions of a particle in LBM.
metlbm=zeros(Nx_lbm*Ny_lbm*Nz_lbm*num_met,1);   %Vector with the type of metabolite species [dimensionless]. 
metlbm(:)=repmat(1:num_met,Nx_lbm*Ny_lbm*Nz_lbm,1);

%Save data at time 0 hr
rhotCA_sp(:,:,:,:,1)=rhoCA_sp(:,:,:,:);
rhotLBM_met(:,:,:,:,1)=rhoLBM_met(:,:,:,:);
phentCA(:,:,:,:,1)=rhoCA_sp(:,:,:,:)*0;
postCA_sp(:,:,:,1)=posCA_sp;
tvec(1)=0; %[hr]

%Create matrix to storage variables 
Vmt=ones(Ny_lbm*Nx_lbm*Nz_lbm,num_met)*(-1e30);
VbioCA_sp=zeros(Ny_ca,Nx_ca,Nz_ca,num_sp);
phen_sp=zeros(Ny_ca,Nx_ca,Nz_ca,num_sp);
frac=zeros(Ny_lbm*Nx_lbm*Nz_lbm,num_met);  

rhoLBM_met=reshape(rhoLBM_met,Ny_lbm*Nx_lbm*Nz_lbm,num_met);
if num_metInert>0
    rho_inert=reshape(rho_inert,Ny_lbm*Nx_lbm*Nz_lbm,num_met);
end

%Compute the matrix with the ID and properties of each cell
[IDpos,Ncells]=CellProperties(posCA_sp,Ny_ca,Nx_ca,Nz_ca,num_sp,rhoCA_sp,ageCA,Q,IDca);

if strfind(subpop_analysis,'YES')==1 %Track subpopulations  
    %Add the subpopulation type
    IDpos=[IDpos,zeros(Ny_ca*Nx_ca*Nz_ca,1)]; 
    IDpos(1:Ncells,6)=rN;
    
    subpoptCA=zeros(Ny_ca,Nx_ca,Nz_ca,num_sp,floor(tsim/(nAt*At))+1);         %Matrix-CA of the phenotype of each cell [dimensionless].
    subpop_sp=zeros(Ny_ca,Nx_ca,Nz_ca,num_sp);  %Matrix-CA to store the subpopulation
    subpop_sp(IDpos(1:Ncells,1))=IDpos(1:Ncells,6)*2+1;
    subpoptCA(:,:,:,:,1)=subpop_sp(:,:,:,:);
end

%Since inert macromolecules (with Rmet>0) are introduced in the system after timeInert,
%for efficiency, we will use SPT (computationally more expensive) only after timeInert
%Then, the volume fraction availiable can be computed with the first term
%of the SPT equation (see main text)
%Compute the matrix with the probability to find available space in each LBM-voxel [dimensionless].
if addEPS==1
    [frac,VavLBM]=PavailableSpace_EPS(fepsCell1,v_met,CellSizeVar,Vcell_cte,vol_lbm,v_sp,fwater,IDpos,Ny_lbm,Nx_lbm,Nz_lbm,num_met,frac_medium,Ncells,frac*0);
else
    [frac,VavLBM]=PavailableSpace(CellSizeVar,Vcell_cte,vol_lbm,v_sp,fwater,IDpos,Ny_lbm,Nx_lbm,Nz_lbm,num_met,frac_medium,Ncells,frac*0);
end

%Parameters used in the LBM simulations
Dmax=max(Dmet);                 %Maximum diffusion coefficient of the metabolites [mm^2 ms^-1].
Prob=Dmet/Dmax;                 %Probability of a metabolite to move to a neighboring LBM-voxel during At_lbm [dimensionless].
time_lbm=At;                    %Time for LBM simulation before updating the cells position (using IbM rules). 

if strfind(pde_method,'LBM')==1     %If we want to use Lattice Boltzman Method to compute the metabolites diffusion 
    aQ1=IDlbm(:);                   %Index vector
    if grid_sh=='s'
        W=1/(2*dim);                    %Weight for the equilibrium function in LBM [dimensionless].
    else 
        W=1/9; 
    end
    At_lbm=(Ax_lbm^2)/(2*dim*Dmax); %Time step for LBM simulations [ms].

else    %If we want to use Crank-Nicholson method to compute the metabolites diffusion 
    if dim==3
        Ax_min=min([Ax_lbm,Ay_lbm,Az_lbm]);
    else
        Ax_min=min([Ax_lbm,Ay_lbm]);
    end
    nCN=50; %To decrease running time, we can use At_lbm larger than the conservative At_lbm=(Ax_min^2)/(2*dim*Dmax), i.e. when nCN=1.
             %A large nCN could cause inaccuracies
    At_lbm=((Ax_min^2)/(2*dim*Dmax))*nCN; %Time step for LBM simulations [ms].
    if At_lbm>At %Correction of the maximum At_lbm: At_lbm should be less than At.  
        At_lbm=At;
    end
    lambda_x=Dmet(metlbm)*(At_lbm)*.5/(Ax_lbm^2);
    lambda_y=Dmet(metlbm)*(At_lbm)*.5/(Ay_lbm^2);
    lambda_z=Dmet(metlbm)*(At_lbm)*.5/(Az_lbm^2);
end

%Vector to store the mean squared displacement (MSD) of metabolites and species
r2t_met=zeros(num_met,1); %MSD for metabolites
r2t_sp=zeros(num_sp,1);   %MSD for species

[indCsat,rD_ind]=type_medium(mediumLBM,num_met,metlbm,rD);

if ~exist('GEMs','var')  %In case Neural Networks will be used to compute the metabolic fluxes, and therefore GEM models are not provided.
    GEMs=[];
end


%Simulation
fprintf(strcat('\nComputing... \n',nameDataFile,'\n\n'))
cont2=2;

for contt=2:floor(tsim/At)+1

    %Innert metabolites are added to the system at time timeInert
    if contt==timeInert
        rhoLBM_met=rhoLBM_met+rho_inert;
        
        %Since inert macromolecules (with Rmet>0) are introduced in the system after timeInert,
        %for efficiency, we will use SPT only after timeInert
        %Compute the matrix with the probability to find available space in each LBM-voxel [dimensionless].
        [dens,dens_int,r,rint,vff]=PavailableSpace_SPT_part1(metaB,CellSizeVar,vol_lbm,v_sp,IDpos,frac_medium,Ncells,Reps1,Rcell_cte,Rmet,Rprot,metlbm,rhoLBM_met,AvN,fepsCell,num_eps,num_met,Ny_lbm,Nx_lbm,Nz_lbm,MMprot);
        [frac,VavLBM]=PavailableSpace_SPT_part2(dens,dens_int,r,rint,Ny_lbm,Nx_lbm,Nz_lbm,num_met,Rmet,metlbm,P_kT,frac*0,vff,vol_lbm);
        if cA==4 || cA==3
            frac=reshape(frac_medium,Ny_lbm*Nx_lbm*Nz_lbm,num_met);
            VavLBM=(frac==0)+frac*vol_lbm;
        end
    end
    
    %Compute the metabolic fluxes and the diffusion of the metabolites
    if contt<timeInert
         if strfind(pde_method,'LBM')==1 
           [rhoLBM_met,flux_bio,flux,IDpos,r2t_met]=cLBM_aTFA(frac_medium,cA,CellSizeVar,Vcell_cte,vol_lbm,v_sp,fwater,FD,Ax_lbm,indCsat,rD_ind,Csat,r2t_met,metlbm,AvN,IDlbm,frac,W,Prob,Q,FBD,aQ1,Nx_lbm,Ny_lbm,Nz_lbm,time_lbm,Vmt,At_lbm,Ncells,num_met,num_sp,rhoLBM_met,VavLBM,maxMcell,IDpos,n,Vmax,Km,VmaxDiff,vshrinkage,num_metInert,contt,cellDist,NNlist,GEMs);

        else %pde_method='CN' 
             [rhoLBM_met,flux_bio,flux,IDpos]=CN_aTFA(frac_medium,cA,CellSizeVar,Vcell_cte,vol_lbm,v_sp,fwater,FD,indCsat,rD_ind,Csat,frac,FBD,Nx_lbm,Ny_lbm,Nz_lbm,time_lbm,Vmt,At_lbm,Ncells,num_met,rhoLBM_met,VavLBM,maxMcell,IDpos,Vmax,Km,VmaxDiff,vshrinkage,num_metInert,contt,a1y,c1y,a1x,c1x,a1z,c1z,b1y,b2y,b1x,b2x,b1z,b2z,ap1x,cp1x,ap1y,cp1y,ap1z,cp1z,CNx,CNy,CNz,lambda_x,lambda_y,lambda_z,n,num_sp,cellDist,NNlist,GEMs);
        end 
    else
        if strfind(pde_method,'LBM')==1 
            [rhoLBM_met,flux_bio,flux,IDpos,r2t_met]=cLBM_aTFA_SPT(P_kT,MMprot,num_eps,fepsCell,Rprot,Rmet,Reps1,metaB,Rcell_cte,frac_medium,cA,CellSizeVar,Vcell_cte,vol_lbm,v_sp,fwater,FD,Ax_lbm,indCsat,rD_ind,Csat,r2t_met,metlbm,AvN,IDlbm,frac,W,Prob,Q,FBD,aQ1,Nx_lbm,Ny_lbm,Nz_lbm,time_lbm,Vmt,At_lbm,Ncells,num_met,num_sp,rhoLBM_met,VavLBM,maxMcell,IDpos,n,Vmax,Km,VmaxDiff,vshrinkage,num_metInert,contt,cellDist,NNlist,GEMs);

        else %pde_method='CN' 
            [rhoLBM_met,flux_bio,flux,IDpos]=CN_aTFA_SPT(P_kT,MMprot,num_eps,fepsCell,Rprot,Rmet,Reps1,metaB,Rcell_cte,frac_medium,cA,CellSizeVar,Vcell_cte,vol_lbm,v_sp,fwater,FD,indCsat,rD_ind,Csat,frac,FBD,Nx_lbm,Ny_lbm,Nz_lbm,time_lbm,Vmt,At_lbm,Ncells,num_met,rhoLBM_met,VavLBM,maxMcell,IDpos,Vmax,Km,VmaxDiff,vshrinkage,num_metInert,contt,a1y,c1y,a1x,c1x,a1z,c1z,b1y,b2y,b1x,b2x,b1z,b2z,ap1x,cp1x,ap1y,cp1y,ap1z,cp1z,CNx,CNy,CNz,lambda_x,lambda_y,lambda_z,n,num_sp,cellDist,NNlist,GEMs);
        end    
    end


    %Compute phenotypes and growth rate for saving later:
    if rem(contt-1,nAt)==0
        
        fprintf ('Iteration: %g / %g\n',contt,floor(tsim/At)+1);
        
        %Identify phenotypes using the metabolic fluxes store in flux
        phen_sp(:)=0;
        for sp=1:num_sp
            ph=ph_matrix{sp};
            aa=find(IDpos(1:Ncells,2)==sp);
            dd=bsxfun(@lt,flux(aa,:),minflux{sp});
            for p=1:num_ph(sp)
                %Identify the metabolite consumption
                phen_sp(IDpos(aa,1))=phen_sp(IDpos(aa,1))+(sum(abs(dd-ph(p,:)),2)==0)*p;
            end
        end  
        
        %Growth rate [h^-1]
        VbioCA_sp(:)=0;
        VbioCA_sp(IDpos(1:Ncells,1))=flux_bio;

        %Biomass [gDW] 
        rhoCA_sp(:)=0;
        rhoCA_sp(IDpos(1:Ncells,1))=IDpos(1:Ncells,3);
                
        %Save data
        rhotCA_sp(:,:,:,:,cont2)=rhoCA_sp(:,:,:,:);
        rhotLBM_met(:,:,:,:,cont2)=reshape(rhoLBM_met,Ny_lbm,Nx_lbm,Nz_lbm,num_met);
        phentCA(:,:,:,:,cont2)=phen_sp(:,:,:,:);
        r2tt_met(:,cont2)=r2t_met;
        r2tt_sp(:,cont2)=r2t_sp;
        postCA_sp(:,:,:,cont2)=posCA_sp;
        tvec(cont2)=(contt-1)*At/(1000*3600); %[hr]
        VbiotCA_sp(:,:,:,:,cont2)=VbioCA_sp(:,:,:,:);
        
        if strfind(subpop_analysis,'YES')==1 %Track subpopulations  
            %Subpopulations in each species
            %subpop_sp=1 (i.e. IDpos(1:Ncells,6)=0 ) for mutant subpopulation e.g. S.enterica methionine-producer and E.coli methionine auxotroph
            %subpop_sp=3 (i.e. IDpos(1:Ncells,6)=1 ) for WT subpopulation
            subpop_sp(:)=0;
            subpop_sp(IDpos(1:Ncells,1))=IDpos(1:Ncells,6)*2+1;  
            subpoptCA(:,:,:,:,cont2)=subpop_sp(:,:,:,:);
        end
        
        cont2=cont2+1;
    end

    %IbM. Spatial distribution of the cells. 
    if Ncells<Nx_ca*Ny_ca*Nz_ca  
        if any(Dsp>0)==1 %Cells grow in planktonic state 
            Pmov=Dsp(IDpos(1:Ncells,2))*2*dim*At/(Ax_ca^2);   %Probability of mobility of species sp if the site ramdomly chosen is empty.  
            move=rand(Ncells,1);
            %Index a11 identifies the cells that can move to a neighbor voxel or will divided (when the cell mass>=Mmaxsp) or die (when the cell mass<=Mminsp) 
            a11=find(IDpos(1:Ncells,3)>=Mmaxsp(IDpos(1:Ncells,2)) | IDpos(1:Ncells,3)<=Mminsp(IDpos(1:Ncells,2)) | move<=Pmov);
        
        else %Cells grow in biofilm
            %Index a11 identifies the cells that will divided (when the cell mass>=Mmaxsp) or die (when the cell mass<=Mminsp)
            a11=find(IDpos(1:Ncells,3)>=Mmaxsp(IDpos(1:Ncells,2)) | IDpos(1:Ncells,3)<=Mminsp(IDpos(1:Ncells,2)));
        end
        a22=length(a11);
        if a22>0
             %Motility and cell division
             if grid_sh=='s'
                 [IDpos,Ncells,posCA_sp,r2t_sp,mediumLBM]=IbMcells(boundaryCA,mediumLBM,r2t_sp,Ncells,IDpos,Ny_ca,Nx_ca,Nz_ca,num_sp,posCA_sp,IDca,rV,Q,a11,a22,dim,Mmaxsp,Mminsp,Ax_ca);
             else 
                 [IDpos,Ncells,posCA_sp,r2t_sp,mediumLBM]=IbMcells_hex(mediumLBM,r2t_sp,Ncells,IDpos,Ny_ca,Nx_ca,Nz_ca,num_sp,posCA_sp,IDca,rV,Q,a11,a22,Mmaxsp,Mminsp,Macc1,Ax_ca);
             end
             [indCsat,rD_ind]=type_medium(mediumLBM,num_met,metlbm,rD);
        end
    end

    %Compute the matrix with the probability to find available space in each LBM-voxel [dimensionless].
    if contt<timeInert
        if addEPS==1
            [frac,VavLBM]=PavailableSpace_EPS(fepsCell1,v_met,CellSizeVar,Vcell_cte,vol_lbm,v_sp,fwater,IDpos,Ny_lbm,Nx_lbm,Nz_lbm,num_met,frac_medium,Ncells,frac*0);
        else
            [frac,VavLBM]=PavailableSpace(CellSizeVar,Vcell_cte,vol_lbm,v_sp,fwater,IDpos,Ny_lbm,Nx_lbm,Nz_lbm,num_met,frac_medium,Ncells,frac*0);
        end
    else
        if cA==1 || cA==2
            [dens,dens_int,r,rint,vff]=PavailableSpace_SPT_part1(metaB,CellSizeVar,vol_lbm,v_sp,IDpos,frac_medium,Ncells,Reps1,Rcell_cte,Rmet,Rprot,metlbm,rhoLBM_met,AvN,fepsCell,num_eps,num_met,Ny_lbm,Nx_lbm,Nz_lbm,MMprot);
            [frac,VavLBM]=PavailableSpace_SPT_part2(dens,dens_int,r,rint,Ny_lbm,Nx_lbm,Nz_lbm,num_met,Rmet,metlbm,P_kT,frac*0,vff,vol_lbm);
        end
    end
end

tcomp1=toc(tStart)/3600; %Total running time

%Save results
str=strcat(ResultsFile,'_cA',num2str(cA),'_Nrep',num2str(rep),'.mat');
save (str);