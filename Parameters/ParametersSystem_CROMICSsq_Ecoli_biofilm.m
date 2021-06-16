%Parameters of the System

clear 


%********************************  ASSUMPTIONS  ***************************************
%1. The cell is a sphere.
%2. The system is divided in cubes for both IbM and LBM. 
%3. In IbM simulation (here defined as CA), each voxel can contains only one cell. Therefore, the side of the voxel (Ax_ca) is equal to the diameter of the cell.   
%4. Lattice Boltzmann method requires a square or cubic lattice, while Crank-Nicholson method can handle a cubic or rectangular lattice.
%5. If at least one of the boundaries boundaryX0_lbm or boundaryXn_lbm is Periodic both boundaries will be taken as Periodic.


%**************************  DATA PROVIDED BY THE USER  ********************************

dim=2;               %Dimensions of the system (Dim=2 for 2D, Dim=3 for 3D). 
                     %Note: only 2D systems are simulated with a hexagonal configuration. 
tsim=8.5*3600*1000;  %Simulation time [ms].
At=1*60*1000;        %Time step [ms]. 
nAt=3;               %The simulation results will be record every nAt*At [ms].  
% pde_method='CN';    %For the semi-implicit Crank-Nicholson, CN can handle non-uniform lattice of rectangular-prism.
pde_method='LBM';    %For the explicit Lattice Boltzmann, LBM can handle square and hexagonal lattices in 2D, and cubic lattices in 3D.
% cellDist=0;        %No cell distinction, cellular behaviour is the same for cells sharing the same voxel. 
cellDist=1;          %Cell distinction, each cell can have different metabolic behaviour even if they share the same voxel. 
grid_sh='s';         %h=Hexagonal grid,  s=Square grid
                     %Only 2D systems are simulated with a hexagonal configuration, 
                     %the maximum number of cells per voxel is equal to 1.
subpop_analysis='NO';%No track subpopulations 
% subpop_analysis='YES';%Track subpopulations 

%Name used to save the dataFile (include the path, it is suggested to use the folder "Data" ):  
DataFile=strcat('./Data/dataForCROMICSsq_Ecoli_biofilm.mat')
%Name that will be used to save the results (include the path, it is suggested to use the folder "Results" ):   
ResultsFile=strcat('./Results/resultsCROMICSsq_Ecoli_biofilm');


%Parameters CA
distanceX=0.1901;   %[mm].
distanceY=distanceX;%[mm].
distanceZ=.0019;    %[mm].
Ax_ca=0.0019;       %Side length of CA-voxel [mm], equivalent to 2 time the maximum cells radius.
Ay_ca=Ax_ca;
Az_ca=Ax_ca;        %We assume a monolayer
%We assume cubic LBM-voxels, then nx = ny = nz and Ax_lbm = Ay_lbm = Az_lbm. 
%Note: in the hexagonal configuration, the maximum number of cells per voxel is equal to 1, i.e. nx=ny=nz=1.
nx=1;               %Ratio Ax_lbm/Ax_ca, i.e. the side length of LBM-voxel Ax_lbm will be divided in "nx" times. 
ny=1;               %Ratio Ay_lbm/Ay_ca, i.e. the side length of LBM-voxel Ay_lbm will be divided in "ny" times. 
nz=1;
Ax_lbm=nx*Ax_ca;    %Side length of LBM-voxel [mm].
Ay_lbm=ny*Ay_ca;    %Side length of LBM-voxel [mm].
Az_lbm=nz*Az_ca;    %Side length of LBM-voxel [mm].
Nx_lbm=round(distanceX/Ax_lbm); %Number of LBM-voxels in the coordinate x of the system. 
Ny_lbm=round(distanceY/Ay_lbm); %Number of LBM-voxels in the coordinate y of the system. 
Nz_lbm=round(distanceZ/Az_lbm); %Number of LBM-voxels in the coordinate z of the system. If the system is 2D, then Nz_ca=1.
Nx_ca=Nx_lbm*nx;                %Number of CA-voxels in the coordinate x of the system. The size of the system is proportional to [Nx_ca*Ax_ca] x [Ny_ca*Ax_ca]. 
Ny_ca=Ny_lbm*ny;                %Number of CA-voxels in the coordinate y of the system.
Nz_ca=Nz_lbm*nz;
vol_lbm=Ax_lbm*Ay_lbm*Az_lbm;   %Volume of a LBM-voxel [mm^3].

%When a 2D system is simulated, we assume a monoloyer of cells 
if dim==2 && Nz_ca==1   %Correction of Ny_lbm to 2D systems
    Nz_lbm=1;
end

%Parameters of the cell species simulated
%All vectors must be [number of species x parameter value]
metaB=1;                %Number of cells in one metabacteria, equivalent to Ax_ca/maximumCellRadius, maximumCellRadius=0.0019mm. 
Mmaxsp=[1.172e-12].*metaB; %Maximum mass of a cell or metabacteria [gDW]. The vector Mmaxsp contains the mass of all species simulated.
Mminsp=[0];             %Minimum mass of a cell or metabacteria [gDW]. 
v_sp=[3070.5];          %Specific volume [mm^3 gDW^-1]
Dsp=[0];                %Diffusion coefficient of the cell [mm^2 ms^-1]. 
vshrinkage=0.016;       %Cell shrinkage rate [h^-1].
mu=(0.489e-12)*metaB;   %Mean of cell mass [gDW]. 
sigma=(0.132e-12)*metaB;%Standard deviation of cell mass [gDW].
%For biofilm simulations, we can choose a random initial thickness of the biofilm:
nCell=floor(4+(14-4)*rand(1,Nx_ca)); %Initial number of cells in the 'y' axis 
%Active uptake substrate (the uptake rate follows the Michaelis-Menten equation [glucose, O2, acetate, CO2]):
Km=[0.015 0 0];   %Matrix [number of species x number of metabolites] of Michaelis-Menten constant of all metabolites and species [mM].
Vmax=[10 0 0];                  %Matrix [num_sp x num_met] of the maximum uptake rate of all metabolites and species [mmol/gDW h].
%Passive uptake substrate (simple diffusion):
VmaxDiff=[0 15 17]; 
rV=[1];  %Volume ratio between the mother and daughter cell, i.e. Vmother/Vdaughter once the cellular division takes place [dimensionless].

%We will use Neural Networks to estimate the metabolic fluxes of the species 
NNlist={'NN15x15_Ecoli_glc_paTFA'};

% %***************
% %If we want to use directly Genome-scale models, then you have to provide
% %the following info:
% load Ecoli_glc_ac_o2_aTFA.mat  %sp1
% GEMs.sp1.model=ttmodel;
% %Names of the exchange reactions whose metabolic flux will be computed biomass_rxn name should be the last one in the vector.
% rxnNames={'DM_lcts_e';'DM_o2_e';'DM_met-L_e';'DM_ac_e';'DM_gal_e'}; %reactions name must be the same for all species.
% [RxnF1,RxnR1]=FandBrxns(ttmodel,rxnNames);
% GEMs.sp1.indexMetF=RxnF1;
% GEMs.sp1.indexMetR=RxnR1;
% rxnNames={'Ec_biomass_iJO1366_core_53p95M'};
% [RxnF1,RxnR1]=FandBrxns(ttmodel,rxnNames);
% GEMs.sp1.indexBioF=RxnF1;
% GEMs.sp1.indexBioR=RxnR1;
% clear model;
%  
% %File modelTFA.m computes the metabolic fluxes using TFA and the GEMs provided  
% NNlist={'modelTFA'};
% %***************


%Parameters of the substrates/metabolites simulated, including the possible metabolic products 
%All vectors must be [num_met x 1]. Note: Place inert metabolites at the end of the vector. 
AvN=6.022140857e23;             %Avogadro constant [molecules mol^-1].
MMmet=[180.15;31.998;59.044];   %Molecular mass of metabolites [g mol^-1]. [glucose, O2, acetate]
Mmet=MMmet/AvN;                 %Mass of a metabolite molecule [g]. 
v_met=730;                      %Specific volume of metabolites [mm^3 g^-1].
Rmet=(MMmet*v_met/(4*pi*AvN/3)).^(1/3); %Radius of a metabolite molecule [mm].
Dmet=[6.7e-7;2e-6;1.21e-6];%Diffusion coefficient of metabolites [mm^2 ms^-1]. [lactose, O2, methionine, acetate, galactose]
num_metInert=0;                 %Indicate the number of inert metabolites.

num_sp=length(Mmaxsp);          %Number of species simulated.
num_met=length(Mmet);           %Number of metabolites simulated.


addEPS=false;   %Add EPS associated to biomass, if addEPS=true provide the following info
% %Data for the EPS (Polysaccharides,proteins,DNA) associated to the biomass 
% %Intracellular components:
% MMprot=[72000]; %Molecular weigth of intracellular macromolecules (proteins) of 72 kDa [g mol^-1]
% Rprot=(MMprot*v_met/(4*pi*AvN/3)).^(1/3); %[mm].
% %EPS: 
% fepsCell1=[0;.1/.9]; %Ratio g_eps/g_cell. Matrix [num_sp x num_eps]
% MMeps=[250000000];  %Average molecular weight of Pls, proteins, and DNA [g mol^-1]. Vector [1 x num_eps]
% v_eps=v_sp(2)*3;    %Specific volume of eps. Vector [1 x num_eps]
% Reps1=(MMeps*v_eps/(4*pi*AvN/3)).^(1/3); %[mm]. Vector [1 x num_eps]
% fepsCell=(fepsCell1./MMeps).*AvN; %Vector [num_sp x num_eps]
% num_eps=length(MMeps);


P_kT=101325/((1.38064852*10^-23)*(298.15)*(10^9)); %P/kT. Pressure [Pa], Boltzmann constant [m^2 kg s^-2 K^-1], Temperature [K].



%Diffusion of species and metabolites in diferente medium
%Ratio between the diffusion coefficient at medium1 and water, 
%i.e. rD=Dsp(medium1)/Dsp(water).
% rD=zeros(num_met,4); 
% rD(:,1)=0;    %Medium1=Solids. Metabolites cannot diffuse through solids
% rD(:,2)=1;    %Medium1=Polysaccharide matrix.
% rD(:,3)=1;    %Medium1=Agar.
% rD(:,4)=0;    %Medium1=Air. Only O2 can diffuse to air.
rD=[0 1 1 0;    %glucose
    0 1 1 1;    %O2
    0 1 1 0];   %acetate


%Vector Cinert contains the boundary concentration of all the inert chemical compounds. 
%In case they are supplied to system at a specific time.
timeInert=0;    %Time at which the inert chemical compounds are added to the system [ms].
Cinert=[0];     %Concentration of the inert compounds that will be supplied at the top boundary of the system [mM]. 
                %If there is not inert metabolites then write Cinert=[0].
                %The concentration will be changed to [mmol/mm^3] below.

%Vector Cmet contains the initial concentraton of all metabolites. Including inert metabolites if they are present  
Cmet=[0.5,0.21,0];%Initial concentration of the metabolites [mM]. [glucose, O2, acetate]. The concentration will be changed to [mmol/mm^3] below.

%Assuming that voxels containing air supply O2 to the colony. We assume that a layer 
%of water cover the colny, then the O2 concentration in air is equal to concentration 
%of the dissolved gas in water
Csat=0.21;    %Concentration of the dissolved gas in water [mM]. 

%Phenotypes of interest
%Provide the matrix with the uptake flux threshold for each metabolites that will be used to identify the phenotypes
%The matrix of each species must be [num_phenotypes x num_sp]. Note: Use the same metabolites order as definend in MMmet.
%Each row of the matrix indicates a different phenotype. 
%For example if the metabolite order is [lactose, O2, methionine,acetate,galactose bio]
%and we want to identify the consumption of lactose and O2 using a flux threshold of -1e-4 mmol gDW^-1 h^-1, 
%i.e. lactose_flux<-1e-4 mmol gDW^-1 and O2<-1e-6 mmol gDW^-1,
%then the matrix row should be: [-1e-4 -1e-6 0 0 0]. 
%Use only one flux threshold for each metabolite, i.e. -1e-6 works for O2 in E. coli.
%E. coli [glucose, O2, acetate]
p_sp1=[-1e-4 -1e-4   0;     %phen_sp=1 for glucose + O2 --> biomass
         0   -1e-4 -1e-4;   %phen_sp=2 for acetate + O2 --> biomass
       -1e-4   0     0;     %phen_sp=3 for glucose --> biomass
         0     0     0 ;  %phen_sp=4 for inactive cell, then cell will shrink (i.e. flux_bio=-vshirinkage)
       -1e-4 -1e-4 -1e-4];   %phen_sp=5 for glucose + acetate + O2 --> biomass
%Combine the species matrixes in one:
ph_matrix={p_sp1};


%Porosity if the medium. Voxels can contain (stationary) solids to simulate a porous medium
porosity=1;                     %Porosity of the medium, i.e. Available_space/Total_space


%Boundary conditions for the cells in the CA-lattice
boundaryCA.X0='B';              %Boundary at top of the system-CA, i.e. for cells [dimensionless].
                                %P=Periodic   B=Bounceback   D=Dirichlet
boundaryCA.Xn='B';             
boundaryCA.Y0='B';               
boundaryCA.Yn='B';
boundaryCA.Z0='B';
boundaryCA.Zn='B';


%Boundary conditions for the metabolites in the LBM-lattice. [glucose, O2, acetate, CO2]
boundaryLBM.X0={'B','B','B'};   %Boundary at top of the system-LBM, i.e. for metabolites [dimensionless]. 
CmetDb.X0=zeros(size(Mmet));    %If the boundary is 'Dirichlet', then provide the fixed concentration for each metabolite in [g/mm^3]. 
                                %This value will be considered only if boundaryX0_lbm='Dirichlet'
boundaryLBM.Xn={'B','B','B'};             
CmetDb.Xn=zeros(size(Mmet));	
boundaryLBM.Y0={'B','B','B'};               
CmetDb.Y0=zeros(size(Mmet));        %[mM]. The concentration will be changed to [mmol/mm^3] below.
boundaryLBM.Yn={'B','B','B'};
CmetDb.Yn=zeros(size(Mmet));
boundaryLBM.Z0={'B','B','B'};
CmetDb.Z0=zeros(size(Mmet));
boundaryLBM.Zn={'B','B','B'};
CmetDb.Zn=zeros(size(Mmet));


%*********  Initial cell distribution in the system for Ecoli_biofilm example  *********
%Medium data. Voxels can contain different medium 1=solids (porous medium), 2=polysaccharide matrix, 3=agar, 4=air. 
medium=ones(Ny_ca,Nx_ca,Nz_ca)*2;   %We assume that all voxels have "polysaccharide matrix" or water where both species and metabolites can diffuse freely.
npore=ceil(Nx_ca*Ny_ca*Nz_ca*(1-porosity));
xyz=randperm(Ny_ca*Nx_ca*Nz_ca);
medium(xyz(1:npore))=1; %Voxel is occupied by solids.

%Initial cells position
rhoCA_sp=zeros(Ny_ca,Nx_ca,Nz_ca,num_sp);       %Matrix-CA of the density distribution for cell species [g].
for j=1:Nx_ca
    rhoCA_sp(Ny_ca-nCell(j):Ny_ca,j,:,1)=(medium(Ny_ca-nCell(j):Ny_ca,j,:)>1).*normrnd(mu,sigma,[nCell(j)+1,1]); %Biomass from a normal distribution
end


% %***************
% %In a well mixed system, the cells are randomly located in the whole system 
% %(for the Brownian motion experiment): 
% rpos=randperm(Nx_ca*Ny_ca*Nz_ca);
% rh=rhoCA_sp(find(rhoCA_sp));
% rhoCA_sp=rhoCA_sp*0;
% rhoCA_sp(rpos(1:length(rh))')=rh;
%***************





%*************  PARAMETERS ESTIMATION BASED ON DATA PROVIDED BY THE USERS  **************

maxMsp=max(Mmaxsp);
sp=find(Mmaxsp==maxMsp);
Rsp=((3*v_sp(sp(1))*maxMsp/(4*pi)).^(1/3)); %Maximum radius of a cell [mm].

%Compute the matrixes to identify species phenotypes 
num_ph=zeros(num_sp,1);
for sp=1:num_sp
    ph=ph_matrix{sp};
    ph_matrix{sp}=(ph<0);
    num_ph(sp)=size(ph,1);
    minflux{sp}=min(ph);     
end

%Conversion of the metabolite concentration from [mM] to [mmol/mm^3]
CmetDb.Xn=CmetDb.Xn*1e-6;   %Metabolite concentration [mmol/mm^3].
CmetDb.X0=CmetDb.X0*1e-6;
CmetDb.Yn=CmetDb.Yn*1e-6;
CmetDb.Y0=CmetDb.Y0*1e-6;
CmetDb.Zn=CmetDb.Zn*1e-6;
CmetDb.Z0=CmetDb.Z0*1e-6;
Cmet=Cmet*1e-6;             %Initial concentration of the metabolites [mmol/mm^3]. 
Csat=Csat*1e-6;             %Oxygen saturation [mmol/mm^3]



%Position matrix containing cells species and inert solids (in porous media) 
posCA_sp=zeros(Ny_ca,Nx_ca,Nz_ca);         %Cell position matrix-CA of all the species.
socc=find(medium==1 | medium==3);          %Voxels occupied by inert solids (porous medium), i.e. medium==1, or by agar, i.e. medium==3, and therefore these voxels cannont contain cells.
posCA_sp(socc)=num_sp+1;                   %Voxel is occupied by inert solids (porous medium).

n=max([nx ny nz]);  %n>1 indicates if one LBM-voxel can allocate more than one cells (or CA-voxel).
                    %n=1 indicates a fine discretization where one LBM-voxel can allocate only one cell (or CA-voxel).  

                    
%Compute the matrixes and vectors used by CN or LBM methods 
if grid_sh=='h'
   %Hexagonal grid
    [Q,IDlbm,IDca,FD,FBD,vol_ca,vol_lbm,posCA_sp,frac_medium,frac_cells,frac_intcells,Macc]=HexLatticeLBM(Nx_lbm,Ny_lbm,Nz_lbm,num_met,num_sp,Ax_ca,n,Nx_ca,Ny_ca,Nz_ca,dim,rhoCA_sp,posCA_sp,boundaryLBM,boundaryCA,CmetDb,v_sp,v_met,medium);

else
    %Square grid
    [Q,IDlbm,IDca,FD,FBD,vol_ca,vol_lbm,posCA_sp,frac_medium,frac_cells,frac_intcells,a1y,c1y,a1x,c1x,a1z,c1z,b1y,b2y,b1x,b2x,b1z,b2z,ap1x,cp1x,ap1y,cp1y,ap1z,cp1z,CNx,CNy,CNz]=SqLatticeCNandLBM(Nx_lbm,Ny_lbm,Nz_lbm,num_met,num_sp,Ax_ca,Ay_ca,Az_ca,Ax_lbm,Ay_lbm,Az_lbm,nx,ny,nz,Nx_ca,Ny_ca,Nz_ca,dim,rhoCA_sp,posCA_sp,boundaryLBM,boundaryCA,CmetDb,v_sp,v_met,medium);
end


                    
%Matrix-LBM with the medium data. Voxels can contain different medium 1=solids (porous medium), 2=polysaccharide matrix, 3=agar, 4=air. 
%mediumLBM is computed from the matrix-CA medium  
mediumLBM=ones(Ny_lbm,Nx_lbm,Nz_lbm)*4;   %All voxels have "air" by default.
for med=3:-1:2
    aa=find(medium==med);
    mediumLBM(IDca(aa,Q+1))=med;
end


%Compute the matrix for the LBM-lattice containing the amount of metabolites in each LBM-voxel
%The amount of metabolite is proportional to the concentration Cmet and the
%available space in the LBM-voxel (i.e. the space not occupied by cells)
rhoLBM_met=zeros(Ny_lbm,Nx_lbm,Nz_lbm,num_met); %Matrix-LBM with the amount of metabolites met in each LBM-voxel[mmol].
for met=1:num_met
    if met==2
        rhoLBM_met(Nx_lbm*Ny_lbm*Nz_lbm*(met-1)+1:Nx_lbm*Ny_lbm*Nz_lbm*met)=(mediumLBM(:)==2).*Cmet(met)*vol_lbm.*(frac_medium(Nx_lbm*Ny_lbm*Nz_lbm*(met-1)+1:Nx_lbm*Ny_lbm*Nz_lbm*met)'-frac_intcells(:));
    else
        rhoLBM_met(Nx_lbm*Ny_lbm*Nz_lbm*(met-1)+1:Nx_lbm*Ny_lbm*Nz_lbm*met)=(mediumLBM(:)==2).*Cmet(met)*vol_lbm.*(frac_medium(Nx_lbm*Ny_lbm*Nz_lbm*(met-1)+1:Nx_lbm*Ny_lbm*Nz_lbm*met)'-frac_cells(:));
    end
end

Csat=Csat*vol_lbm;  %Oxygen saturation [mmol] in a LBM-voxel.


%Compute a matrix for the LBM-lattice containing the amount of inert metabolites that will be added to the system from the top boundary. 
if num_metInert>0
    rho_inert=zeros(Ny_lbm,Nx_lbm,Nz_lbm,num_met);
    for met=num_met-num_metInert+1:num_met
        %Inert chemical compounds
        rho_inert(1,:,:,met)=Cinert(met-(num_met-num_metInert))*(1e-6)*vol_lbm;
    end
end

ageCA=rhoCA_sp>0;  %Matrix-Ca with the age of the cells (i.e. )
r2t_met=zeros(num_met,1);
r2t_sp=zeros(num_sp,1);


%If a cell grows but it couldn't divide because there are no free voxels (e.g. it the cell is surounded by solid walls). 
%Then the cell will not be able to grow or consumed nutrient, when its mass is equal or higher than maxMcell
maxMcell=0.9*vol_ca./v_sp(1);  %Maximum volume that a cell can occupy in a voxel-CA. Assumed 

%Delete variables not need for the simulation
var={'nCell','medium','npore','xyz','frac_cells','frac_intcells'};
clear(var{:});


save (DataFile)
