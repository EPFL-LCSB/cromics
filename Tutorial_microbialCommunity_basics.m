%% CROwding-Modeling of In-silico Community Systems (CROMICS)

% Paper:
% 1. Angeles-Martinez and Hatzimanikatis. "Spatio-temporal modeling of the crowding conditions and metabolic variability in microbial communities".

%This tutorial demonstrate the use of CROMICS for the spatio-temporal simulation 
%of microbial communities.

%CROMICS is an individual-based model wherein the effect of the crowding conditions 
%is explicitly incorporated. CROMICS combines methods and techniques such as 
%stoichiometric-based models (in particular thermodynamic flux analysis (TFA) [Henry et al., 2007])
%to compute the metabolic fluxes of each individual cell in the system, and
%scaled particle theory (SPT) [Lebowitz et al., 1965] to estimate the available space for the
%motion of the metabolites. The metabolites diffusion can be simulated using either
%a Crank-Nicholson approach (CN) or a crowding adaptation of the Lattice Boltzmann method (cLBM). 

%For this tutorial, we use as illustrative example the mutualistic
%consortium formed by E. coli AmetB and two subpopulations of S. enterica 
%(WT cells and methionine-secreting mutants). 
%Other examples are provided in the file Tutorial_moreExamples.m

%All simulated data reported in the paper are available at https://doi.org/10.5281/zenodo.5009045.



%% Setting up the paths 

clear
clc

%Add path to CROMICS in the file addPathCROMICS.m
addPathCROMICS 

%CROMICS can use neural networks (NN) or the genome-scale models (GEM) of the species
%to compute the metabolic fluxes. 
%If you want to use thermo-curated GEMs and TFA, install matTFA (https://github.com/EPFL-LCSB/matTFA)
%prior running CROMICS,and add the matTFA path in the file addPathCROMICS.m
%Alternatively, GEMs can be solved using flux balance analysis. For this, add the 
%COBRAtoolbox path in the file addPathCROMICS.m


%% Create Data and Results folders

%Dataset used in CROMICS simulations are stored in "Data" folder
%Check if Data/ exists, if not create it  
if ~exist('Data', 'dir')
    mkdir('Data')
end

%Results computed by CROMICS are stored in "Results" folder
%Check if Results/ exists, if not create it  
if ~exist('Results', 'dir')
    mkdir('Results')
end


%% Preparing the community model for CROMICS simulation 

%All system parameters required for the simulation (including method preferences 
%such as CN or cLBM, metabolites and species) should be defined in a 'ParameterDataFile.m' 
%located in:     ./Parameters/ParameterDataFile.m
%By running 'ParameterDataFile.m' a new 'ParameterDataFile.mat' will be created,
%wich will be used by CROMICS

%For example, the parameters used to simulate the microbial consortium
%formed by 50% of E. coli AmetB and 50% of S. enterica with metabolic variability
%are given in        ./Parameters/ParametersSystem_CROMICS_Ecoli50_Senterica50_crowAz.m
%By running the file:
ParametersSystem_CROMICS_Ecoli50_Senterica50_crowAz
%the new file 'dataForCROMICS_Ecoli50_Senterica50_30spots_crow0.4_Nmut70.mat' is created 
%and store in:   ./Data/dataForCROMICS_Ecoli50_Senterica50_30spots_crow0.4_Nmut70.mat

%If you already have a ParameterDataFile.mat you want to use, then skip this section,
%To run CROMICS you will need to provide only the name of the ParameterDataFile.mat 



%% Perform CROMICS simulations

%Provide the name of the ParameterDataFile.mat of the system you want to simulate, 
%For the experiment 50% of E. coli AmetB and 50% of S. enterica with
%phenotypic variability
nameDataFile='dataForCROMICS_Ecoli50_Senterica50_30spots_crow0.4_Nmut70.mat'

%Select the number of repetitions of the simulation.
repete=3; 

%Select the crowding assumption for each simulation repetition
%cA=1 when the radii of the cell is proportional to the cell mass 
%cA=4 when the cells are volumeless and Dmet=Dmet_water 
cA=1;

%Call CROMICS and parallelize the simulation
%Note that the parallelization requires Parallel Computing Toolbox
parfor i=1:repete 
    fprintf('repetition %g\n',i);
    
    %*******
    %If you want to randomly modify the initial positions of the subpopulations using bacterial spots given in the original file "nameDataFile" , if not skip this section
    %This section was used for the subpopulation experiment:
    %Ecoli50%_Senterica50%, where Senterica_mutant is 1% of the total Senterica population
     [nameDataFile_new]=randPositionSubpop(nameDataFile,i);
    %*******

    %Run CROMICS using the ParameterDataFile.mat:
    %When the radii of all metabolites are zero, i.e. inert macromolecules are
    %not present in the system, then we suggest to use the function 'CROMICSsq_launch.m'. 
    %In 'CROMICSsq_launch.m' the activity coefficient (given by the SPT equation) 
    %has been simplied to activity_coefficient = 1 -(space_occupied_by_cells). 
    feval(strcat('CROMICSsq_launch'),i,cA,nameDataFile_new);
    
    %If you want to simulate the diffusion of inert macromolecules (whose radii are
    %greater than 0), then use the function 'CROMICSsq_SPT_launch.m'
    %In 'CROMICSsq_SPT_launch.m' the activity coefficient is computed using the 
    %complete SPT equation that depends on the radii and abundance of cells 
    %and macromolecules. Thus, the simulation using 'CROMICSsq_SPT_launch.m' is 
    %computationally more expensive than in the simplified version 'CROMICSsq_launch.m'
    %described above). For thi case use the following line:
%     feval(strcat('CROMICSsq_SPT_launch'),i,cA,nameDataFile_new);
end

%Results are saved in:    ./Results/ResultsFile.mat
%using the name 'ResultsFile' that was provided by the user in the 'ParameterDataFile.m'

%Some examples for the visualization of results are given in the file
%./ResultsVisualization/ResultsVisualization.m

