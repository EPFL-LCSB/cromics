%% CROwding-Modeling of In-silico Community Systems (CROMICS)

% Paper:
% 1. Angeles-Martinez and Hatzimanikatis. "The influence of the crowding assumptions in biofilm simulations".

%This tutorial demonstrate the use of CROMICS for the spatio-temporal simulation 
%of biofilms.

%CROMICS is an individual-based model wherein the effect of the crowding conditions 
%is explicitly incorporated. CROMICS combines methods and techniques such as 
%stoichiometric-based models (in particular thermodynamic flux analysis (TFA) [Henry et al., 2007])
%to compute the metabolic fluxes of each individual cell in the system, and
%scaled particle theory (SPT) [Lebowitz et al., 1965] to estimate the available space for the
%motion of the metabolites. The metabolites diffusion can be simulated using either
%a Crank-Nicholson approach (CN) or a crowding adaptation of the Lattice Boltzmann method (cLBM). 

%All simulated data reported in the paper are available at https://doi.org/10.5281/zenodo.5005946.



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



%% Experiment: Active layer depth. System is filled with E. coli cells, constant glucose supply
    
clear
clc

%Files below contain all the parameters requiered for the simulation
%           ./Parameters/ParametersSystem_CROMICSsq_Ecoli_biofilm_cteSupplyGlc.m
%By running these files, a ParameterDataFile.mat is created and stored in:
%           ./Data/dataForCROMICSsq_Ecoli_biofilm_cteSupply2.25.mat
% which will be used as follows:

%Provide the name of the ParameterDataFile.mat, for example:
nameDataFile='dataForCROMICSsq_Ecoli_biofilm_cteSupply2.25.mat'

%Select the number of repetitions of the simulation.
repete=3;       %Number of repetitions of the simulation.
  
%Select the crowding assumption for each simulation repetition
%cA=1 when the radii of the cell is proportional to the cell mass 
%cA=2 when the radii of the cell is constant 
%cA=3 when the cells are volumeless and Dmet=0.25Dmet_water 
%cA=4 when the cells are volumeless and Dmet=Dmet_water 
cA=1; 

%Call CROMICS and parallelize the simulation
%Note that the parallelization requires Parallel Computing Toolbox
parfor i=1:repete 
    fprintf('repetition %g\n',i);
    
    %Run CROMICS using the ParameterDataFile.mat:
    %When the radii of all metabolites are zero, i.e. inert macromolecules are
    %not present in the system, then we suggest to use the function 'CROMICSsq_hex_launch.m'. 
    %In 'CROMICSsq_launch.m' the activity coefficient (given by the SPT equation) 
    %has been simplied to activity_coefficient = 1 -(space_occupied_by_cells). 
    feval(strcat('CROMICSsq_hex_launch'),i,cA,nameDataFile);
    
end

%Results are saved in:    ./Results/ResultsFile.mat
%using the name 'ResultsFile' that was provided by the user in the 'ParameterDataFile.m'

%Some examples for the visualization of results are given in the file:
%./ResultsVisualization/ResultsVisualization.m



%% Experiment: Ecoli biofilm using a square lattice
    
clear
clc

%Files below contain all the parameters requiered for the simulation
%           ./Parameters/ParametersSystem_CROMICSsq_Ecoli_biofilm.m
%By running these files, a ParameterDataFile.mat is created, which will be
%used as follows:

%Provide the name of the ParameterDataFile.mat, for example:
nameDataFile='dataForCROMICSsq_Ecoli_biofilm.mat'

%Select the number of repetitions of the simulation.
repete=3;       %Number of repetitions of the simulation.
  
%Select the crowding assumption for each simulation repetition
%cA=1 when the radii of the cell is proportional to the cell mass 
%cA=2 when the radii of the cell is constant 
%cA=3 when the cells are volumeless and Dmet=0.25Dmet_water 
%cA=4 when the cells are volumeless and Dmet=Dmet_water 
cA=1; 

%Call CROMICS and parallelize the simulation
%Note that the parallelization requires Parallel Computing Toolbox
parfor i=1:repete 
    fprintf('repetition %g\n',i);
    
    %Run CROMICS using the ParameterDataFile.mat:
    %When the radii of all metabolites are zero, i.e. inert macromolecules are
    %not present in the system, then we suggest to use the function 'CROMICSsq_hex_launch.m'. 
    %In 'CROMICSsq_launch.m' the activity coefficient (given by the SPT equation) 
    %has been simplied to activity_coefficient = 1 -(space_occupied_by_cells). 
    feval(strcat('CROMICSsq_hex_launch'),i,cA,nameDataFile);
    
end

%Results are saved in:    ./Results/ResultsFile.mat
%using the name 'ResultsFile' that was provided by the user in the 'ParameterDataFile.m'

%Some examples for the visualization of results are given in the file:
%./ResultsVisualization/ResultsVisualization.m



%% Experiment: Ecoli biofilm using a hexagonal lattice
    
clear
clc

%Files below contain all the parameters requiered for the simulation
%           ./Parameters/ParametersSystem_CROMICShex_Ecoli_biofilm.m
%By running these files, a ParameterDataFile.mat is created, which will be
%used as follows:

%Provide the name of the ParameterDataFile.mat, for example:
nameDataFile='dataForCROMICShex_Ecoli_biofilm.mat'

%Select the number of repetitions of the simulation.
repete=3;       %Number of repetitions of the simulation.
  
%Select the crowding assumption for each simulation repetition
%cA=1 when the radii of the cell is proportional to the cell mass 
%cA=2 when the radii of the cell is constant 
%cA=3 when the cells are volumeless and Dmet=0.25Dmet_water 
%cA=4 when the cells are volumeless and Dmet=Dmet_water 
cA=1; 

%Call CROMICS and parallelize the simulation
%Note that the parallelization requires Parallel Computing Toolbox
parfor i=1:repete 
    fprintf('repetition %g\n',i);
    
    %Run CROMICS using the ParameterDataFile.mat:
    %When the radii of all metabolites are zero, i.e. inert macromolecules are
    %not present in the system, then we suggest to use the function 'CROMICSsq_hex_launch.m'. 
    %In 'CROMICSsq_launch.m' the activity coefficient (given by the SPT equation) 
    %has been simplied to activity_coefficient = 1 -(space_occupied_by_cells). 
    feval(strcat('CROMICSsq_hex_launch'),i,cA,nameDataFile);
    
end

%Results are saved in:    ./Results/ResultsFile.mat
%using the name 'ResultsFile' that was provided by the user in the 'ParameterDataFile.m'

%Some examples for the visualization of results are given in the file:
%./ResultsVisualization/ResultsVisualization.m



%% Experiment: Ecoli in planktonic state using a square lattice
    
clear
clc

%Files below contain all the parameters requiered for the simulation
%           ./Parameters/ParametersSystem_CROMICSsq_Ecoli_planktonic.m
%By running these files, a ParameterDataFile.mat is created, which will be
%used as follows:

%Provide the name of the ParameterDataFile.mat, for example:
nameDataFile='dataForCROMICSsq_Ecoli_planktonic.mat'

%Select the number of repetitions of the simulation.
repete=3;       %Number of repetitions of the simulation.
  
%Select the crowding assumption for each simulation repetition
%cA=1 when the radii of the cell is proportional to the cell mass 
%cA=2 when the radii of the cell is constant 
%cA=3 when the cells are volumeless and Dmet=0.25Dmet_water 
%cA=4 when the cells are volumeless and Dmet=Dmet_water 
cA=1; 

%Call CROMICS and parallelize the simulation
%Note that the parallelization requires Parallel Computing Toolbox
parfor i=1:repete 
    fprintf('repetition %g\n',i);
    
    %Run CROMICS using the ParameterDataFile.mat:
    %When the radii of all metabolites are zero, i.e. inert macromolecules are
    %not present in the system, then we suggest to use the function 'CROMICSsq_hex_launch.m'. 
    %In 'CROMICSsq_launch.m' the activity coefficient (given by the SPT equation) 
    %has been simplied to activity_coefficient = 1 -(space_occupied_by_cells). 
    feval(strcat('CROMICSsq_hex_launch'),i,cA,nameDataFile);
    
end

%Results are saved in:    ./Results/ResultsFile.mat
%using the name 'ResultsFile' that was provided by the user in the 'ParameterDataFile.m'

%Some examples for the visualization of results are given in the file
%./ResultsVisualization/ResultsVisualization.m



