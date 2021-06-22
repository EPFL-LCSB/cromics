%% More Experiments to simulate using CROMICS

% Paper:
% 1. Angeles-Martinez and Hatzimanikatis. "Spatio-temporal modeling of the crowding conditions and metabolic variability in microbial communities".

%All simulated data reported in the paper are available at https://doi.org/10.5281/zenodo.5009045.



%% Add path to CROMICS in the file addPathCROMICS.m
addPathCROMICS 


%% Experiment: Ecoli99%_Senterica1% and Ecoli99%_Senterica1%
    
clear
clc

%Files below contain all the parameters requiered for the simulation
%           ./Parameters/ParametersSystem_CROMICS_Ecoli1_Senterica99.m 
%           ./Parameters/ParametersSystem_CROMICS_Ecoli99_Senterica1.m
%By running these files, a ParameterDataFile.mat is created, which will be
%used as follows:

%Provide the name of the ParameterDataFile.mat, for example:
%  nameDataFile='dataForCROMICS_Ecoli99_Senterica1_crow0.2.mat'
nameDataFile='dataForCROMICS_Ecoli1_Senterica99_crow0.2.mat'
%Select the number of repetitions of the simulation.
repete=3;       %Number of repetitions of the simulation.
  
%Select the crowding assumption for each simulation repetition
%cA=1 when the radii of the cell is proportional to the cell mass 
%cA=4 when the cells are volumeless and Dmet=Dmet_water 
cA=1; 

%Call CROMICS and parallelize the simulation
%Note that the parallelization requires Parallel Computing Toolbox
parfor i=1:repete 
    fprintf('repetition %g\n',i);
    
    %*******
    %If you want to randomly modify the initial positions of the bacterial spots from the original file "nameDataFile", if not skip this section
    %This section was used for the COMETS experiment: Ecoli99%_Senterica1% and Ecoli99%_Senterica1%
    [nameDataFile_new]=randPositionBacSpots(nameDataFile,i);
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

%Some examples for the visualization of results are given in the file:
%./ResultsVisualization/ResultsVisualization.m


%% Experiment: Methionine abundance and Ecoli growth at different distances from Senterica    

clear
clc

%Files below contain all the parameters requiered for the simulation
%           ./Parameters/ParametersSystem_CROMICS_Senterica_center_crowAz.m 
%By running these files, a ParameterDataFile.mat is created, which will be
%used as follows:

%Provide the name of the ParameterDataFile.mat, for example:
%For active cells:
%  nameDataFile='dataForCROMICS_Senterica_center_Vocc0.2_Vmax1_meth.mat'
%For inactive cells:
nameDataFile='dataForCROMICS_Senterica_center_Vocc0.2_Vmax0_meth.mat'
%Select the number of repetitions of the simulation.
repete=3;       %Number of repetitions of the simulation.
  
%Select the crowding assumption for each simulation repetition
%cA=1 when the radii of the cell is proportional to the cell mass 
%cA=4 when the cells are volumeless and Dmet=Dmet_water 
cA=1; 

%Call CROMICS and parallelize the simulation
%Note that the parallelization requires Parallel Computing Toolbox
parfor i=1:repete 
    fprintf('repetition %g\n',i);
    
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

%Some examples for the visualization of results are given in the file:
%./ResultsVisualization/ResultsVisualization.m




%% Experiment: Biofilm EPS-secreting Ecoli and non-EPS-secreting Ecoli 


clear
clc

%Files below contain all the parameters requiered for the simulation
%           ./Parameters/ParametersSystem_CROMICS_SPT_Ecoli_Biofilm_10eps.m 
%           ./Parameters/ParametersSystem_CROMICS_SPT_Ecoli_Biofilm_30eps.m
%By running these files, a ParameterDataFile.mat is created, which will be
%used as follows:

%Provide the name of the ParameterDataFile, here are 2 example:
% nameDataFile='dataForCROMICS_Ecoli_10eps.mat';
nameDataFile='dataForCROMICS_Ecoli_30eps.mat';

%Select the number of repetitions of the simulation.
repete=3;       %Number of repetitions of the simulation

%Select the crowding assumption for each simulation repetition
%cA=1 when the radii of the cell is proportional to the cell mass 
%cA=4 when the cells are volumeless and Dmet=Dmet_water 
cA=1;

t1=cputime;

%Call CROMICS and parallelize the simulation
%Note that the parallelization requires Parallel Computing Toolbox
parfor i=1:repete
    fprintf('repetition %g\n',i);
     %Run CROMICS using the ParameterDataFile.mat:
    %When the radii of all metabolites are zero, i.e. inert macromolecules are
    %not present in the system, then we suggest to use the function 'CROMICSsq_launch.m'. 
    %In 'CROMICSsq_launch.m' the activity coefficient (given by the SPT equation) 
    %has been simplied to activity_coefficient = 1 -(space_occupied_by_cells). 
    feval(strcat('CROMICSsq_launch'),i,cA,nameDataFile);
    
    %If you want to simulate the diffusion of inert macromolecules (whose radii are
    %greater than 0), then use the function 'CROMICSsq_SPT_launch.m'
    %In 'CROMICSsq_SPT_launch.m' the activity coefficient is computed using the 
    %complete SPT equation that depends on the radii and abundance of cells 
    %and macromolecules. Thus, the simulation using 'CROMICSsq_SPT_launch.m' is 
    %computationally more expensive than in the simplified version 'CROMICSsq_launch.m'
    %described above). For thi case use the following line:
%     feval(strcat('CROMICSsq_SPT_launch'),i,cA,nameDataFile);

end

%Results are saved in:    ./Results/ResultsFile.mat
%using the name 'ResultsFile' that was provided by the user in the 'ParameterDataFile.m'

%Some examples for the visualization of results are given in the file
%./ResultsVisualization/ResultsVisualization.m

