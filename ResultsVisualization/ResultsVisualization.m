%% Visualization of results. 

%This file contains some examples to plot the results for their analysis.
%You will need provide the name of the ResultsFile.mat
%If resultsFile.mat is located outside the CROMICS folder, then provide also the path


%% Plot the abundance of species, metabolites and phenotypes, as well as the species relative abundance 

clc;
clear; 
clf;

%Provide the name of ResultsFile.mat
resultsFile = 'resultsCROMICS_Ecoli1_Senterica99_cA1_Nrep1'; 
load (resultsFile)
%Once the resultsFile is loaded, 'plotFunction.m' can be used to visualize the
%results. Some example are given below.

%The figure will be divided in a 'mplot' by 'nplot' grid.
nplot=2; 
mplot=2;

%Plot the total biomass of all species [g].
fig=1; %The figure will be located in position 'fig' in a 'mplot' by 'nplot' grid.
plotBiomass(mplot,nplot,fig,num_sp,tvec,rhotCA_sp);

%Plots the relative abundance of all species [g/g].
fig=fig+1;
plotSpeciesRelAbundance(mplot,nplot,fig,num_sp,tvec,rhotCA_sp);

%Plot the phenotype abundance of species 'sp'  [individuals].
fig=fig+1;
sp=1; %Species of interest
plotPhenotypes(mplot,nplot,fig,sp,tvec,phentCA);

%Plot the total abundance of all metabolites [mmol] 
fig=fig+1;
plotMetabolites(mplot,nplot,fig,num_met,tvec,rhotLBM_met);



%% Plot snapshots of phenotypes, metabolites and subpopulations

clc;
clear; 
clf;

%Provide the name of resultsFile.mat
resultsFile = 'resultsCROMICS_Ecoli50_Senterica50_subpop_cA1_Nrep1.mat'; 
load (resultsFile)
%Once the resultsFile is loaded, plotFunction can be used to visualize the results.

%The figure will be divided in a 'mplot' by 'nplot' grid.
nplot=2; 
mplot=2;

%Plot a snapshot of the metabolites distribution at the point 'h' that was store in the results matrixes 
fig=1;  %The figure will be located in position 'fig' in a 'mplot' by 'nplot' grid.
met=3;  %Metabolite to plot 
h=70;   %The snapshot corresponds to point 'h' that was store in the results matrixes 
%Color used to plot the metabolite gradient in the snapshot.
%If cmap2 is not provided a default color will be used.
cmap2=ones(1000,3); 
cmap2(:,1)=(1:-1/1000:1/1000)';
cmap2(:,3)=(1:-1/1000:1/1000)';
plotSnapshotMet_2D(mplot,nplot,fig,rhotLBM_met,tvec,h,met,cmap2);

%Plot a snapshot of the phenotypes distribution of species 'sp' at the point 'h' that was store in the results matrixes 
fig=fig+1; %The figure will be located in position 'fig' in a 'mplot' by 'nplot' grid.
sp=1;      %Species to plot 
h=70;      %The snapshot corresponds to point 'h' that was store in the results matrixes 
plotSnapshotPhen_2D(mplot,nplot,fig,phentCA,tvec,h,sp);

%Plot a snapshot of the subpopulation distribution of species 'sp' at the point 'h' that was store in the results matrixes 
fig=fig+1; %The figure will be located in position 'fig' in a 'mplot' by 'nplot' grid.
sp=2;      %Species to plot 
h=70;     %The snapshot corresponds to point 'h' that was store in the results matrixes 
plotSnapshotSubpop_2D(mplot,nplot,fig,subpoptCA,tvec,h,sp);

%Plot a snapshot of the growth rate of species 'sp' at the point 'h' that was store in the results matrixes 
fig=fig+1; %The figure will be located in position 'fig' in a 'mplot' by 'nplot' grid.
sp=1;      %Species to plot 
h=70;     %The snapshot corresponds to point 'h' that was store in the results matrixes 
plotSnapshotGrowthRate_2D(mplot,nplot,fig,VbiotCA_sp,tvec,h,sp);



%% Plot the evolution of the relative abundance of populations and subpopulations

clc;
clear; 
clf;

%Provide the name of resultsFile.mat
resultsFile = 'resultsCROMICS_Ecoli50_Senterica50_subpop_cA1_Nrep1.mat'; 
load (resultsFile)
%Once the resultsFile is loaded, plotFunction can be used to visualize the results.

%Plot the evolution relative abundance of populations and subpopulations in the whole system 
figure (1)
nplot=2; %The figure will be located in position 'fig' in a 'mplot' by 'nplot' grid.
mplot=2;
fig=1;  
plotFracSubpop(mplot,nplot,fig,rhotCA_sp,subpoptCA,tvec);

%Plot the evolution relative abundance of populations and subpopulations in different regions.
%The 2D system is divided in 25 regions. 
figure(2)
plotFracSubpop25Regions_2D(rhotCA_sp,subpoptCA,tvec);



%% Plot a 3D biofilm and snapshot of one layer of the biofilm

clc;
clear; 
clf;

%Provide the name of resultsFile.mat
resultsFile = 'resultsCROMICS_Ecoli_multiB_Nrep1_cA1.mat'; 
load (resultsFile)
%Once the resultsFile is loaded, 'plotFunction.m' can be used to visualize the results.

%The figure will be divided in a 'mplot' by 'nplot' grid.
nplot=2; 
mplot=2;

%Plot a 3D biofilm composed by 2 species
fig=1;  %The figure will be located in position 'fig' in a 'mplot' by 'nplot' grid.
h=361; %The snapshot corresponds to point 'h' that was store in the results matrixes 
plotSnapshotCells_3D(mplot,nplot,fig,postCA_sp,tvec,h,num_sp);

%Plot one layer of the 3D biofilm composed by 2 species
fig=fig+1;
h=361; %The snapshot corresponds to point 'h' that was store in the results matrixes 
plotSnapshotCells_3D_oneLayer(mplot,nplot,fig,postCA_sp,tvec,h,num_sp);

%Plot the space fraction occupied by cells and macrolecules in one layer of the 3D biofilm composed by 2 species
fig=fig+1;
h=361; %The snapshot corresponds to point 'h' that was store in the results matrixes 
plotSnapshotSpaceOccupied_3D_oneLayer(mplot,nplot,fig,rhotCA_sp,tvec,v_sp,vol_ca,h);

%Plot phenotypes in one layer of the 3D biofilm composed by 2 species
fig=fig+1;
h=361; %The snapshot corresponds to point 'h' that was store in the results matrixes 
plotSnapshotPhenotypes_3D_oneLayer(mplot,nplot,fig,phentCA,tvec,h);



%% Plot the mean square displacement (MSD) of the metabolites 

clc;
clear; 
clf;

%Provide the name of resultsFile.mat
resultsFile = 'resultsCROMICS_Ecoli_multiB_Nrep1_cA1.mat';
load (resultsFile)
%Once the resultsFile is loaded, 'plotFunction.m' can be used to visualize the results.

%The figure will be divided in a 'mplot' by 'nplot' grid.
nplot=2; 
mplot=2;
fig=1;  %The figure will be located in position 'fig' in a 'mplot' by 'nplot' grid.
plotMSDmetabolites(mplot,nplot,fig,r2tt_met,num_met,tvec);




