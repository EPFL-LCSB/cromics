% Flux sampling for NN

clear

load Ecoli_core_lcts_TFA.mat

%Names of the exchange reactions whose metabolic flux will be computed
%biomass_rxn name should be the last one in the vector.
rxnNames={'DM_lcts_e';'DM_o2_e';'DM_met-L_e';'DM_ac_e';'DM_gal_e';'Ec_biomass_iJO1366_core_53p95M'};
rxnInterest=5;

%Range of fluxes in which samples will be taken. 
%No sample is taken if Vmin(rxn) = Vmax(rxn) = 0 for reaction rxn
Vmin=[0;0;0;0;0;0];  %[mmol gDW^-1 h^-1]
Vmax=[10;10;0.25;0;10;0];%[mmol gDW^-1 h^-1] %v3

%Number of samples to generate
nsamples=45000;

%Name of the file to save the input and output data
nameFile='TestPoints_Ecoli_core_lcts_pTFA_45K.mat'

%If there are more than 1 carbon source in the medium, but the cell cannot 
%uptake both substrate at the same time (catabolic repression)%same time, 
%then define here the unique Carbon source.
%If all the substrates can be consumed at the same time then
%uniqueCsource={};
uniqueCsource={'DM_lcts_e'};

%Select the type of sample distribution: 
% sampleDistribution='mixed';   %Use a mixed distributions: exponential (with mean 0.5, and a uniform distribution from Vmin to Vmax)
sampleDistribution='uniform'; %Use a uniform distributions (from Vmin to Vmax)



%***************



t1=cputime;
model=ttmodel;


%Take the absolute value for the uptake substrate.
Vmax=abs(Vmax); 
Vmin=abs(Vmin);
rxn=length(Vmin);



% Generating randomly the input data
% 2 "batches" : exp + const distribution
% Randomly generate either (with equal probability) :
% (1) a value from an exponential distribution with a mean of 0.5
% (2) a value from a uniform distribution between Vmin and Vmax.
if strfind(sampleDistribution,'mix')==1  
    numCs=length(uniqueCsource);
    numCs=(numCs==0)+numCs;
    nsamples2=round(nsamples/numCs);
    input_tab=[];
    for nC=1:numCs
        uCs=uniqueCsource;
        uCs(nC)=[];
        input=[];
        for r=1:rxn
            if Vmax(r)-Vmin(r)>0 && isempty(cell2mat(strfind(uCs,rxnNames(r))))==1 %Uptake flux
                prob = rand(nsamples2,1);
                a = (prob<=0.3) .* (Vmin(r) + exprnd(0.5,nsamples2,1)) + (prob>0.3) .* (Vmin(r)+((Vmax(r)-Vmin(r)).*rand(nsamples2,1)));
                a(a>Vmax(r)) = Vmax(r);
            elseif Vmax(r)==Vmin(r) && Vmin(r)>0 %Maximum uptake flux is fixed to single value
                a=ones(nsamples2,1)*Vmin(r);
            else %Biomass or metabolite production
                a=zeros(nsamples2,1);
            end
            input =[input,a];
        end
        input_tab=[input_tab;input];
    end
    
    
elseif strfind(sampleDistribution,'uni')==1 %We take a uniform distribution  
    numCs=length(uniqueCsource);
    numCs=(numCs==0)+numCs;
    nsamples2=round(nsamples/numCs);
    input_tab=[];
    for nC=1:numCs
        uCs=uniqueCsource;
        uCs(nC)=[];
        input=[];
        for r=1:rxn
            if Vmax(r)-Vmin(r)>0 && isempty(cell2mat(strfind(uCs,rxnNames(r))))==1 %Uptake flux
                a = (Vmin(r)+((Vmax(r)-Vmin(r)).*rand(nsamples2,1)));
                a(a>Vmax(r)) = Vmax(r);
            elseif Vmax(r)==Vmin(r) && Vmin(r)>0 %Maximum uptake flux is fixed to single value
                a=ones(nsamples2,1)*Vmin(r);
            else %Biomass or metabolite production
                a=zeros(nsamples2,1);
            end
            input =[input,a];
        end
        input_tab=[input_tab;input];
    end
   
end


%Identify the index of Forward and Reverse reactions of the fluxes we are
%interested in
RxnF=zeros(rxn,1);
RxnR=zeros(rxn,1);
for r=1:rxn
    if isempty(cell2mat(strfind(rxnNames(r),'bio')))==0
        bio = find(ismember(model.varNames,strcat('F_',rxnNames(r))));
        brxn = r;
    end
    
    RxnF(r)=find(ismember(model.varNames,strcat('F_',rxnNames(r))));
    RxnR(r)=find(ismember(model.varNames,strcat('R_',rxnNames(r))));
end


% Output variable
target_solTable = zeros(nsamples,rxn);


model.f(:)=0;
model.f(bio)=1;
model_orig=model;


tini=tic;

parfor ii=1:nsamples
  
    model=model_orig;

    % Modify the upper limits for exchange fluxes
    model.var_ub(RxnR)=input_tab(ii,:)'; 

    % Max biomass 
    [flux]=solveProblem(model,rxn,RxnF,RxnR);

    % We do the second optimization only if the biomass > 0
    model.var_lb(bio)=flux(brxn)*.99; 
    model.var_ub(bio)=flux(brxn); 

    % pFBA
    model.f(1:2*size(model.rxns,1))=-1;  %Minimize the sum of fluxes 
    [flux]=solveProblem(model,rxn,RxnF,RxnR);
        
    target_solTable(ii,:)=flux;

end

time=toc(tini);
input_tab = -input_tab; % We prefer negative inputs for our data

% Saving the results in a data file
save(nameFile, 'input_tab', 'target_solTable', 'nsamples','time','rxnNames');



function [flux]=solveProblem(model,rxn,RxnF,RxnR)
    sol=solveTFAmodelCplex(model);
    flux=zeros(1,rxn);
    
    if isempty(sol.val)==0 && isnan(sol.val)==0
        flux(:)=sol.x(RxnF)-sol.x(RxnR); 
    end
end