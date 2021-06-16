function [RxnF,RxnR]=FandBrxns(model,rxnNames)

%Identify the index of Forward and Reverse reactions of the fluxes of interest
%(required to perform TFA)
rxn=length(rxnNames);
RxnF=zeros(rxn,1);
RxnR=zeros(rxn,1);
for r=1:rxn
    RxnF(r)=find(ismember(model.varNames,strcat('F_',rxnNames(r))));
    RxnR(r)=find(ismember(model.varNames,strcat('R_',rxnNames(r))));
end

