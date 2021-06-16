function [RxnF]=Frxns(model,rxnNames)

%Identify the index of the reactions of interest
%(required to perform FBA)
rxn=length(rxnNames);
RxnF=zeros(rxn,1);
for r=1:rxn
    RxnF(r)=find(ismember(model.rxns,rxnNames(r)));
end

