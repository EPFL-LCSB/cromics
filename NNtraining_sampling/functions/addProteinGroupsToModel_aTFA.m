function model=addProteinGroupsToModel_aTFA(model,optC)
%Function adapted from the original CAFBA
%Mori M, Hwa T, Martin OC, De Martino A, Marinari E. Constrained Allocation Flux Balance Analysis. PLoS Computational Biology. 2016; 12:e1004913.

% This function adds the 'protGroup' field to a COBRA model.
% Inputs:
%    model : COBRA metabolic model
%    optC  : Chooses the content of the phiC group. Here we report some
%            examples:
%            {'default'} (equivalent to 'glc')
%            'glc'
%            'lcts'
%            'fru'
%            'akg'
%            'g6p'
%            'g6p+glcn'
%           The full list can be found in the setPhiC() subfunction.
% Output:
%    model : same model as input, with the additional 'protGroup' field
%
%

if ~exist('optC','var')
    optC='glc';
end
phiEopt='subSys';

% Reactions
Biomass=find(model.f);

if ~isfield(model,'protGroup')
    protGroup(1).name='phiC';
    protGroup(1).phi0=0;
    [model,rxns_C] = setPhiC(model,optC);
    protGroup(1).rxns=rxns_C;
    
    % Defines all protGroup entries
    protGroup(2).name='phiE';
    protGroup(2).phi0=0.45;
    protGroup(2).rxns=set_PhiE_list(model,phiEopt);
    protGroup(2).rxns=removeRxnsOverlap(protGroup(2).rxns,protGroup(1).rxns);
    
    protGroup(3).name='phiR';
    protGroup(3).phi0=0.066;
    protGroup(3).rxns=Biomass;
    protGroup(3).rxns=removeRxnsOverlap(protGroup(3).rxns,protGroup(1).rxns);
    protGroup(3).rxns=removeRxnsOverlap(protGroup(3).rxns,protGroup(2).rxns);
    
    protGroup(4).name='phiQ and others'; % weights=0
    protGroup(4).phi0=0;
    protGroup(4).rxns=(1:length(model.varNames))';
    protGroup(4).rxns=removeRxnsOverlap(protGroup(4).rxns,protGroup(1).rxns);
    protGroup(4).rxns=removeRxnsOverlap(protGroup(4).rxns,protGroup(2).rxns);
    protGroup(4).rxns=removeRxnsOverlap(protGroup(4).rxns,protGroup(3).rxns);
    
    model.protGroup=protGroup;
else
    % only updates the first prot group
    [model,rxns_C] = setPhiC(model,optC);
    model.protGroup(1).rxns=rxns_C;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION - Set the PhiC rxns vector and bounds %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model,rxns_C] = setPhiC(model,optC)
% This function also updates the bound for the exchange reactions

if strcmp(optC,'glc')
    rxns_C=find(strcmp(model.varNames,strcat('R_','DM_glc_e'))); % glc transport
elseif strcmp(optC,'lcts')
    rxns_C=find(strcmp(model.varNames,strcat('R_','DM_lcts_e'))); % lcts transport
elseif strcmp(optC,'ac')
    rxns_C=find(strcmp(model.varNames,strcat('R_','DM_ac_e'))); % lcts transport
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION - Set the PhiE rxns vector  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phiEVec=set_PhiE_list(model,phiEopt)
nRxns=length(model.rxns);

phiEVec=[];
subSysList={'Transport Extracellular'; ...
            'Transport, Inner Membrane'; ...
            'Transport, Outer Membrane'; ...
            'Transport, Outer Membrane Porin'; ...
            'Putative Transporters'; ...
            'Exchange'; ...
            'Others'; ...
            '';...
            'Transport Inner Membrane'; ...
            'Transport Outer Membrane'; ...
            'Transport Outer Membrane Porin'; ...
            'Alternate Carbon Metabolism';... 
            };
%             'Alternate Carbon Metabolism';...  %LAM: I have eliminated 

p1=find(~ismember(model.subSystems,subSysList));
% %*** LAM:
% aa=find(ismember(model.rxns,'ATPM'));
% p1(find(p1==aa))=[];
% Biomass=find(model.f);
% p1(find(p1==Biomass))=[];
% %***
phiEVec=find(ismember(model.varNames,[strcat('R_',model.rxns(p1));strcat('F_',model.rxns(p1))]));

end

function vec3=removeRxnsOverlap(vec1,vec2)
% Subtracts from vec1 all entries  which also appear in vec2
vec3=vec1(~ismember(vec1,vec2));
end

