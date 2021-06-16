function model = setUpModel_aTFA(model,subst)
%Function adapted from the original CAFBA
%Mori M, Hwa T, Martin OC, De Martino A, Marinari E. Constrained Allocation Flux Balance Analysis. PLoS Computational Biology. 2016; 12:e1004913.


% Add protein groups to the model. The groups are numbered in this way:
% 1 : C
% 2 : E
% 3 : R
% 4 : Q
model=addProteinGroupsToModel_aTFA(model,subst);

% For example, this is the code to model lactose:
%   model=addProteinGroupsToModel_aTFA(model,'lac');
% Look inside the function's source code for more informations.

% We now set the offsets for the different groups.
model.protGroup(1).phi0 = 0; %if this value is equal to zero then the substrate "subst" selected will not affect the simulation. 
model.protGroup(2).phi0 = 0;
model.protGroup(3).phi0 = 0.066; %ro=rmin=0.0661 micro_g_RP/micro_g_TP
model.protGroup(4).phi0 = 0.45;

% Now we have to set the weights for the different groups.
% (note that the function creates the field model.weights if needed)
model=setWeights(model,1,0);
model=setWeights(model,2,0.00083);%1.55*10^-3);%0.00083);
% model=setWeights(model,2,1.55*10^-3);%0.00083);
model=setWeights(model,3,0.169);
model=setWeights(model,4,0);

% SETTING UP THE ALLOCATION CONSTRAINT
model.A=[model.A ; model.w'];

% right hand side
mySum=0;
for n=1:length(model.protGroup)
    mySum=mySum+model.protGroup(n).phi0;
end
model.rhs=[model.rhs;1-mySum];
model.constraintNames=[model.constraintNames;'allocation'];
model.constraintType=[model.constraintType;'='];

