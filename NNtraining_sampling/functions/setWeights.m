function model = setWeights(model,group,myW,rescaleFlag)
%Original function of CAFBA
%Mori M, Hwa T, Martin OC, De Martino A, Marinari E. Constrained Allocation Flux Balance Analysis. PLoS Computational Biology. 2016; 12:e1004913.

% This function can be used to set the values of the weights of the
% different proteome sectors, instead of manually editing the 'w' field in
% the model.
%
% Inputs:
% 
%   model : a valid CAFBA model (including the protGroup field)
%   myW : value for the weights in the selected group (default: 0)
%   group: proteome group to be varied (integer number from 1 to 4)
%   rescaleFlag: this is useful in the case of randomly generated weights.
%          If 0 (default), the function sets all the weights in the group
%          to the value specified by myW. If 1, instead, the function
%          rescales all weights in the proteome sector proportionally,
%          so as to fix their average value to myW without affecting the
%          ratios of different weights, w_i/w_j for any i,j in the group.
%
% Output:
%
%   model: a valid CAFBA model (including the protGroup and w fields).
%
% See also the setWeights_rand() function.

% Check the inputs
if ~exist('model','var')
    error('No model was provided. Exiting from setWeights().');
end;
if ~isfield(model,'protGroup') % extra check on the model
    error('protGroup field missing. Exiting from setWeights().');
end;
if ~isfield(model,'w') % extra check on the model
    model.w=zeros(length(model.rxns),1);
end

% Defaults
if ~exist('myW','var') || isempty(myW)
    myW=0;
end;
if ~exist('group','var') || isempty(group)
    group=1; % C group
end;
if ~exist('rescaleFlag','var')
    rescaleFlag=0;
end




% All rxns in group
if group>0
    % Just the group
    if rescaleFlag==0
        % All weights are uniformly set to 'myW'
        for ri=1:length(model.protGroup(group).rxns)
            r=model.protGroup(group).rxns(ri);
            model.w(r)=myW;
        end;
    elseif rescaleFlag==1 
        % All weights in the group are rescaled such that
        % the new average is 'myW'
        average=0;
        for ri=1:length(model.protGroup(group).rxns)
            r=model.protGroup(group).rxns(ri);
            average=average+model.w(r);
        end;
        average=average/length(model.protGroup(group).rxns);
        for ri=1:length(model.protGroup(group).rxns)
            r=model.protGroup(group).rxns(ri);
            model.w(r)=(myW/average)*model.w(r);
        end;        
    else
        % All weights in the group are rescaled by a factor 'myW'
        for ri=1:length(model.protGroup(group).rxns)
            r=model.protGroup(group).rxns(ri);            
            model.w(r)=myW*model.w(r);
        end
    end
else
    error('group>length(model.protGroup) or group<0 in setWeights()');
end

end

