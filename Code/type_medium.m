function [indCsat,rD_ind]=type_medium(mediumLBM,num_met,metlbm,rD)

indCsat=find(mediumLBM(:)==4);      %Index of the LBM-voxels that are occupied by air.
med=repmat(mediumLBM(:),num_met,1); %Medium type in each LBM-voxel [dimensionless].

%rD(metlbm,med) Ratio between the diffusion coefficient at medium1 and
%water, i.e. rD=Dsp(medium1)/Dsp(water):
rDind=sub2ind(size(rD),metlbm,med); 
rD_ind=rD(rDind);      