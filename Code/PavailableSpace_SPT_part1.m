function [dens,dens_int,r,rint,vff]=PavailableSpace_SPT_part1(metaB,CellSizeVar,vol_lbm,v_sp,IDpos,frac_medium,Ncells,Reps1,Rcell_cte,Rmet,Rprot,metlbm,rhoLBM_met,AvN,fepsCell,num_eps,num_met,Ny_lbm,Nx_lbm,Nz_lbm,MMprot)
%Compute vectors dens [molecules mm^3] and r [mm] to estimate the SPT-probability 

if num_eps>0
    Reps=repmat(Reps1,Ncells,1);
    d_eps=(fepsCell(IDpos(1:Ncells,2),:).*IDpos(1:Ncells,3))./(vol_lbm*frac_medium(IDpos(1:Ncells,4))); %Density of EPS molecules in a LBM-voxel [molecules mm^-3].
    vff_eps=repmat(IDpos(1:Ncells,4),1,num_eps);
else
    Reps=[];
    d_eps=[];
    vff_eps=[];
end

if CellSizeVar==0
    r=[Reps(:);Rcell_cte(IDpos(1:Ncells,2));Rmet(metlbm)]; %Vector with the radii of all molecules in the system.
else
    Rcell=(3*(IDpos(1:Ncells,3)./metaB(IDpos(1:Ncells,2))).*v_sp(IDpos(1:Ncells,2))/(4*pi)).^(1/3);
    r=[Reps(:);Rcell;Rmet(metlbm)];
end

rint=[Reps(:);repmat(Rprot,Ncells,1);Rmet(metlbm)];
d_int=IDpos(1:Ncells,3).*(v_sp(IDpos(1:Ncells,2))>0)*(AvN/MMprot)./(vol_lbm*frac_medium(IDpos(1:Ncells,4))); %Density of intracellular macromolecules in a LBM-voxel [molecules mm^-3].

d_cell=metaB(IDpos(1:Ncells,2)).*(ones(Ncells,1)./(vol_lbm*frac_medium(IDpos(1:Ncells,4)))); %Density of cells in a LBM-voxel [cells mm^-3].
d_met=rhoLBM_met(:).*(Rmet(metlbm)>0)*AvN./(1000*vol_lbm*frac_medium(:)); %Density of metabolites molecules in a LBM-voxel [molecules mm^-3].
dens=[d_eps(:);d_cell;d_met];
dens_int=[d_eps(:);d_int;d_met];

vff=[vff_eps(:);IDpos(1:Ncells,4);repmat((1:Ny_lbm*Nx_lbm*Nz_lbm)',num_met,1)];



   
    