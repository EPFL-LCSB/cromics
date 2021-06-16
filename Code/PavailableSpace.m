function [frac,VavLBM]=PavailableSpace(CellSizeVar,Vcell_cte,vol_lbm,v_sp,fwater,IDpos,Ny_lbm,Nx_lbm,Nz_lbm,num_met,frac_medium,Ncells,frac)
%Compute the probability to find available space in each LBM-voxel.


if CellSizeVar==0 %Cell size is constant
    ff=Vcell_cte(IDpos(1:Ncells,2))/vol_lbm;
else %Cell size is proprtional to the cell biomass
    ff=IDpos(1:Ncells,3).*v_sp(IDpos(1:Ncells,2))/vol_lbm;
end

%Compute the space fraction occupied by cells
frac1=frac(:,1);
frac1(1:max(IDpos(1:Ncells,4)))=accumarray(IDpos(1:Ncells,4),ff);
frac(:)=repmat(frac1,num_met,1);
frac(Nx_lbm*Ny_lbm*Nz_lbm+1:Nx_lbm*Ny_lbm*Nz_lbm*2)=frac1*(1-fwater);

%Compute the space fraction available for the metabolites
frac(:)=frac_medium(:)-frac(:); %Matrix with the probability to find available space in the LBM-grid [dimensionless].
frac(:)=(frac(:)>0).*frac(:);

VavLBM=(frac==0)+frac*vol_lbm; %Available volume in a LBM-voxel [mm^3].
