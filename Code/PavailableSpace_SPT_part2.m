function [frac,VavLBM]=PavailableSpace_SPT_part2(dens,dens_int,r,rint,Ny_lbm,Nx_lbm,Nz_lbm,num_met,Rmet,metlbm,P_kT,frac,vff,vol_lbm)
%Compute the probability to find available space in each LBM-voxel.

%SPT variables:
s0=((pi/6)*dens);
s1=((2*pi/6)*dens.*r);
s2=((4*pi/6)*dens.*(r.^2));
s3=((8*pi/6)*dens.*(r.^3));
S0=repmat(accumarray(vff,s0),num_met,1);
S1=repmat(accumarray(vff,s1),num_met,1);
S2=repmat(accumarray(vff,s2),num_met,1);
S3=repmat(accumarray(vff,s3),num_met,1);

%For O2:
% s0=((pi/6)*dens_int); %We assume that the cell mass [gDW] = mass of intracellular proteins, therefore dens_int=dens --> S0=S0_int 
s1=((2*pi/6)*dens_int.*rint);
s2=((4*pi/6)*dens_int.*(rint.^2));
s3=((8*pi/6)*dens_int.*(rint.^3));
% S0_int=accumarray(vff,s0);
S1_int=accumarray(vff,s1);
S2_int=accumarray(vff,s2);
S3_int=accumarray(vff,s3);
% S0(Nx_lbm*Ny_lbm*Nz_lbm+1:Nx_lbm*Ny_lbm*Nz_lbm*2)=S0_int;
S1(Nx_lbm*Ny_lbm*Nz_lbm+1:Nx_lbm*Ny_lbm*Nz_lbm*2)=S1_int;
S2(Nx_lbm*Ny_lbm*Nz_lbm+1:Nx_lbm*Ny_lbm*Nz_lbm*2)=S2_int;
S3(Nx_lbm*Ny_lbm*Nz_lbm+1:Nx_lbm*Ny_lbm*Nz_lbm*2)=S3_int;

ln_act=-log(1-S3)+(6*S2./(1-S3)).*Rmet(metlbm)+((12*S1./(1-S3))+(18*(S2.^2)./((1-S3).^2))).*(Rmet(metlbm).^2)+((8*S0./(1-S3))+(24*S1.*S2./((1-S3).^2))+(24*(S2.^3)./((1-S3).^3))).*(Rmet(metlbm).^3);
frac(:)=exp(-ln_act); %Matrix with the available space fraction in the LBM-grid [dimensionless].

VavLBM=frac*vol_lbm; %Available volume in a LBM-voxel [mm^3].

    
 
aa=find(S3(:)>=1);
frac(aa)=0;
VavLBM(aa)=0;
 