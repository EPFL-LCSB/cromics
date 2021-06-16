
function [nameDataFile_new]=randPositionBacSpots(nameDataFile,rep)

%Here, we randomly modify the initial positions of the bacterial spots from the original file "nameDataFile" 

load (nameDataFile)
rhoCA_sp(:)=0;       %Matrix-CA of the density distribution for cell species [g].

%Medium data. Voxels can contain different medium 1=solids (porous medium), 2=polysaccharide matrix, 3=agar, 4=air. 
medium=ones(Ny_ca,Nx_ca,Nz_ca)*4;   %All voxels have "air" by default.
npore=ceil(Nx_ca*Ny_ca*Nz_ca*(1-porosity));
xyz=randperm(Ny_ca*Nx_ca*Nz_ca);
medium(xyz(1:npore))=1;      %Voxel is occupied by solids.

nN=20;
nN1=nN^2;           %Square bacterial spot of 20 by 20 CA-voxels as used in the experiment E.coli1_S.enterica99 
num=[ceil((nN1)*((100-Ecoli)/100)),(nN1)-ceil((nN1)*((100-Ecoli)/100))]; %Number of cells of each species located in one bacterial spot
Nxx=Nx_ca/nN;       %Assume a lattice division to allocate the bacterial spots. 

aa=randperm(Nxx*Nxx)';
ii=1;
bac=0;
while bac<bacterialspots
    [i,j]=ind2sub([Nxx,Nxx],aa(ii));
    
    if  isempty(find(rhoCA_sp((i*nN)-(nN-1):i*nN,(j*nN)-(nN-1):j*nN,:,:)))
        %Here, we randomly allocate the cells of E. coli and S.enterica in a
        %bacterial spot
        rhoProv=zeros(nN,nN,1,num_sp);
        mu=massB*(Ecoli/100)/num(2);
        rhoProv(:,:,1,1)=normrnd(mu,mu*s_m,[nN,nN]);
        rr=randperm(nN*nN)';
        rhoProv(rr(1:num(1)))=0;
        sig=sum(rhoProv(:))-massB*(Ecoli/100);
        check=rhoProv-sig/num(2);
        rhoProv=(rhoProv>0).*(check>0).*check;
        mu=massB*((100-Ecoli)/100)/num(1);
        rhoProv(rr(1:num(1))+(nN*nN))=normrnd(mu,mu*s_m,[num(1),1]);    
        sig=sum(sum(rhoProv(:,:,1,2)))-massB*((100-Ecoli)/100);
        check=rhoProv(:,:,1,2)-sig/num(1);
        rhoProv(:,:,1,2)=(rhoProv(:,:,1,2)>0).*(check>0).*check;
    
        rhoCA_sp((i*nN)-(nN-1):i*nN,(j*nN)-(nN-1):j*nN,:,:)=rhoProv;
        medium((i*nN)-(nN-1):i*nN,(j*nN)-(nN-1):j*nN,:)=2;
        bac=bac+1;
    end
    ii=ii+1;
end

posCA_sp(:)=0;              %Cell position matrix-CA of all the species.
socc=find(medium==1 | medium==3);               %Voxels occupied by inert solids (porous medium), i.e. medium==1, or by agar, i.e. medium==3, and therefore these voxels cannont contain cells.
posCA_sp(socc)=num_sp+1;                        %Voxel is occupied by inert solids (porous medium).

%Square grid
[Q,IDlbm,IDca,FD,FBD,vol_ca,vol_lbm,posCA_sp,frac_medium,frac_cells,frac_intcells,a1y,c1y,a1x,c1x,a1z,c1z,b1y,b2y,b1x,b2x,b1z,b2z,ap1x,cp1x,ap1y,cp1y,ap1z,cp1z,CNx,CNy,CNz]=SqLatticeCNandLBM(Nx_lbm,Ny_lbm,Nz_lbm,num_met,num_sp,Ax_ca,Ay_ca,Az_ca,Ax_lbm,Ay_lbm,Az_lbm,nx,ny,nz,Nx_ca,Ny_ca,Nz_ca,dim,rhoCA_sp,posCA_sp,boundaryLBM,boundaryCA,CmetDb,v_sp,v_met,medium);

%Compute mediumLBM, i.e. the medium (2=EPS, 3=agar, 4=air) contained ina each voxel-LBM  
mediumLBM=ones(Ny_lbm,Nx_lbm,Nz_lbm)*2;   %All voxels have "air" by default.
for med=3:-1:2
    aa=find(medium==med);
    mediumLBM(IDca(aa,Q+1))=med;
end

%Compute rhoLBM_met
rhoLBM_met(:)=0; %Matrix-LBM of the density distribution for met [g].
for met=1:num_met
    if met==2
        rhoLBM_met(Nx_lbm*Ny_lbm*Nz_lbm*(met-1)+1:Nx_lbm*Ny_lbm*Nz_lbm*met)=(mediumLBM(:)==2).*Cmet(met)*vol_lbm.*(frac_medium(Nx_lbm*Ny_lbm*Nz_lbm*(met-1)+1:Nx_lbm*Ny_lbm*Nz_lbm*met)'-frac_intcells(:));
    else
        rhoLBM_met(Nx_lbm*Ny_lbm*Nz_lbm*(met-1)+1:Nx_lbm*Ny_lbm*Nz_lbm*met)=(mediumLBM(:)==2).*Cmet(met)*vol_lbm.*(frac_medium(Nx_lbm*Ny_lbm*Nz_lbm*(met-1)+1:Nx_lbm*Ny_lbm*Nz_lbm*met)'-frac_cells(:));
    end
end

ageCA=rhoCA_sp>0;

var={'nCell','medium','npore','xyz','frac_cells','frac_intcells'};
clear(var{:});

%Save the new file nameDataFile_new
nameDataFile_new=strcat('./Data/',strrep(nameDataFile,'.mat',strcat('rep',num2str(rep),'.mat')));
save (nameDataFile_new);

