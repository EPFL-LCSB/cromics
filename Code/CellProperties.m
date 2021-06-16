function [IDpos,Ncells]=CellProperties(posCA_sp,Ny_ca,Nx_ca,Nz_ca,num_sp,rhoCA_sp,ageCA,Q,IDca)
%Create an identification matrix (ID) with properties and position of each cell.

ijk=find(posCA_sp<=num_sp & posCA_sp>0);
Ncells=length(ijk);
IDpos=zeros(Ny_ca*Nx_ca*Nz_ca,5);


if Ncells>0
    %IDpos=[posCA,sp,X,posLBM,age,abio]   Matrix containig the info of the microbial cell 
    IDpos(1:Ncells,2)=posCA_sp(posCA_sp<=num_sp & posCA_sp>0); %Microbial species type.
    IDpos(1:Ncells,1)=Ny_ca*Nx_ca*Nz_ca*(IDpos(1:Ncells,2)-1)+ijk; %Index position of the cell in the CA-grid.
    IDpos(1:Ncells,3)=rhoCA_sp(IDpos(1:Ncells,1)); %Biomass of each cell [gDW].
    IDpos(1:Ncells,4)=IDca(ijk,Q+1); %Corresponding index position of the cell in the LBM-grid.
    IDpos(1:Ncells,5)=ageCA(IDpos(1:Ncells,1));
end

[~,hh]=sort(IDpos(:,2),'descend');
IDpos=IDpos(hh,:);
