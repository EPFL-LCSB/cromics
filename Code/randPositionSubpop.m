function [nameDataFile_new]=randPositionSubpop(nameDataFile,rep)

%Here, we randomly modify the initial positions of the subpopulations using bacterial spots given in the original file "nameDataFile" 
%rN is the vector that identify if a cell is a mutant (rN=1) or wild type (rN=0) 

load (nameDataFile)

[IDpos,Ncells]=CellProperties(posCA_sp,Ny_ca,Nx_ca,Nz_ca,num_sp,rhoCA_sp,ageCA,Q,IDca);

%E. coli, all cells are methionine auxotroph mutants
sp=1;
ff=find(rhoCA_sp(:,:,:,sp));
rN1=ones(length(ff),1); %rN1=1 for mutants 

%S. enterica (mutants are only in 1 spots)
sp=2;
vv=bacSpot_pos;
ff=find(rhoCA_sp(:,:,:,sp));
rN2=zeros(length(ff),1);
xx_pos=zeros(length(ff),1);
xx_pos(:)=vv(ff);
voxMutant=randperm(bacterialspots,length(Nmutant));

for vv=1:length(Nmutant)
    f1=find(xx_pos==voxMutant(vv));
        %Mutants are randomly distributed in the bacterial spot
        gg=f1(randperm(length(f1)));
        rN2(gg(1:Nmutant(vv)))=1;
end

rN=[rN2;rN1];

%Save the new file nameDataFile_new
nameDataFile_new=strcat('./Data/',strrep(nameDataFile,'.mat',strcat('rep',num2str(rep),'.mat')));
save (nameDataFile_new);

