function [IDpos,Ncells_new,posCA_sp,r2t_sp,mediumLBM]=IbMcells_hex(mediumLBM,r2t_sp,Ncells_new,IDpos,Ny_ca,Nx_ca,Nz_ca,num_sp,posCA_sp,IDca,rV,Q,a1,a2,Mmaxsp,Mminsp,Macc1,Ax_ca)

       
%Motility and division of the cells
r2t=zeros(num_sp,1); %Vector of the mean squared displacement of the species 
Macc=zeros(Q-1,1); %Vector of the number of voxels occupied by cells in the vertical, horizontal, and diagonal directions around a daughter cell 
aa=a1(randperm(length(a1))); %Permute the order of the cells that will divide, move or die

for ii=1:a2
    rm=aa(ii);   %Randomly select a cell/colony
    sp=IDpos(rm,2);
    ijk=IDpos(rm,1)-Ny_ca*Nx_ca*Nz_ca*(sp-1);
 
    
    %Cell division
    if IDpos(rm,3)>=Mmaxsp(sp)   %The cell achieves its maximum biomass, and it will divide.
 
        %Define some variables:
        recalcFBA=0; %If recalcFBA=1 we have to recalculate FBA to consider that the cell cannot grow due to the space limitations, e.g. the cells is surrounded by solid walls
        spaceCA_orig=posCA_sp;
        IDpos_orig=IDpos;
        Ncells_neworig=Ncells_new;

        data_tc=IDpos(rm,:);    %Store data of the target cell, i.e. the one we want to find a new place. 
        data_tc([3,5])=[IDpos(rm,3)/(rV(sp)+1),1]; 

        IDpos(rm,3)=IDpos(rm,3)-data_tc(3);%Update the biomass of the mother cell [gDW].
        IDpos(rm,5)=IDpos(rm,5)+1;    %Update age of the mother cell [generations].
        sp_tc=sp;
        
        %Find empty neighbour sites where to place the target cell
        occup=posCA_sp(IDca(ijk,2:Q)); %IDca containe the nearest neighbours
        dir=find(occup==0);

        if isempty(dir)==0  %There are empty neighbour voxels (or sites) where to place the target cell.
            rm1=unidrnd(length(dir),1,1);   %Select a random empty site.
            drm1=dir(rm1)+1;
            posCA_sp(IDca(ijk,drm1))=sp_tc; %Update info with the target cell data.
            mediumLBM(IDca(ijk,Q+drm1))=2;
            data_tc([1,4])=[IDca(ijk,drm1)+Ny_ca*Nx_ca*Nz_ca*(sp_tc-1), IDca(ijk,Q+drm1)];
            Ncells_new=Ncells_new+1;
            IDpos(Ncells_new,:)=data_tc;
        
            
        else %Shoving
        %There is not empty voxels, therefore the dividing cell has to push the other cells to create an empty site to allocate the daughter cell
        
            %Identify if there are solid walls around the cell
            dir=((occup<num_sp+1).*(IDca(ijk,2:Q)~=IDca(ijk,1))); %i.e. occup=num_sp+1 (inert solids), IDca(ijk,2:Q)==IDca(ijk,1) (bounceback boundary)
            if all(dir==0)==1 %|| any(storeijk>20)==1 %The cell is surrounded by inert solids or it has been chosen more than 20 times in 1 IbM step.
                recalcFBA=1; %There is not space where to locate the daughter/displaced cell, therefore the cell cannot be divided. 
                             %We assume that the cell is deformed and can fullfill the whole voxel.  
                posCA_sp=spaceCA_orig;
                IDpos=IDpos_orig;
                Ncells_new=Ncells_neworig;
            end  

            if recalcFBA==0 %There ir no wall around, so the mother cell can push other cells. 
                            %Look for a new (and for now occupied) site where to locate the daughter/displaced cell 
                [i,j,k]=ind2sub([Ny_ca,Nx_ca,Nz_ca],ijk);

                
                %The idea here is to identify the direction (in a vonNeumann neighbourhood) with the nearest empty site. 
                
                %Vector Macc store the number of CA-voxels occupied by cells that are found between the mother cells and the nearest empty voxel.
                %If no empty site are found, or a wall is found before then Macc takes a default big number  
                Macc(:)=(10^30);
                for idir=1:length(dir)
                    if dir(idir)>0
                        dd=Macc1(ijk,dir(idir),:);
                        vec=posCA_sp(dd(find(dd)));
                        nc=find(vec==0 | vec>num_sp,1,'first');
                        if vec(nc)==0
                            Macc(dir(idir))=nc;
                            v=genvarname(strcat('vpos',num2str(dir(idir))));
                            eval([v '= dd(nc:-1:1);'])
                        end
                    end
                end

                [~,is]=sort(Macc); 
                if Macc(is(1))==10^30 %There is not space where to locate the daughter/displaced cell, therefore the cell cannot be divided. 
                    posCA_sp=spaceCA_orig;
                    IDpos=IDpos_orig;
                    Ncells_new=Ncells_neworig;
                else %Choose the direction with the nearest empty site. If the nearest empty site can found in more than 2 directions, choose one of them at random 
                    ri=randi(length(find(Macc==Macc(is(1)))),1); 
                    wd=is(ri);
                    vv=eval(strcat('vpos',num2str(wd)));

                    target=vv(1);
                    mediumLBM(IDca(target,Q+1))=2;
                    for p=2:length(vv)-1
                        %Update the position of the cell displaced.
                        sp_new=posCA_sp(vv(p));  
                        bb=find(IDpos(:,1)==vv(p)+Ny_ca*Nx_ca*Nz_ca*(sp_new-1));
                        IDpos(bb,1)=target+Ny_ca*Nx_ca*Nz_ca*(sp_new-1);
                        IDpos(bb,4)=IDca(target,Q+1);
                        posCA_sp(target)=sp_new;
                        target=vv(p);
                    end
                    
                    %Place the daugther cell (new cell)
                    posCA_sp(target)=sp_tc; 
                    Ncells_new=Ncells_new+1;
                    IDpos(Ncells_new,:)=data_tc;
                    IDpos(Ncells_new,1)=target+Ny_ca*Nx_ca*Nz_ca*(sp_tc-1);
                    IDpos(Ncells_new,4)=IDca(target,Q+1);
                    
                end
            end
        end
        
    elseif IDpos(rm,3)<=Mminsp(sp) %Cell will die 
        IDpos(rm,:)=0;
        posCA_sp(ijk)=0;
        
    else %Cell can move
        %If the cell is not motile, therefore it shows a Brownian motion
        direction=unidrnd(Q-1,1,1)+1;   %Choose the direction where the molecule will move  
        ijkt=IDca(ijk,direction);
        ijkt_lbm=IDca(ijk,Q+direction);

         if posCA_sp(ijkt)==0 && ijkt~=ijk  %If the target voxel is empty, then the cell will move, otherwise the cell will remain in its current position 
            if mediumLBM(ijkt_lbm)==2 || isempty(find(mediumLBM(IDca(ijkt,Q+2:2*Q))==3))==0 %Motion is allowed to areas with polysaccharide matrix or in surface areas with agar 
                mediumLBM(ijkt_lbm)=2;
                posCA_sp(ijk)=0;
                posCA_sp(ijkt)=sp;

                r2t(sp)=r2t(sp)+(Ax_ca^2); %Update the mean squared displacement of the species [mm^2].
                IDpos(rm,1)=ijkt+Ny_ca*Nx_ca*Nz_ca*(sp-1);
                IDpos(rm,4)=IDca(ijkt,(Q)+(direction));
                
            end
        end
    end
end

%Sort the IDpos matrix by species 
[~,hh]=sort(IDpos(:,2),'descend');
IDpos=IDpos(hh,:);
Ncells_new=length(find(posCA_sp(:)>0 & posCA_sp(:)<=num_sp));


% Compute the mean square displacement of the microbial species
for sp=1:num_sp
    cells1=length(find(posCA_sp==sp));
    r2t_sp(sp)=r2t_sp(sp)+(r2t(sp)/cells1);
end
