function [Q,IDlbm,IDca,FD,FBD,vol_ca,vol_lbm,posCA_sp,frac_medium,frac_cells,frac_intcells,Macc]=HexLatticeLBM(Nx_lbm,Ny_lbm,Nz_lbm,num_met,num_sp,Ax_ca,n,Nx_ca,Ny_ca,Nz_ca,dim,rhoCA_sp,posCA_sp,boundaryLBM,boundaryCA,CmetDb,v_sp,v_met,medium)
%Hexagonal grid

%Compute the volume of a CA-voxel and LBM-voxel
Area_lbm=2*sqrt(3)*((n*0.5*Ax_ca)^2); 
Area_ca=2*sqrt(3)*((0.5*Ax_ca)^2);
vol_lbm=Area_lbm*Ax_ca; %Volume of a LBM-voxel [mm^3].
vol_ca=Area_ca*Ax_ca;   %Volume of a CA-voxel [mm^3].

Q=(dim*3)+1;            %Posible directions of a particle.


%Compute index of neighbouring voxels in CA-lattice (stored in IDca) and
%LBM-lattice (stored in IDlbm)
IDlbm=zeros(Nx_lbm*Ny_lbm*Nz_lbm*num_met,Q);    %Matrix with the neighbouring positions of each voxel [dimensionless].
IDca=zeros(Nx_ca*Ny_ca*Nz_ca,Q);
rhoBD=zeros(Nx_lbm*Ny_lbm*Nz_lbm*num_met,1);    %Matrix with the amount of metabolites that enter to the system through the Dirichlet boundary [mmol].


ncol=max([Nx_ca,Ny_ca]);
Macc=zeros(Nx_ca*Ny_ca*Nz_ca,Q-1,ncol);
vM=repmat(-1:-1:-ceil((ncol-1)/2),2);

for ij=1:Nx_ca*Ny_ca
    row=ij-(ceil(ij/Ny_ca)-1)*Ny_ca;
    col=ceil(ij/Ny_ca);

    %Compute matrix Macc that contains the indexes of voxels located in one direction "d" in the hexagonal lattice. d=1...nvox
    if rem(col,2)==1 %odd column
        %Macc in direction 2
        col_final=(col+ncol>Nx_ca)*Nx_ca+(col+ncol<=Nx_ca)*(col+ncol);
        cf=col_final-col;
        vec=ij+(1:cf)*Ny_ca+[0,vM(1:cf-1)];
        rowvec=vec-(ceil(vec/Ny_ca)-1)*Ny_ca;
        con=find(rowvec==1);
        csal=find(vec<=Nx_ca*Ny_ca);
        if isempty(csal)==0 && isempty(con)==0 
            if length(con)>1
                cc=con(2);
            else
                cc=con(1);
            end
            cc1=min([cc csal(end)]);
            Macc(ij,1,1:cc1)=vec(1:cc1); 
        elseif isempty(con)==0 
            if length(con)>1
                cc=con(2);
            else
                cc=con(1);
            end
            Macc(ij,1,1:cc)=vec(1:cc); 
        elseif isempty(csal)==0 
            Macc(ij,1,1:csal(end))=vec(1:csal(end));
        end

        %Macc in direction 4
        col_final=(col-ncol>=1)*(col-ncol)+(col-ncol<1)*1;
        cf=col-col_final;
        vec=ij-(1:cf)*Ny_ca+[0,vM(1:cf-1)];
        rowvec=vec-(ceil(vec/Ny_ca)-1)*Ny_ca;
        con=find(rowvec==1);
        csal=find(vec>=1);
        if isempty(csal)==0 && isempty(con)==0 
            if length(con)>1
                cc=con(2);
            else
                cc=con(1);
            end
            cc1=min([cc csal(end)]);
            Macc(ij,3,1:cc1)=vec(1:cc1); 
        elseif isempty(con)==0 
            if length(con)>1
                cc=con(2);
            else
                cc=con(1);
            end
            Macc(ij,3,1:cc)=vec(1:cc); 
        elseif isempty(csal)==0 
            Macc(ij,3,1:csal(end))=vec(1:csal(end));
        end

        %Macc in direction 5
        col_final=(col-ncol>=1)*(col-ncol)+(col-ncol<1)*1;
        cf=col-col_final;
        vec=ij-(1:cf)*Ny_ca-vM(1:cf);
        rowvec=vec-(ceil(vec/Ny_ca)-1)*Ny_ca;
        con=find(rowvec==1);
        csal=find(vec>=1);
        if isempty(csal)==0 && isempty(con)==0 
            if length(con)>1
                cc=con(2);
            else
                cc=con(1);
            end
            cc1=min([cc csal(end)]);
            Macc(ij,4,1:cc1)=vec(1:cc1); 
        elseif isempty(con)==0 
            if length(con)>1
                cc=con(2);
            else
                cc=con(1);
            end
            Macc(ij,4,1:cc)=vec(1:cc); 
        elseif isempty(csal)==0 
            Macc(ij,4,1:csal(end))=vec(1:csal(end));
        end

        %Macc in direction 7
        col_final=(col+ncol>Nx_ca)*Nx_ca+(col+ncol<=Nx_ca)*(col+ncol);
        cf=col_final-col;
        vec=ij+(1:cf)*Ny_ca-vM(1:cf);
        rowvec=vec-(ceil(vec/Ny_ca)-1)*Ny_ca;
        con=find(rowvec==1);
        csal=find(vec<=Nx_ca*Ny_ca);
        if isempty(csal)==0 && isempty(con)==0 
            if length(con)>1
                cc=con(2);
            else
                cc=con(1);
            end
            cc1=min([cc csal(end)]);
            Macc(ij,6,1:cc1)=vec(1:cc1); 
        elseif isempty(con)==0 
            if length(con)>1
                cc=con(2);
            else
                cc=con(1);
            end
            Macc(ij,6,1:cc)=vec(1:cc); 
        elseif isempty(csal)==0 
            Macc(ij,6,1:csal(end))=vec(1:csal(end));
        end

    else %even column

        %Macc in direction 2
        col_final=(col+ncol>Nx_ca)*Nx_ca+(col+ncol<=Nx_ca)*(col+ncol);
        cf=col_final-col;
        vec=ij+(1:cf)*Ny_ca+vM(1:cf);
        rowvec=vec-(ceil(vec/Ny_ca)-1)*Ny_ca;
        con=find(rowvec==1);
        csal=find(vec<=Nx_ca*Ny_ca);
        if isempty(csal)==0 && isempty(con)==0 
            if length(con)>1
                cc=con(2);
            else
                cc=con(1);
            end
            cc1=min([cc csal(end)]);
            Macc(ij,1,1:cc1)=vec(1:cc1); 
        elseif isempty(con)==0 
            if length(con)>1
                cc=con(2);
            else
                cc=con(1);
            end
            Macc(ij,1,1:cc)=vec(1:cc); 
        elseif isempty(csal)==0 
            Macc(ij,1,1:csal(end))=vec(1:csal(end));
        end

        %Macc in direction 4
        col_final=(col-ncol>=1)*(col-ncol)+(col-ncol<1)*1;
        cf=col-col_final;
        vec=ij-(1:cf)*Ny_ca+vM(1:cf);
        rowvec=vec-(ceil(vec/Ny_ca)-1)*Ny_ca;
        con=find(rowvec==1);
        csal=find(vec>=1);
        if isempty(csal)==0 && isempty(con)==0 
            if length(con)>1
                cc=con(2);
            else
                cc=con(1);
            end
            cc1=min([cc csal(end)]);
            Macc(ij,3,1:cc1)=vec(1:cc1); 
        elseif isempty(con)==0 
            if length(con)>1
                cc=con(2);
            else
                cc=con(1);
            end
            Macc(ij,3,1:cc)=vec(1:cc); 
        elseif isempty(csal)==0 
            Macc(ij,3,1:csal(end))=vec(1:csal(end));
        end

        %Macc in direction 5
        col_final=(col-ncol>=1)*(col-ncol)+(col-ncol<1)*1;
        cf=col-col_final;
        vec=ij-(1:cf)*Ny_ca-[0,vM(1:cf-1)];
        rowvec=vec-(ceil(vec/Ny_ca)-1)*Ny_ca;
        con=find(rowvec==Ny_ca);
        csal=find(vec>=1);
        if isempty(csal)==0 && isempty(con)==0 
            if length(con)>1
                cc=con(2);
            else
                cc=con(1);
            end
            cc1=min([cc csal(end)]);
            Macc(ij,4,1:cc1)=vec(1:cc1); 
        elseif isempty(con)==0 
            if length(con)>1
                cc=con(2);
            else
                cc=con(1);
            end
            Macc(ij,4,1:cc)=vec(1:cc); 
        elseif isempty(csal)==0 
            Macc(ij,4,1:csal(end))=vec(1:csal(end));
        end

        %Macc in direction 7
        col_final=(col+ncol>Nx_ca)*Nx_ca+(col+ncol<=Nx_ca)*(col+ncol);
        cf=col_final-col;
        vec=ij+(1:cf)*Ny_ca-[0,vM(1:cf-1)];
        rowvec=vec-(ceil(vec/Ny_ca)-1)*Ny_ca;
        con=find(rowvec==1);
        csal=find(vec<=Nx_ca*Ny_ca);
        if isempty(csal)==0 && isempty(con)==0 
            if length(con)>1
                cc=con(2);
            else
                cc=con(1);
            end
            cc1=min([cc csal(end)]);
            Macc(ij,6,1:cc1)=vec(1:cc1); 
        elseif isempty(con)==0 
            if length(con)>1
                cc=con(2);
            else
                cc=con(1);
            end
            Macc(ij,6,1:cc)=vec(1:cc); 
        elseif isempty(csal)==0 
            Macc(ij,6,1:csal(end))=vec(1:csal(end));
        end
    end

    %Macc in direction 3
    row_final=(row-1<=ncol)*(row-1)+(row-1>ncol)*ncol;
    vec=ij-1:-1:ij-row_final;
    Macc(ij,2,1:length(vec))=vec;

    %Macc in direction 6
    row_final=(Ny_ca-row<=ncol)*(Ny_ca-row)+(Ny_ca-row>ncol)*ncol;
    vec=ij+1:ij+row_final;
    Macc(ij,5,1:length(vec))=vec;

    
   

    %Identify the neighboring CA-voxels for each voxel ijk:
    if rem(col,2)==0 
        IDca(ij,1)=ij;
        %Direction 2
        if ij+Ny_ca-1<=Ny_ca*Nx_ca && row>1 
            IDca(ij,2)=ij+Ny_ca-1;
            if posCA_sp(ij+Ny_ca-1)==num_sp+1 %Solid
                IDca(ij,2)=ij;
            end
        elseif row==1 && col==Nx_ca && boundaryCA.Y0=='P' || row==1 && col==Nx_ca && boundaryCA.Yn=='P'     
            IDca(ij,2)=Ny_ca;
        elseif col==Nx_ca && boundaryCA.X0=='P' || col==Nx_ca && boundaryCA.Xn=='P'     
            IDca(ij,2)=ij-(Nx_ca-1)*Ny_ca-1;
        elseif row==1 && boundaryCA.X0=='P' || row==1 && boundaryCA.Xn=='P'
            IDca(ij,2)=ij+(Ny_ca*2)-1;
        else %Bounceback boundary
            IDca(ij,2)=ij;
        end

        %Direction 3
        if row>1 
            IDca(ij,3)=ij-1;
            if posCA_sp(ij-1)==num_sp+1 %Solid
                IDca(ij,3)=ij;
            end
        elseif row==1 && boundaryCA.Y0=='P' || row==1 && boundaryCA.Yn=='P'
            IDca(ij,3)=ij+Ny_ca-1;
        else %Bounceback boundary
            IDca(ij,3)=ij;
        end

        %Direction 4
        if ij-Ny_ca-1>0 && row>1 
            IDca(ij,4)=ij-Ny_ca-1;
            if posCA_sp(ij-Ny_ca-1)==num_sp+1 %Solid
                IDca(ij,4)=ij;
            end
        elseif row==1 && boundaryCA.Y0=='P' || row==1 && boundaryCA.Yn=='P'
            IDca(ij,4)=ij-1;
        else %Bounceback boundary
            IDca(ij,4)=ij;
        end

        %Direction 5
        IDca(ij,5)=ij-Ny_ca;
        if posCA_sp(ij-Ny_ca)==num_sp+1 %Solid
            IDca(ij,5)=ij;
        end

        %Direction 6
        if row<Ny_ca 
            IDca(ij,6)=ij+1;
            if posCA_sp(ij+1)==num_sp+1 %Solid
                IDca(ij,6)=ij;
            end
        elseif row==Ny_ca && boundaryCA.Y0=='P' || row==Ny_ca && boundaryCA.Yn=='P'
            IDca(ij,6)=ij-Ny_ca+1;
        else %Bounceback boundary
            IDca(ij,6)=ij;
        end

        %Direction 7
        if col<Nx_ca
            IDca(ij,7)=ij+Ny_ca;
            if posCA_sp(ij+Ny_ca)==num_sp+1 %Solid
                IDca(ij,7)=ij;
            end
        elseif col==Nx_ca && boundaryCA.X0=='P' || row==Ny_ca && boundaryCA.Xn=='P'
            IDca(ij,7)=ij-Ny_ca*(col-1);
        else
            IDca(ij,7)=ij;
        end

    else
        IDca(ij,1)=ij;

        %Direction 2
        IDca(ij,2)=ij+Ny_ca;
        if posCA_sp(ij+Ny_ca)==num_sp+1 %Solid
            IDca(ij,2)=ij;
        end

        %Direction 3
        if row>1 
            IDca(ij,3)=ij-1;
            if posCA_sp(ij-1)==num_sp+1 %Solid
                IDca(ij,3)=ij;
            end
        elseif row==1 && boundaryCA.Y0=='P' || row==1 && boundaryCA.Yn=='P'
            IDca(ij,3)=ij+Ny_ca-1;
        else %Bounceback boundary
            IDca(ij,3)=ij;
        end

        %Direction 4
        if col>1 
            IDca(ij,4)=ij-Ny_ca;
            if posCA_sp(ij-Ny_ca)==num_sp+1 %Solid
                IDca(ij,4)=ij;
            end
        elseif col==1 && boundaryCA.X0=='P' || col==1 && boundaryCA.Xn=='P'
            IDca(ij,4)=ij+(Nx_ca-1)*Ny_ca;
        else %Bounceback boundary
            IDca(ij,4)=ij;
        end

        %Direction 5
        if col>1 && row<Ny_ca
            IDca(ij,5)=ij-Ny_ca+1;
            if posCA_sp(ij-Ny_ca+1)==num_sp+1 %Solid
                IDca(ij,5)=ij;
            end
        elseif ij==Ny_ca && boundaryCA.Y0=='P' || ij==Ny_ca && boundaryCA.Yn=='P'
            IDca(ij,5)=(Nx_ca-1)*Ny_ca+1;    
        elseif col==1 && boundaryCA.X0=='P' || col==1 && boundaryCA.Xn=='P'
            IDca(ij,5)=(Nx_ca-1)*Ny_ca+ij+1;
        elseif row==Ny_ca && boundaryCA.Y0=='P' || row==Ny_ca && boundaryCA.Yn=='P'
            IDca(ij,5)=(col-2)*Ny_ca+1;
        else
            IDca(ij,5)=ij;
        end

        %Direction 6
        if row<Ny_ca 
            IDca(ij,6)=ij+1;
            if posCA_sp(ij+1)==num_sp+1 %Solid
                IDca(ij,6)=ij;
            end
        elseif row==Ny_ca && boundaryCA.Y0=='P' || row==Ny_ca && boundaryCA.Yn=='P'
            IDca(ij,6)=ij-Ny_ca+1;
        else %Bounceback boundary
            IDca(ij,6)=ij;
        end

        %Direction 7
        if row<Ny_ca
            IDca(ij,7)=ij+Ny_ca+1;
            if posCA_sp(ij+Ny_ca+1)==num_sp+1 %Solid
                IDca(ij,7)=ij;
            end
        elseif row==Ny_ca && boundaryCA.Y0=='P' || row==Ny_ca && boundaryCA.Yn=='P'
            IDca(ij,7)=ij+1;
        else
            IDca(ij,7)=ij;
        end
    end


    
    %Identify the neighboring LBM-voxels for each voxel ijk:
    row=ij-(ceil(ij/Ny_lbm)-1)*Ny_lbm;
    col=ceil(ij/Ny_lbm);

    for met=1:num_met
        idvec=Nx_lbm*Ny_lbm*(met-1);

        if rem(col,2)==0 

            IDlbm(ij+idvec,1)=ij+idvec;

            %Direction 2
            if ij+Ny_lbm-1<=Ny_lbm*Nx_lbm && row>1 
                IDlbm(ij+idvec,2)=ij+Ny_lbm-1+idvec;
            elseif row==1 && col==Nx_lbm && boundaryLBM.Y0{met}=='P' || row==1 && col==Nx_lbm && boundaryLBM.Yn{met}=='P'     
                IDlbm(ij+idvec,2)=Ny_lbm+idvec;
            elseif row==1 && boundaryLBM.Y0{met}=='D'   
                IDlbm(ij+idvec,2)=Nx_lbm*Ny_lbm*num_met+2;
            elseif col==Nx_lbm && boundaryLBM.X0{met}=='P' || col==Nx_lbm && boundaryLBM.Xn{met}=='P'     
                IDlbm(ij+idvec,2)=ij-(Nx_lbm-1)*Ny_lbm-1+idvec;
            elseif col==Nx_lbm && boundaryLBM.Xn{met}=='D'  
                IDlbm(ij+idvec,2)=Nx_lbm*Ny_lbm*num_met+2;
            elseif row==1 && boundaryLBM.X0{met}=='P' || row==1 && boundaryLBM.Xn{met}=='P'
                IDlbm(ij+idvec,2)=ij+(Ny_lbm*2)-1+idvec;
            else %Bounceback boundary
                IDlbm(ij+idvec,2)=Nx_lbm*Ny_lbm*num_met+1;
            end

            %Direction 3
            if row>1 
                IDlbm(ij+idvec,3)=ij-1+idvec;
            elseif row==1 && boundaryLBM.Y0{met}=='P' || row==1 && boundaryLBM.Yn{met}=='P'
                IDlbm(ij+idvec,3)=ij+Ny_lbm-1+idvec;
            elseif row==1 && boundaryLBM.Y0{met}=='D' 
                IDlbm(ij+idvec,3)=Nx_lbm*Ny_lbm*num_met+2;
            else %Bounceback boundary
                IDlbm(ij+idvec,3)=Nx_lbm*Ny_lbm*num_met+1;
            end

            %Direction 4
            if ij-Ny_lbm-1>0 && row>1 
                IDlbm(ij+idvec,4)=ij-Ny_lbm-1+idvec;
            elseif row==1 && boundaryLBM.Y0{met}=='P' || row==1 && boundaryLBM.Yn{met}=='P'
                IDlbm(ij+idvec,4)=ij-1+idvec;
            elseif row==1 && boundaryLBM.Y0{met}=='D' 
                IDlbm(ij+idvec,4)=Nx_lbm*Ny_lbm*num_met+2;
            else %Bounceback boundary
                IDlbm(ij+idvec,4)=Nx_lbm*Ny_lbm*num_met+1;
            end

            %Direction 5
            IDlbm(ij+idvec,5)=ij-Ny_lbm+idvec;

            %Direction 6
            if row<Ny_lbm 
                IDlbm(ij+idvec,6)=ij+1+idvec;
            elseif row==Ny_lbm && boundaryLBM.Y0{met}=='P' || row==Ny_lbm && boundaryLBM.Yn{met}=='P'
                IDlbm(ij+idvec,6)=ij-Ny_lbm+1+idvec;
            elseif row==Ny_lbm && boundaryLBM.Yn{met}=='D' 
                IDlbm(ij+idvec,6)=Nx_lbm*Ny_lbm*num_met+2;
            else %Bounceback boundary
                IDlbm(ij+idvec,6)=Nx_lbm*Ny_lbm*num_met+1;
            end

            %Direction 7
            if col<Nx_lbm
                IDlbm(ij+idvec,7)=ij+Ny_lbm+idvec;
            elseif col==Nx_lbm && boundaryLBM.X0{met}=='P' || row==Ny_lbm && boundaryLBM.Xn{met}=='P'
                IDlbm(ij+idvec,7)=ij-Ny_lbm*(col-1)+idvec;
            elseif col==Nx_lbm && boundaryLBM.Xn{met}=='D' 
                IDlbm(ij+idvec,7)=Nx_lbm*Ny_lbm*num_met+2;
            else
                IDlbm(ij+idvec,7)=Nx_lbm*Ny_lbm*num_met+1;
            end

        else

            IDlbm(ij+idvec,1)=ij+idvec;

            %Direction 2
            IDlbm(ij+idvec,2)=ij+Ny_lbm+idvec;

            %Direction 3
            if row>1 
                IDlbm(ij+idvec,3)=ij-1+idvec;
            elseif row==1 && boundaryLBM.Y0{met}=='P' || row==1 && boundaryLBM.Yn{met}=='P'
                IDlbm(ij+idvec,3)=ij+Ny_lbm-1+idvec;
            elseif row==1 && boundaryLBM.Y0{met}=='D' 
                IDlbm(ij+idvec,3)=Nx_lbm*Ny_lbm*num_met+2;
            else %Bounceback boundary
                IDlbm(ij+idvec,3)=Nx_lbm*Ny_lbm*num_met+1;
            end

            %Direction 4
            if col>1 
                IDlbm(ij+idvec,4)=ij-Ny_lbm+idvec;
            elseif col==1 && boundaryLBM.X0{met}=='P' || col==1 && boundaryLBM.Xn{met}=='P'
                IDlbm(ij+idvec,4)=ij+(Nx_lbm-1)*Ny_lbm+idvec;
            elseif col==1 && boundaryLBM.X0{met}=='D' 
                IDlbm(ij+idvec,4)=Nx_lbm*Ny_lbm*num_met+2;
            else %Bounceback boundary
                IDlbm(ij+idvec,4)=Nx_lbm*Ny_lbm*num_met+1;
            end

            %Direction 5
            if col>1 && row<Ny_lbm
                IDlbm(ij+idvec,5)=ij-Ny_lbm+1+idvec;
            elseif ij==Ny_lbm && boundaryLBM.Y0{met}=='P' || ij==Ny_lbm && boundaryLBM.Yn{met}=='P'
                IDlbm(ij+idvec,5)=(Nx_lbm-1)*Ny_lbm+1+idvec; 
            elseif ij==Ny_lbm && boundaryLBM.Yn{met}=='D' 
                IDlbm(ij+idvec,5)=Nx_lbm*Ny_lbm*num_met+2; 
            elseif col==1 && boundaryLBM.X0{met}=='P' || col==1 && boundaryLBM.Xn{met}=='P'
                IDlbm(ij+idvec,5)=(Nx_lbm-1)*Ny_lbm+ij+1+idvec;
            elseif col==1 && boundaryLBM.X0{met}=='D' 
                IDlbm(ij+idvec,5)=Nx_lbm*Ny_lbm*num_met+2;
            elseif row==Ny_lbm && boundaryLBM.Y0{met}=='P' || row==Ny_lbm && boundaryLBM.Yn{met}=='P'
                IDlbm(ij+idvec,5)=(col-2)*Ny_lbm+1+idvec;
            elseif row==Ny_lbm && boundaryLBM.Yn{met}=='D' 
                IDlbm(ij+idvec,5)=Nx_lbm*Ny_lbm*num_met+2;
            else
                IDlbm(ij+idvec,5)=Nx_lbm*Ny_lbm*num_met+1;
            end

            %Direction 6
            if row<Ny_lbm 
                IDlbm(ij+idvec,6)=ij+1+idvec;
            elseif row==Ny_lbm && boundaryLBM.Y0{met}=='P' || row==Ny_lbm && boundaryLBM.Yn{met}=='P'
                IDlbm(ij+idvec,6)=ij-Ny_lbm+1+idvec;
            elseif row==Ny_lbm && boundaryLBM.Yn{met}=='D' 
                IDlbm(ij+idvec,6)=Nx_lbm*Ny_lbm*num_met+2;
            else %Bounceback boundary
                IDlbm(ij+idvec,6)=Nx_lbm*Ny_lbm*num_met+1;
            end

            %Direction 7
            if row<Ny_lbm
                IDlbm(ij+idvec,7)=ij+Ny_lbm+1+idvec;
            elseif row==Ny_lbm && boundaryLBM.Y0{met}=='P' || row==Ny_lbm && boundaryLBM.Yn{met}=='P'
                IDlbm(ij+idvec,7)=ij+1+idvec;
            elseif row==Ny_lbm && boundaryLBM.Yn{met}=='D' 
                IDlbm(ij+idvec,7)=Nx_lbm*Ny_lbm*num_met+2;
            else
                IDlbm(ij+idvec,7)=Nx_lbm*Ny_lbm*num_met+1;
            end
        end

    end
end






FDx=zeros(Nx_lbm*Ny_lbm*Nz_lbm*num_met,1);      %Matrix with the amount of metabolites that enter to the system through the Dirichlet boundary [mmol].
FDy=zeros(Nx_lbm*Ny_lbm*Nz_lbm*num_met,1);      %Matrix with the amount of metabolites that enter to the system through the Dirichlet boundary [mmol].
FDz=zeros(Nx_lbm*Ny_lbm*Nz_lbm*num_met,1);      %Matrix with the amount of metabolites that enter to the system through the Dirichlet boundary [mmol].
rhoBDx=zeros(Nx_lbm*Ny_lbm*Nz_lbm*num_met,1);   %Matrix with the amount of metabolites that enter to the system through the Dirichlet boundary [mmol].
rhoBDy=zeros(Nx_lbm*Ny_lbm*Nz_lbm*num_met,1);   %Matrix with the amount of metabolites that enter to the system through the Dirichlet boundary [mmol].
rhoBDz=zeros(Nx_lbm*Ny_lbm*Nz_lbm*num_met,1);   %Matrix with the amount of metabolites that enter to the system through the Dirichlet boundary [mmol].



[i,j,k]=ind2sub([Ny_lbm,Nx_lbm,Nz_lbm],(1:Ny_lbm*Nx_lbm*Nz_lbm)');
for met=1:num_met
    
    %Vector of metabolites supply through Dirichlet boundaries 
    idvec=Nx_lbm*Ny_lbm*Nz_lbm*(met-1)+(1:Nx_lbm*Ny_lbm*Nz_lbm)';
    %Identifying the Dirichlet boundaries 
    vecDx=1-((i-1)/(Ny_lbm-1)); %gradient in Y direction
    rhoBDx(idvec)=(CmetDb.X0(met)*vol_lbm).*(j==1 & rhoBDx(idvec)<CmetDb.X0(met)*vol_lbm & strcmp(boundaryLBM.X0{met},'D'))+(CmetDb.Xn(met)*vol_lbm).*(j==Nx_lbm & rhoBDx(idvec)<CmetDb.Xn(met)*vol_lbm & strcmp(boundaryLBM.Xn{met},'D'))+(CmetDb.X0(met)*vol_lbm).*vecDx.*(j==1 & rhoBDx(idvec)<CmetDb.X0(met)*vol_lbm & strcmp(boundaryLBM.X0{met},'Dg'))+(CmetDb.Xn(met)*vol_lbm).*vecDx.*(j==Nx_lbm & rhoBDx(idvec)<CmetDb.Xn(met)*vol_lbm & strcmp(boundaryLBM.Xn{met},'Dg'));
    FDx(idvec)=(j==1 & strcmp(boundaryLBM.X0{met},'D'))+(j==Nx_lbm & strcmp(boundaryLBM.Xn{met},'D'))+(j==1 & strcmp(boundaryLBM.X0{met},'Dg'))+(j==Nx_lbm & strcmp(boundaryLBM.Xn{met},'Dg'));
    
    vecDy=1-((j-1)/(Nx_lbm-1)); %gradient in X direction
    rhoBDy(idvec)=(CmetDb.Y0(met)*vol_lbm).*(i==1 & rhoBDy(idvec)<CmetDb.Y0(met)*vol_lbm & strcmp(boundaryLBM.Y0{met},'D'))+(CmetDb.Yn(met)*vol_lbm).*(i==Ny_lbm & rhoBDy(idvec)<CmetDb.Yn(met)*vol_lbm & strcmp(boundaryLBM.Yn{met},'D'))+(CmetDb.Y0(met)*vol_lbm).*vecDy.*(i==1 & rhoBDy(idvec)<CmetDb.Y0(met)*vol_lbm & strcmp(boundaryLBM.Y0{met},'Dg'))+(CmetDb.Yn(met)*vol_lbm).*vecDy.*(i==Ny_lbm & rhoBDy(idvec)<CmetDb.Yn(met)*vol_lbm & strcmp(boundaryLBM.Yn{met},'Dg'));
    FDy(idvec)=(i==1 & strcmp(boundaryLBM.Y0{met},'D'))+(i==Ny_lbm & strcmp(boundaryLBM.Yn{met},'D'))+(i==1 & strcmp(boundaryLBM.Y0{met},'Dg'))+(i==Ny_lbm & strcmp(boundaryLBM.Yn{met},'Dg'));
  
    
end

FD=FDx+FDy+FDz;
index=find(FD);
FD(index)=1; %Index vector that ponits (when FD=1) voxels located in the Dirichlet boundary.

FBD=rhoBDx+rhoBDy+rhoBDz; %Vector with the amount molecules that enter in the system through the Dirichlet boundary [mmol].


%IDca contains the index of the CA-voxel and the corresponding LBM-voxel.
%For hexagonal grid, both CA- and LBM- lattices are the same.
IDca=[IDca,IDca];


%Compute posCA_sp, a position matrix-CA containing the distribution of species and solid in the CA-lattice 
frac_cells=zeros(Ny_lbm,Nx_lbm,Nz_lbm);     %Volume fraction of the LBM-voxel that is occupied by cells [dimensionless]
frac_intcells=zeros(Ny_lbm,Nx_lbm,Nz_lbm);  %Volume fraction of the LBM-voxel that is occupied by intracellular proteins [dimensionless]
for sp=1:num_sp
    xyz=find(rhoCA_sp(Nx_ca*Ny_ca*Nz_ca*(sp-1)+1:Nx_ca*Ny_ca*Nz_ca*sp));
    posCA_sp(xyz)=sp;
    ijk=IDca(xyz,Q+1);
    X=rhoCA_sp(xyz)'*v_sp(sp)/vol_lbm;
    frac_cells(1:max(ijk))=frac_cells(1:max(ijk))+accumarray(ijk,X(:))';
    Xintcel=rhoCA_sp(xyz)'*v_met/vol_lbm;
    frac_intcells(1:max(ijk))=frac_intcells(1:max(ijk))+accumarray(ijk,Xintcel(:))';
end


%Compute frac_medium
frac_medium=zeros(Ny_lbm,Nx_lbm,Nz_lbm,num_met);    %Matrix-LBM with the volume fraction NOT ocuppied by the porous medium [dimensionless].
xyz=find(medium==1); %solids
ijk=IDca(xyz,Q+1);
frac_solids=zeros(Ny_lbm,Nx_lbm,Nz_lbm);
frac_solids(1:max(ijk))=accumarray(ijk,(xyz>0).*vol_ca/vol_lbm);
frac_medium(:)=repmat(1-frac_solids(:),num_met,1);

