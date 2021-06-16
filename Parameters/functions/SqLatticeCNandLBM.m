function [Q,IDlbm,IDca,FD,FBD,vol_ca,vol_lbm,posCA_sp,frac_medium,frac_cells,frac_intcells,a1y,c1y,a1x,c1x,a1z,c1z,b1y,b2y,b1x,b2x,b1z,b2z,ap1x,cp1x,ap1y,cp1y,ap1z,cp1z,CNx,CNy,CNz]=SqLatticeCNandLBM(Nx_lbm,Ny_lbm,Nz_lbm,num_met,num_sp,Ax_ca,Ay_ca,Az_ca,Ax_lbm,Ay_lbm,Az_lbm,nx,ny,nz,Nx_ca,Ny_ca,Nz_ca,dim,rhoCA_sp,posCA_sp,boundaryLBM,boundaryCA,CmetDb,v_sp,v_met,medium)
%Square or cubic grid

%Compute the volume of a CA-voxel and LBM-voxel
if dim==2
    Area_lbm=Ax_lbm*Ay_lbm;
    vol_lbm=Area_lbm*Az_lbm; %Volume of a LBM-voxel [mm^3].
else
    Area_lbm=Ax_lbm*Az_lbm;
    vol_lbm=Area_lbm*Ay_lbm;    
end
vol_ca=Ax_ca*Ay_ca*Az_ca;    %Volume of a CA-voxel [mm^3].

Q=(dim*2)+1;   %Posible directions of a particle.


%Compute index of neighbouring voxels in CA-lattice and store it in IDca
IDca=zeros(Nx_ca*Ny_ca*Nz_ca,Q*2); %Directions are: [(i,j,k),(i-1,j,k);(i+1,j,k);(i,j-1,k);(i,j+1,k);(i,j,k-1);(i,j,k+1)]

IDca(:,1)=(1:Nx_ca*Ny_ca*Nz_ca)';   
[i,j,k]=ind2sub([Ny_ca,Nx_ca,Nz_ca],IDca(:,1));
IDca(:,Q+1)=sub2ind([Ny_lbm,Nx_lbm,Nz_lbm],ceil(i/ny),ceil(j/nx),ceil(k/nz));        

jd=j; kd=k;
ijk_prov=sub2ind([Ny_ca,Nx_ca,Nz_ca],(i-1>0).*(i-1)+(i-1<=0).*i,j,k);
id=(i-1).*(i-1>0 & i-1<=Ny_ca & posCA_sp(ijk_prov)~=num_sp+1)+(i).*(i-1>0 & i-1<=Ny_ca & posCA_sp(ijk_prov)==num_sp+1)+(Ny_ca).*((i-1)==0 & boundaryCA.Y0=='P' | (i-1)==0 & boundaryCA.Yn=='P')+(i).*((i-1)==0 & boundaryCA.Y0=='B' | (i-1)==0 & boundaryCA.Yn=='B');
IDca(:,2)=sub2ind([Ny_ca,Nx_ca,Nz_ca],id,jd,kd);      
IDca(:,Q+2)=sub2ind([Ny_lbm,Nx_lbm,Nz_lbm],ceil(id/ny),ceil(jd/nx),ceil(kd/nz));  
    
jd=j; kd=k;
ijk_prov=sub2ind([Ny_ca,Nx_ca,Nz_ca],(i+1<=Ny_ca).*(i+1)+(i+1>Ny_ca).*i,j,k);
id=(i+1).*(i+1>0 & i+1<=Ny_ca & posCA_sp(ijk_prov)~=num_sp+1)+(i).*(i+1>0 & i+1<=Ny_ca & posCA_sp(ijk_prov)==num_sp+1)+(1).*((i+1)>Ny_ca & boundaryCA.Y0=='P' | (i+1)>Ny_ca & boundaryCA.Yn=='P')+(i).*((i+1)>Ny_ca & boundaryCA.Y0=='B' | (i+1)>Ny_ca & boundaryCA.Yn=='B');
IDca(:,3)=sub2ind([Ny_ca,Nx_ca,Nz_ca],id,jd,kd);      
IDca(:,Q+3)=sub2ind([Ny_lbm,Nx_lbm,Nz_lbm],ceil(id/ny),ceil(jd/nx),ceil(kd/nz));  
     
id=i; kd=k;
ijk_prov=sub2ind([Ny_ca,Nx_ca,Nz_ca],i,(j-1>0).*(j-1)+(j-1<=0).*j,k);
jd=(j-1).*(j-1>0 & j-1<=Nx_ca & posCA_sp(ijk_prov)~=num_sp+1)+(j).*(j-1>0 & j-1<=Nx_ca & posCA_sp(ijk_prov)==num_sp+1)+(Nx_ca).*((j-1)==0 & boundaryCA.X0=='P' | (j-1)==0 & boundaryCA.Xn=='P')+(j).*((j-1)==0 & boundaryCA.X0=='B' | (j-1)==0 & boundaryCA.Xn=='B');
IDca(:,4)=sub2ind([Ny_ca,Nx_ca,Nz_ca],id,jd,kd);      
IDca(:,Q+4)=sub2ind([Ny_lbm,Nx_lbm,Nz_lbm],ceil(id/ny),ceil(jd/nx),ceil(kd/nz)); 

id=i; kd=k;
ijk_prov=sub2ind([Ny_ca,Nx_ca,Nz_ca],i,(j+1<=Nx_ca).*(j+1)+(j+1>Nx_ca).*j,k);
jd=(j+1).*(j+1>0 & j+1<=Nx_ca & posCA_sp(ijk_prov)~=num_sp+1)+(j).*(j+1>0 & j+1<=Nx_ca & posCA_sp(ijk_prov)==num_sp+1)+(1).*((j+1)>Nx_ca & boundaryCA.X0=='P' | (j+1)>Nx_ca & boundaryCA.Xn=='P')+(j).*((j+1)>Nx_ca & boundaryCA.X0=='B' | (j+1)>Nx_ca & boundaryCA.Xn=='B');
IDca(:,5)=sub2ind([Ny_ca,Nx_ca,Nz_ca],id,jd,kd);      
IDca(:,Q+5)=sub2ind([Ny_lbm,Nx_lbm,Nz_lbm],ceil(id/ny),ceil(jd/nx),ceil(kd/nz));  

if dim==3
    id=i; jd=j; 
    ijk_prov=sub2ind([Ny_ca,Nx_ca,Nz_ca],i,j,(k-1>0).*(k-1)+(k-1<=0).*k);
    kd=(k-1).*(k-1>0 & k-1<=Nz_ca & posCA_sp(ijk_prov)~=num_sp+1)+(k).*(k-1>0 & k-1<=Nz_ca & posCA_sp(ijk_prov)==num_sp+1)+(Nz_ca).*((k-1)==0 & boundaryCA.Z0=='P' | (k-1)==0 & boundaryCA.Zn=='P')+(k).*((k-1)==0 & boundaryCA.Z0=='B' | (k-1)==0 & boundaryCA.Zn=='B');
    IDca(:,6)=sub2ind([Ny_ca,Nx_ca,Nz_ca],id,jd,kd);      
    IDca(:,Q+6)=sub2ind([Ny_lbm,Nx_lbm,Nz_lbm],ceil(id/ny),ceil(jd/nx),ceil(kd/nz));  

    id=i; jd=j; 
    ijk_prov=sub2ind([Ny_ca,Nx_ca,Nz_ca],i,j,(k+1<=Nz_ca).*(k+1)+(k+1>Nz_ca).*k);
    kd=(k+1).*(k+1>0 & k+1<=Nz_ca & posCA_sp(ijk_prov)~=num_sp+1)+(k).*(k+1>0 & k+1<=Nz_ca & posCA_sp(ijk_prov)==num_sp+1)+(1).*((k+1)>Nz_ca & boundaryCA.Z0=='P' | (k+1)>Nz_ca & boundaryCA.Zn=='P')+(k).*((k+1)>Nz_ca & boundaryCA.Z0=='B' | (k+1)>Nz_ca & boundaryCA.Zn=='B');
    IDca(:,7)=sub2ind([Ny_ca,Nx_ca,Nz_ca],id,jd,kd);      
    IDca(:,Q+7)=sub2ind([Ny_lbm,Nx_lbm,Nz_lbm],ceil(id/ny),ceil(jd/nx),ceil(kd/nz));  
end



%Compute index of neighbouring voxel in LBM-lattice and store it in IDlbm
IDlbm=zeros(Nx_lbm*Ny_lbm*Nz_lbm*num_met,Q);    %Matrix with the neighbouring positions of each voxel [dimensionless].
IDlbm(:,1)=(1:Nx_lbm*Ny_lbm*Nz_lbm*num_met)';
IDlbm(:,2)=(1:Nx_lbm*Ny_lbm*Nz_lbm*num_met)'+Ny_lbm;
IDlbm(:,3)=(1:Nx_lbm*Ny_lbm*Nz_lbm*num_met)'-1;
IDlbm(:,4)=(1:Nx_lbm*Ny_lbm*Nz_lbm*num_met)'-Ny_lbm;
IDlbm(:,5)=(1:Nx_lbm*Ny_lbm*Nz_lbm*num_met)'+1;
if dim==3
    IDlbm(:,6)=(1:Nx_lbm*Ny_lbm*Nz_lbm*num_met)'+Ny_lbm*Nx_lbm;
    IDlbm(:,7)=(1:Nx_lbm*Ny_lbm*Nz_lbm*num_met)'-Ny_lbm*Nx_lbm;
end



FDx=zeros(Nx_lbm*Ny_lbm*Nz_lbm*num_met,1);      %Matrix with the amount of metabolites that enter to the system through the Dirichlet boundary [mmol].
FDy=zeros(Nx_lbm*Ny_lbm*Nz_lbm*num_met,1);      %Matrix with the amount of metabolites that enter to the system through the Dirichlet boundary [mmol].
FDz=zeros(Nx_lbm*Ny_lbm*Nz_lbm*num_met,1);      %Matrix with the amount of metabolites that enter to the system through the Dirichlet boundary [mmol].
rhoBDx=zeros(Nx_lbm*Ny_lbm*Nz_lbm*num_met,1);   %Matrix with the amount of metabolites that enter to the system through the Dirichlet boundary [mmol].
rhoBDy=zeros(Nx_lbm*Ny_lbm*Nz_lbm*num_met,1);   %Matrix with the amount of metabolites that enter to the system through the Dirichlet boundary [mmol].
rhoBDz=zeros(Nx_lbm*Ny_lbm*Nz_lbm*num_met,1);   %Matrix with the amount of metabolites that enter to the system through the Dirichlet boundary [mmol].
periodic=zeros(Nx_lbm*Ny_lbm*Nz_lbm*num_met,Q);    %Matrix with the neighbouring positions of each voxel [dimensionless], used to compute the vectors required by Crank-Nicholson method 



[i,j,k]=ind2sub([Ny_lbm,Nx_lbm,Nz_lbm],(1:Ny_lbm*Nx_lbm*Nz_lbm)');
for met=1:num_met
    
    if strcmp(boundaryLBM.X0{met},'P')==1 || strcmp(boundaryLBM.Xn{met},'P')==1
        boundaryLBM.X0{met}='P';
        boundaryLBM.Xn{met}='P';
    end
    if strcmp(boundaryLBM.Y0{met},'P')==1 || strcmp(boundaryLBM.Yn{met},'P')==1
        boundaryLBM.Y0{met}='P';
        boundaryLBM.Yn{met}='P';
    end
    if strcmp(boundaryLBM.Z0{met},'P')==1 || strcmp(boundaryLBM.Zn{met},'P')==1
        boundaryLBM.Z0{met}='P';
        boundaryLBM.Zn{met}='P';
    end

    
    
    ycont1=(1:Nx_lbm*Ny_lbm*Nz_lbm)'+Nx_lbm*Ny_lbm*Nz_lbm*(met-1);
     
    j_1=Nx_lbm.*(j==1 & strcmp(boundaryLBM.X0{met},'P') | j==1 & strcmp(boundaryLBM.Xn{met},'P'))+0.*(j==1 & strcmp(boundaryLBM.X0{met},'B'))+(-1).*(j==1 & strcmp(boundaryLBM.X0{met},'D'))+(-1).*(j==1 & strcmp(boundaryLBM.X0{met},'Dg'))+(j-1).*(j~=1);
    j1=1.*(j==Nx_lbm & strcmp(boundaryLBM.X0{met},'P') | j==Nx_lbm & strcmp(boundaryLBM.Xn{met},'P'))+0.*(j==Nx_lbm & strcmp(boundaryLBM.Xn{met},'B'))+(-1).*(j==Nx_lbm & strcmp(boundaryLBM.Xn{met},'D'))+(-1).*(j==Nx_lbm & strcmp(boundaryLBM.Xn{met},'Dg'))+(j+1).*(j~=Nx_lbm);
    
    periodic(ycont1,4)=(ycont1+Ny_lbm*(Nx_lbm-1)).*(j==1 & strcmp(boundaryLBM.X0{met},'P') | j==1 & strcmp(boundaryLBM.Xn{met},'P'));
    periodic(ycont1,2)=(ycont1-(Ny_lbm*(Nx_lbm-1))).*(j==Nx_lbm & strcmp(boundaryLBM.X0{met},'P') | j==Nx_lbm & strcmp(boundaryLBM.Xn{met},'P'));   
    
    i_1=Ny_lbm.*(i==1 & strcmp(boundaryLBM.Y0{met},'P') | i==1 & strcmp(boundaryLBM.Yn{met},'P'))+0.*(i==1 & strcmp(boundaryLBM.Y0{met},'B'))+(-1).*(i==1 & strcmp(boundaryLBM.Y0{met},'D'))+(-1).*(i==1 & strcmp(boundaryLBM.Y0{met},'Dg'))+(i-1).*(i~=1);
    i1=1.*(i==Ny_lbm & strcmp(boundaryLBM.Y0{met},'P') | i==Ny_lbm & strcmp(boundaryLBM.Yn{met},'P'))+0.*(i==Ny_lbm & strcmp(boundaryLBM.Yn{met},'B'))+(-1).*(i==Ny_lbm & strcmp(boundaryLBM.Yn{met},'D'))+(-1).*(i==Ny_lbm & strcmp(boundaryLBM.Yn{met},'Dg'))+(i+1).*(i~=Ny_lbm);
    
    periodic(ycont1,3)=(ycont1+Ny_lbm-1).*(i==1 & strcmp(boundaryLBM.Y0{met},'P') | i==1 & strcmp(boundaryLBM.Yn{met},'P'));
    periodic(ycont1,5)=(ycont1-(Ny_lbm-1)).*(i==Ny_lbm & strcmp(boundaryLBM.Y0{met},'P') | i==Ny_lbm & strcmp(boundaryLBM.Yn{met},'P'));
    
    
    if dim==3
        k_1=Nz_lbm.*(k==1 & strcmp(boundaryLBM.Z0{met},'P') | k==1 & strcmp(boundaryLBM.Zn{met},'P'))+0.*(k==1 & strcmp(boundaryLBM.Z0{met},'B'))+(-1).*(k==1 & strcmp(boundaryLBM.Z0{met},'D'))+(-1).*(k==1 & strcmp(boundaryLBM.Z0{met},'Dg'))+(k-1).*(k~=1);
        k1=1.*(k==Nz_lbm & strcmp(boundaryLBM.Z0{met},'P') | k==Nz_lbm & strcmp(boundaryLBM.Zn{met},'P'))+0.*(k==Nz_lbm & strcmp(boundaryLBM.Zn{met},'B'))+(-1).*(k==Nz_lbm & strcmp(boundaryLBM.Zn{met},'D'))+(-1).*(k==Nz_lbm & strcmp(boundaryLBM.Zn{met},'Dg'))+(k+1).*(k~=Nz_lbm);
        
        periodic(ycont1,7)=(ycont1+Nx_lbm*Ny_lbm*(Nz_lbm-1)).*(k==1 & strcmp(boundaryLBM.Z0{met},'P') | k==1 & strcmp(boundaryLBM.Zn{met},'P'));
        periodic(ycont1,6)=(ycont1-(Nx_lbm*Ny_lbm*(Nz_lbm-1))).*(k==Nz_lbm & strcmp(boundaryLBM.Z0{met},'P') | k==Nz_lbm & strcmp(boundaryLBM.Zn{met},'P'));
    end
    
    
   
    for a=2:Q
        if a==2
            inext=i;  %Indexes of the neighbour voxel (inext,jnext,knext) [dimensionless].
            jnext=j1;
            knext=k;
        elseif a==3
            inext=i_1;
            jnext=j;
            knext=k;
        elseif a==4
            inext=i;
            jnext=j_1;
            knext=k;
        elseif a==5
            inext=i1;
            jnext=j;
            knext=k;
        elseif a==6
            inext=i;
            jnext=j;
            knext=k1;
        elseif a==7
            inext=i;
            jnext=j;
            knext=k_1;
        end

        if dim==2
            if strcmp(boundaryLBM.Z0{met},'D') || strcmp(boundaryLBM.Zn{met},'D') || strcmp(boundaryLBM.Z0{met},'Dg') || strcmp(boundaryLBM.Zn{met},'Dg')
                knext=-1*k;
            else
                knext=k;
            end
        end
        
        
        

                                                                                                               %Periodic                                                                                          %Bounceback - Neumann                                                       %Dirichlet     
            IDlbm(ycont1,a)=(Nx_lbm*Ny_lbm*Nz_lbm*(met-1)+sub2ind([Ny_lbm,Nx_lbm,Nz_lbm],(inext>0).*inext+(inext<=0),(jnext>0).*jnext+(jnext<=0),(knext>0).*knext+(knext<=0))).*(inext>0 & jnext>0 & knext>0) + (Nx_lbm*Ny_lbm*Nz_lbm*num_met+1).*(inext==0 | jnext==0 | knext==0).*(inext>=0 & jnext>=0 & knext>=0) + (Nx_lbm*Ny_lbm*Nz_lbm*num_met+2).*(inext==-1 | jnext==-1 | knext==-1);

  
        
%         if dim==3
%                                                                                                                %Periodic                                                                                          %Bounceback - Neumann                                                       %Dirichlet     
%             IDlbm(ycont1,a)=(Nx_lbm*Ny_lbm*Nz_lbm*(met-1)+sub2ind([Ny_lbm,Nx_lbm,Nz_lbm],(inext>0).*inext+(inext<=0),(jnext>0).*jnext+(jnext<=0),(knext>0).*knext+(knext<=0))).*(inext>0 & jnext>0 & knext>0) + (Nx_lbm*Ny_lbm*Nz_lbm*num_met+1).*(inext==0 | jnext==0 | knext==0) + (Nx_lbm*Ny_lbm*Nz_lbm*num_met+2).*(inext==-1 | jnext==-1 | knext==-1);
%         else
%                                                                                                                 %Periodic                                                                      %Bounceback  - Neumann                        %Dirichlet     
%             IDlbm(ycont1,a)=(Nx_lbm*Ny_lbm*Nz_lbm*(met-1)+sub2ind([Ny_lbm,Nx_lbm,Nz_lbm],(inext>0).*inext+(inext<=0),(jnext>0).*jnext+(jnext<=0),ones(size(inext)))).*(inext>0 & jnext>0 ) + (Nx_lbm*Ny_lbm*Nz_lbm*num_met+1).*(inext==0 | jnext==0) + (Nx_lbm*Ny_lbm*Nz_lbm*num_met+2).*(inext==-1 | jnext==-1 );
%         end
    end

    
    %Vector of metabolites supply through Dirichlet boundaries 
    idvec=Nx_lbm*Ny_lbm*Nz_lbm*(met-1)+(1:Nx_lbm*Ny_lbm*Nz_lbm)';
    %Identifying the Dirichlet boundaries 
    vecDx=1-((i-1)/(Ny_lbm-1)); %gradient in Y direction
    rhoBDx(idvec)=(CmetDb.X0(met)*vol_lbm).*(j==1 & rhoBDx(idvec)<CmetDb.X0(met)*vol_lbm & strcmp(boundaryLBM.X0{met},'D'))+(CmetDb.Xn(met)*vol_lbm).*(j==Nx_lbm & rhoBDx(idvec)<CmetDb.Xn(met)*vol_lbm & strcmp(boundaryLBM.Xn{met},'D'))+(CmetDb.X0(met)*vol_lbm).*vecDx.*(j==1 & rhoBDx(idvec)<CmetDb.X0(met)*vol_lbm & strcmp(boundaryLBM.X0{met},'Dg'))+(CmetDb.Xn(met)*vol_lbm).*vecDx.*(j==Nx_lbm & rhoBDx(idvec)<CmetDb.Xn(met)*vol_lbm & strcmp(boundaryLBM.Xn{met},'Dg'));
    FDx(idvec)=(j==1 & strcmp(boundaryLBM.X0{met},'D'))+(j==Nx_lbm & strcmp(boundaryLBM.Xn{met},'D'))+(j==1 & strcmp(boundaryLBM.X0{met},'Dg'))+(j==Nx_lbm & strcmp(boundaryLBM.Xn{met},'Dg'));
    
    vecDy=1-((j-1)/(Nx_lbm-1)); %gradient in X direction
    rhoBDy(idvec)=(CmetDb.Y0(met)*vol_lbm).*(i==1 & rhoBDy(idvec)<CmetDb.Y0(met)*vol_lbm & strcmp(boundaryLBM.Y0{met},'D'))+(CmetDb.Yn(met)*vol_lbm).*(i==Ny_lbm & rhoBDy(idvec)<CmetDb.Yn(met)*vol_lbm & strcmp(boundaryLBM.Yn{met},'D'))+(CmetDb.Y0(met)*vol_lbm).*vecDy.*(i==1 & rhoBDy(idvec)<CmetDb.Y0(met)*vol_lbm & strcmp(boundaryLBM.Y0{met},'Dg'))+(CmetDb.Yn(met)*vol_lbm).*vecDy.*(i==Ny_lbm & rhoBDy(idvec)<CmetDb.Yn(met)*vol_lbm & strcmp(boundaryLBM.Yn{met},'Dg'));
    FDy(idvec)=(i==1 & strcmp(boundaryLBM.Y0{met},'D'))+(i==Ny_lbm & strcmp(boundaryLBM.Yn{met},'D'))+(i==1 & strcmp(boundaryLBM.Y0{met},'Dg'))+(i==Ny_lbm & strcmp(boundaryLBM.Yn{met},'Dg'));
  
    if dim==3
        vecDz=1-((k-1)/(Nz_lbm-1)); %gradient in Z direction
        rhoBDz(idvec)=(CmetDb.Z0(met)*vol_lbm).*(k==1 & rhoBDz(idvec)<CmetDb.Z0(met)*vol_lbm & strcmp(boundaryLBM.Z0{met},'D'))+(CmetDb.Zn(met)*vol_lbm).*(k==Nz_lbm & rhoBDz(idvec)<CmetDb.Zn(met)*vol_lbm & strcmp(boundaryLBM.Zn{met},'D'))+(CmetDb.Z0(met)*vol_lbm).*vecDz.*(k==1 & rhoBDz(idvec)<CmetDb.Z0(met)*vol_lbm & strcmp(boundaryLBM.Z0{met},'Dg'))+(CmetDb.Zn(met)*vol_lbm).*vecDz.*(k==Nz_lbm & rhoBDz(idvec)<CmetDb.Zn(met)*vol_lbm & strcmp(boundaryLBM.Zn{met},'Dg'));
        FDz(idvec)=(k==1 & strcmp(boundaryLBM.Z0{met},'D'))+(k==Nz_lbm & strcmp(boundaryLBM.Zn{met},'D'))+(k==1 & strcmp(boundaryLBM.Z0{met},'Dg'))+(k==Nz_lbm & strcmp(boundaryLBM.Zn{met},'Dg'));   
    
    elseif dim==2
        Cz=max([CmetDb.Z0(met),CmetDb.Zn(met)]);
        rhoBDz(idvec)=(Cz*vol_lbm).*(strcmp(boundaryLBM.Z0{met},'D') | strcmp(boundaryLBM.Zn{met},'D') | strcmp(boundaryLBM.Z0{met},'Dg') | strcmp(boundaryLBM.Zn{met},'Dg'));
        FDz(idvec)=(strcmp(boundaryLBM.Z0{met},'D') | strcmp(boundaryLBM.Zn{met},'D') | strcmp(boundaryLBM.Z0{met},'Dg') | strcmp(boundaryLBM.Zn{met},'Dg'));   
    end
    
end

FD=FDx+FDy+FDz;
index=find(FD);
FD(index)=1; %Index vector that ponits (when FD=1) voxels located in the Dirichlet boundary.

FBD=rhoBDx+rhoBDy+rhoBDz; %Vector with the amount molecules that enter in the system through the Dirichlet boundary [mmol].



%Vector ID for the semi-implicity Crank-Nicholson
Morder=zeros(Ny_lbm,Nx_lbm,Nz_lbm,num_met);
Morder(:)=1:Ny_lbm*Nx_lbm*Nz_lbm*num_met;
CNy=Morder(:);
b1y=IDlbm(CNy,3);
b2y=IDlbm(CNy,5);
ap1y=periodic(CNy,3)+(periodic(CNy,3)==0).*(Nx_lbm*Ny_lbm*Nz_lbm*num_met+1);
cp1y=periodic(CNy,5)+(periodic(CNy,5)==0).*(Nx_lbm*Ny_lbm*Nz_lbm*num_met+1);
a1y=zeros(Nx_lbm*Ny_lbm*Nz_lbm*num_met,1);
c1y=zeros(Nx_lbm*Ny_lbm*Nz_lbm*num_met,1);
a1y(:)=2:Nx_lbm*Ny_lbm*Nz_lbm*num_met+1;
c1y(:)=(1:Nx_lbm*Ny_lbm*Nz_lbm*num_met)-1;
c1y(1)=Ny_lbm*Nx_lbm*Nz_lbm*num_met+1;
i=Ny_lbm+1:Ny_lbm:Ny_lbm*Nx_lbm*Nz_lbm*num_met;  
a1y(i-1)=Ny_lbm*Nx_lbm*Nz_lbm*num_met+1;
c1y(i)=Ny_lbm*Nx_lbm*Nz_lbm*num_met+1;
CNx=permute(Morder,[2,1,3,4]);
CNx=CNx(:);
b1x=IDlbm(CNx,4);
b2x=IDlbm(CNx,2);
ap1x=periodic(CNx,4)+(periodic(CNx,4)==0).*(Nx_lbm*Ny_lbm*Nz_lbm*num_met+1);
cp1x=periodic(CNx,2)+(periodic(CNx,2)==0).*(Nx_lbm*Ny_lbm*Nz_lbm*num_met+1);
cn_x=[CNx;Nx_lbm*Ny_lbm*Nz_lbm*num_met+1;Nx_lbm*Ny_lbm*Nz_lbm*num_met+2];
a1x=zeros(Nx_lbm*Ny_lbm*Nz_lbm*num_met,1);
c1x=zeros(Nx_lbm*Ny_lbm*Nz_lbm*num_met,1);
a1x(:)=2:Nx_lbm*Ny_lbm*Nz_lbm*num_met+1;
c1x(:)=(1:Nx_lbm*Ny_lbm*Nz_lbm*num_met)-1;
c1x(1)=Ny_lbm*Nx_lbm*Nz_lbm*num_met+1;
i=Nx_lbm+1:Nx_lbm:Ny_lbm*Nx_lbm*Nz_lbm*num_met;  
a1x(i-1)=Ny_lbm*Nx_lbm*Nz_lbm*num_met+1;
c1x(i)=Ny_lbm*Nx_lbm*Nz_lbm*num_met+1;
a1x=cn_x(a1x);
c1x=cn_x(c1x);
if dim==3
    CNz=permute(Morder,[3,1,2,4]);
    CNz=CNz(:);
    b1z=IDlbm(CNz,7);
    b2z=IDlbm(CNz,6);
    ap1z=periodic(CNz,7)+(periodic(CNz,7)==0).*(Nx_lbm*Ny_lbm*Nz_lbm*num_met+1);
    cp1z=periodic(CNz,6)+(periodic(CNz,6)==0).*(Nx_lbm*Ny_lbm*Nz_lbm*num_met+1);
    cn_z=[CNz;Nx_lbm*Ny_lbm*Nz_lbm*num_met+1;Nx_lbm*Ny_lbm*Nz_lbm*num_met+2];
    a1z=zeros(Nx_lbm*Ny_lbm*Nz_lbm*num_met,1);
    c1z=zeros(Nx_lbm*Ny_lbm*Nz_lbm*num_met,1);
    a1z(:)=2:Nx_lbm*Ny_lbm*Nz_lbm*num_met+1;
    c1z(:)=(1:Nx_lbm*Ny_lbm*Nz_lbm*num_met)-1;
    c1z(1)=Ny_lbm*Nx_lbm*Nz_lbm*num_met+1;
    i=Nz_lbm+1:Nz_lbm:Ny_lbm*Nx_lbm*Nz_lbm*num_met;  
    a1z(i-1)=Ny_lbm*Nx_lbm*Nz_lbm*num_met+1;
    c1z(i)=Ny_lbm*Nx_lbm*Nz_lbm*num_met+1;
    a1z=cn_z(a1z);
    c1z=cn_z(c1z);

else
    b1z=[];
    b2z=[];
    CNz=[];
    ap1z=[];
    cp1z=[];
    a1z=[];
    c1z=[];
end




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

