function []=plotSnapshotPhenotypes_3D_oneLayer(mplot,nplot,fig,phentCA,tvec,h)
%Plot phenotypes in one layer of the 3D biofilm composed by 2 species

%Color used to plot space fraction
cmap4=[1 1 1;0 0 .3;0 0.35 .85;0 0.6 .9;0.6 .95 1];

[Nx_ca,Ny_ca,Nz_ca,num_sp,~]=size(phentCA);

for sp=1:num_sp    
    if sp==1    
        gg=phentCA(floor(Ny_ca/2),:,:,sp,h);
    else 
        gg=gg+(phentCA(floor(Ny_ca/2),:,:,sp,h));
    end
end
maxCmet=max(gg(:));
Z=reshape(gg,[Nx_ca,Nz_ca])';
a4=subplot(mplot,nplot,fig);
imagesc(Z);
colormap(a4,cmap4)
colorbar;
caxis([0 maxCmet]);
dd=strcat('Phenotypes', ',t = ',num2str(tvec(h)),'h');
title(dd);

