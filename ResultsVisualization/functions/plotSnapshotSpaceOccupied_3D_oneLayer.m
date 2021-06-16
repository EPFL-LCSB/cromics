function []=plotSnapshotSpaceOccupied_3D_oneLayer(mplot,nplot,fig,rhotCA_sp,tvec,v_sp,vol_ca,h)
%Plot the space fraction occupied by cells and macrolecules in one layer of the 3D biofilm composed by 2 species


%Color used to plot space fraction
cmap2=ones(1000,3);
cmap2(:,1)=[1:-1/1000:1/1000]';
cmap2(:,2)=[1:-1/1000:1/1000]';
cmap2(:,3)=[1:-1/1000:1/1000]';

[Nx_ca,Ny_ca,Nz_ca,num_sp,~]=size(rhotCA_sp);

for sp=1:num_sp    
    if sp==1
        gg=rhotCA_sp(floor(Ny_ca/2),:,:,1,h).*v_sp(sp)/vol_ca;
    else
        gg=gg+rhotCA_sp(floor(Ny_ca/2),:,:,2,h).*v_sp(sp)/vol_ca;
    end
end
maxCmet=max(gg(:));
Z=reshape(gg,[Nx_ca,Nz_ca])';
a4=subplot(mplot,nplot,fig);
imagesc(Z);
colormap(a4,cmap2)
colorbar;
caxis([0 maxCmet]);
dd=strcat('Space occupied', ',t = ',num2str(tvec(h)),'h');
title(dd);

