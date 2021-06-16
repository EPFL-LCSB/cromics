function []=plotSnapshotPhen_2D(mplot,nplot,fig,phentCA,tvec,h,sp)

%Default color used in the snapshot
cmap4=[1 1 1;0 1 0;1 0 1;0 0 1;0 0 0;0 1 1;1 0 1;1 1 0;0 .5 .7;.6 0 .5;.2 .9 0];


prh=phentCA;
gg=prh(:,:,:,sp,:);
maxCmet=max(gg(:));
Z=prh(:,:,1,sp,h);
a4=subplot(mplot,nplot,fig);
imagesc(Z);
colormap(a4,cmap4)
colorbar;
caxis([0 maxCmet]);
dd=strcat('Phenotypes of sp', num2str(sp),', t = ',num2str(tvec(h)),'h');
title(dd);
set(gca,'xtick',[]);
set(gca,'ytick',[]);


            
