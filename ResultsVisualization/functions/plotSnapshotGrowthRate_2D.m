function []=plotSnapshotGrowthRate_2D(mplot,nplot,fig,VbiotCA_sp,tvec,h,sp)
%Plot the snapshot of the growth rate of species "sp" at the point "h"

%Color used in the snapshot
cmap3=ones(1000,3);
cmap3(:,2)=[1:-1/1000:1/1000]';
cmap3(:,3)=[1:-1/1000:1/1000]';

%Plot the snapshot
gg=VbiotCA_sp(:,:,:,sp,h);
maxCmet=max(gg(:));
Z=gg;
a4=subplot(mplot,nplot,fig);
imagesc(Z);
colormap(a4,cmap3)
colorbar;
caxis([0 maxCmet]);
dd=strcat('Growth rate sp',num2str(sp), ' [h^-1] ,t = ',num2str(tvec(h)),'h');
title(dd);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
        
hold on
       