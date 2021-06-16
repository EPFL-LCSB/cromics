function []=plotSnapshotSubpop_2D(mplot,nplot,fig,subpoptCA,tvec,h,sp)

%Color used in the snapshot
cmap4=[1 1 1;0 1 0;1 0 1;0 0 1;0 0 0;0 1 1;1 0 1;1 1 0;0 .5 .7;.6 0 .5;.2 .9 0];

su=subpoptCA(:,:,:,sp,h);
gg=su;
maxCmet=max(gg(:));
Z=su;
a4=subplot(mplot,nplot,fig);
imagesc(Z);
colormap(a4,cmap4)
colorbar;
caxis([0 maxCmet]);
hold on;
dd=strcat('Snapshot subpopulations', ',t = ',num2str(tvec(h)),'h');
title(dd);
