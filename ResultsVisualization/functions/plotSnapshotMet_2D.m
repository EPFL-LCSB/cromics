function []=plotSnapshotMet_2D(mplot,nplot,fig,rhotLBM_met,tvec,h,met,cmap2)

if nargin ==7
    %Default Color used in the snapshot
    cmap2=ones(1000,3);
    cmap2(:,1)=(1:-1/1000:1/1000)';
    cmap2(:,2)=(1:-1/1000:1/1000)';
end


%Plot the snapshot
rrh=rhotLBM_met;
uu=rrh(:,:,1,met,h);
ff=max(uu(:));
Z=rrh(:,:,1,met,h); %[mM].
a2=subplot(mplot,nplot,fig);
imagesc(Z);
colormap(a2,cmap2)
colorbar
caxis([0 ff]);
dd=strcat('Met ',num2str(met),' [mmol], t = ',num2str(tvec(h)),'h');
title(dd);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
        
hold on
       