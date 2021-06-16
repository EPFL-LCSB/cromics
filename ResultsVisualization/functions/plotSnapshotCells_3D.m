function []=plotSnapshotCells_3D(mplot,nplot,fig,postCA_sp,tvec,h,num_sp)
%Plot a 3D biofilm 

subplot(mplot,nplot,fig);
%Color used to plot cells
color1=[0 0 1;1 0 0;0 1 0;0 0 0;1 0 1;1 1 0;0 1 1];
color11=[0 0 0.7;0.8 0 0;0 0.8 0;0.1 0.1 0.1;0.8 0 0.8;0.8 0.8 0;0 0.8 0.8];

[Nx_ca,Ny_ca,Nz_ca,~]=size(postCA_sp);

for sp=1:num_sp
    gg=find(postCA_sp(:,:,:,h)==sp);
    [jj,ii,kk]=ind2sub([Ny_ca,Nx_ca,Nz_ca],gg);
    kk2=Nz_ca-kk;
    plot3(ii,jj,kk2,'o','MarkerSize',5,'MarkerEdgeColor',color1(sp,:),'MarkerFaceColor',color11(sp,:));
    axis([0 Ny_ca 0 Nx_ca 0 Nz_ca])
    hold on
end

dd=strcat('Biofilm', ',t = ',num2str(tvec(h)),'h');
title(dd);