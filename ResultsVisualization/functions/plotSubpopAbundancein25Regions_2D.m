function []=plotFracSubpop25Regions_2D(rhotCA_sp,subpoptCA,tvec)

mplot=5;
nplot=5;

%abundance of subpopulations in different regions
[Nx_ca,Ny_ca,Nz_ca,num_sp,~]=size(rhotCA_sp);
step=floor(Ny_ca/5);
fig=1;
for i=1:step:Ny_ca
    for j=1:step:Nx_ca
        rho=rhotCA_sp(i:i+step-1,j:j+step-1,:,:,tt);
        subpop=subpoptCA(i:i+step-1,j:j+step-1,:,sp,tt);
        plotFracSubpop(mplot,nplot,fig,rho,subpop,tvec);
        fig=fig+1;
    end
end

