function []=plotFracSubpop(mplot,nplot,fig,rhotCA_sp,subpoptCA,tvec)
%Plot the evolution of the relative abundance of the populations and
%subpopulations

subplot(mplot,nplot,fig);
line={'b','r','g','k','m','y','c'};

[Nx_ca,Ny_ca,Nz_ca,num_sp,~]=size(rhotCA_sp);
frac_sp=zeros(length(tvec),num_sp);
pop=unique(subpoptCA);
frac_spop=zeros(length(tvec),num_sp,length(pop));

mass_tot=sum(reshape(rhotCA_sp(:,:,:,:,:),[Nx_ca*Ny_ca*Nz_ca*num_sp,length(tvec)]))';
pp=1;
for sp=1:num_sp
    for tt=1:length(tvec)
        rho=rhotCA_sp(:,:,:,sp,tt);
        frac_sp(tt,sp)=sum(rho(:))/mass_tot(tt);
        for spop=1:length(pop)-1
            aa=find(subpoptCA(:,:,:,sp,tt)==pop(spop+1));
            if ~isempty(aa)
                frac_spop(tt,sp,spop)=sum(rho(aa))/mass_tot(tt);
            end
            if tt==1 && ~isempty(find(subpoptCA(:,:,:,sp,:)==pop(spop+1)))
                ff{pp}=strcat('sp',num2str(sp),' subpop',num2str(pop(spop+1)));
                pp=pp+1;
            end
        end
        
    end
    
end
pp=1;
for sp=1:num_sp
    for spop=1:length(pop)-1
        if ~isempty(find(frac_spop(:,sp,spop)))
            plot(tvec,frac_spop(:,sp,spop),line{pp});
            pp=pp+1;
            hold on
        end
    end
end
set(gca,'visible','on');
ylim([0 1])
xlim([0 tvec(end)])
legend (ff{:});
hold on
