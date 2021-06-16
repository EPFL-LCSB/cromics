function []=plotBiomass(mplot,nplot,fig,num_sp,tvec,rhotCA_sp)

subplot(mplot,nplot,fig);
line={'b','r','g','k','m','y','c','--b','--r','--g','--k','--m','--y','--c'};
Csp1=zeros(length(tvec),num_sp);    %Total species biomass [gDW] at different time.
   
for sp=1:num_sp
   for tt=1:length(tvec)
         Csp1(tt,sp)=sum(sum(sum(rhotCA_sp(:,:,:,sp,tt)))); %[gDW].
    end
    plot(tvec,Csp1(:,sp),line{sp});
    ff{sp}=strcat('sp',num2str(sp));
    hold on
end

title('Biomass');
xlabel('Time [h]');
ylabel('Biomass [gDW]')
legend (ff{:});
