function []=plotMetAbundance(mplot,nplot,fig,num_met,tvec,rhotLBM_met)

subplot(mplot,nplot,fig);
line={'b','r','g','k','m','y','c','--b','--r','--g','--k','--m','--y','--c'};

Cm1=zeros(length(tvec),num_met);    %Total amount of metabolites [mmol] at different time.
for met=1:num_met
   for tt=1:length(tvec)
         Cm1(tt,met)=sum(sum(sum(rhotLBM_met(:,:,:,met,tt)))); %[mmol].
    end
    plot(tvec,Cm1(:,met),line{met});
    ff{met}=strcat('met',num2str(met));
    hold on
end

title('Metabolites abundance');
xlabel('Time [h]');
ylabel('Metabolites abundance [mmol]')
legend (ff{:});