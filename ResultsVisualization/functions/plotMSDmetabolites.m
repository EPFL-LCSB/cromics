function []=plotMSDmetabolites(mplot,nplot,fig,r2tt_met,num_met,tvec)
%Plot the mean square displacement of the metabolites  

subplot(mplot,nplot,fig);
line={'b','r','g','k','m','y','c','--b','--r','--g','--k','--m','--y','--c'};

for met=1:num_met
    plot(tvec,r2tt_met(met,:)',line{met});
    ff{met}=strcat('met',num2str(met));
    hold on
end
title('MSDmetabolites vs time');
xlabel('Time [h]');
ylabel('MSD [mm^2]');
legend (ff{:});
