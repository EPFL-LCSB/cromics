function []=plotPhenotypes(mplot,nplot,fig,sp,tvec,phentCA)
%Plot the number of phenotypes of species "sp"

subplot(mplot,nplot,fig);
line={'b','r','g','k','m','y','c','--b','--r','--g','--k','--m','--y','--c'};

Nph1=zeros(length(tvec),max(unique(phentCA(:,:,:,sp,:))));         %Total number of cells with phenotypes phen [individuals] at different time.
pp=1;
for ph=1:max(unique(phentCA(:,:,:,sp,:)))
    if ~isempty(find(phentCA(:,:,:,sp,:)==ph))
        for tt=1:length(tvec)
            Nph1(tt,ph)=length(find(phentCA(:,:,:,sp,tt)==ph)); %[individuals].
        end
        plot(tvec,Nph1(:,ph),line{ph});
        ff{pp}=strcat('Phen',num2str(ph));
        hold on
        pp=pp+1;
    end
end

title(strcat('Phenotypes of species',num2str(sp)));
xlabel('Time [h]');
ylabel('Phenotypes [individuals]')
legend (ff{:});


