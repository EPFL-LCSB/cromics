function []=plotSpeciesRelAbundance(mplot,nplot,fig,num_sp,tvec,rhotCA_sp)

subplot(mplot,nplot,fig);
line={'b','r','g','k','m','y','c','--b','--r','--g','--k','--m','--y','--c'};
SpeciesRatio=zeros(length(tvec),num_sp); %Species abundance ratio
  
for sp=1:num_sp
    for tt=1:length(tvec)
        SpeciesRatio(tt,sp)=sum(sum(sum(rhotCA_sp(:,:,:,sp,tt))))/sum(sum(sum(sum(rhotCA_sp(:,:,:,:,tt))))); %Dimensionless
    end
    plot(tvec,SpeciesRatio(:,sp),line{sp});
    ff{sp}=strcat('sp',num2str(sp));
    hold on
end

title('Species relative abundance');
xlabel('Time [h]');
ylabel('Relative abundance [dimensionless]')
legend (ff{:});
