function [Vmt,flux_bio,flux]=MetabolicFluxes(Vmt,At_lbm,Ncells,rhoLBM_met,VavLBM,lbmgrid,maxMcell,IDpos,n,Vmax,Km,VmaxDiff,vshrinkage,num_metInert,rtot,NNlist,num_sp,cellDist,posBio,ia,ic,GEMs,varargin)
    %Compute the metabolic fluxes of each cell

    
    %Compute the affective concentration of the metabolites [mM].
    CmetX=rhoLBM_met*(10^6)./VavLBM; 
    
    %If cell biomass reaches the maximum cell mass (maxMcell), then cell is not allowed 
    %to grow or consume nutrients, therefore it not considered in the metabolic flux step   
    biomass=(IDpos(1:Ncells,3)<=maxMcell).*IDpos(1:Ncells,3); %Biomass of each cell [gDW].

    %If there is no cell distiction (i.e. cellDist=0), then cells of the same species 
    %under the same environmental concentration behave in the same way.
    if cellDist==0 && n>1 
        bbio=accumarray(posBio,biomass);
        biomass=bbio(posBio(ia));
    end
    
    lim1=Vmax(IDpos(ia,2),:).*CmetX(IDpos(ia,4),:)./(CmetX(IDpos(ia,4),:)+Km(IDpos(ia,2),:));
    lim1(isnan(lim1))=0;
    lim1=lim1+VmaxDiff(IDpos(ia,2),:);

    if n>1 %When a LBM-voxel can allocate more than 1 CA-voxel
        rhomet1=lim1.*biomass*(At_lbm/(1000*3600));
        %Total amount of metabolites that would be consumed if cells in a LBM-voxel 
        %work maximum physiological uptake flux 
        rtot(1:max(lbmgrid(:)))=accumarray(lbmgrid(:),rhomet1(:)); 
        limflux=-min(lim1,(rhoLBM_met(IDpos(ia,4),:).*lim1./rtot(lbmgrid))).*(biomass>0);  
    else %When the size of 1 LBM-voxel is equal to CA-voxel
        limflux=-min(lim1,(rhoLBM_met(IDpos(ia,4),:)./(biomass*(At_lbm/(1000*3600))))).*(biomass>0);  
    end

      
    %Compute metabolic fluxes using Neural Networks
    [flux,flux_bio]=MetFlux_NN(NNlist,vshrinkage,limflux,num_metInert,num_sp,IDpos(ia,:),GEMs);

    if n>1
        a=flux.*biomass/(1000*3600);
        Vmt(1:max(lbmgrid(:)))=accumarray(lbmgrid(:),a(:));
        
        if cellDist==0 
            flux=flux(ic,:);
            flux_bio=flux_bio(ic);
        end 
    else
        Vmt(lbmgrid)=flux.*biomass/(1000*3600);
    end
end



%Compute the metabolic fluxes using the functions in NNlist
function [flux,flux_bio]=MetFlux_NN(NNlist,vshrinkage,limflux,num_metInert,num_sp,IDpos,GEMs)
         
    flux=limflux(:,1:end);
    flux_bio=limflux(:,end);
    
    for sp=1:num_sp
        funstr=str2func(NNlist{sp});
        vec=find(IDpos(:,2)==sp);
        if size(IDpos,2)>5 %IDpos(:,6) identifies the subpopulations
            [flux(vec,:),flux_bio(vec)]=funstr(vshrinkage,limflux(vec,:),length(vec),num_metInert,sp,GEMs,IDpos(vec,6));
        else
            [flux(vec,:),flux_bio(vec)]=funstr(vshrinkage,limflux(vec,:),length(vec),num_metInert,sp,GEMs);
        end
    end
end
