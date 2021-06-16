function [flux,flux_bio]=modelTFA(vshrinkage,limflux,Ncells,num_metInert,sp,GEMs,~)

    GEM=eval(strcat('GEMs.sp',num2str(sp)));
    rxn=find(GEM.indexMetF);

    RxnF=[GEM.indexMetF(rxn);GEM.indexBioF];
    RxnR=[GEM.indexMetR(rxn);GEM.indexBioR];

    flux=zeros(Ncells,num_metInert+size(limflux,2));
    flux_bio=zeros(Ncells,1);


    parfor ii=1:Ncells

        model=GEM.model;

        %Modify the upper limits for exchange fluxes
        model.var_ub(RxnR(1:end-1))=-limflux(ii,rxn); 
        model.var_ub(RxnF(1:end-1))=1000; 


    %****************** Max biomass ***********************
        model.f(:)=0;
        model.f(RxnF(end))=1;
        flux1=solveProblem(model,RxnF,RxnR);

    % %****************** parsimonious TFA ***********************
        model.var_lb(RxnF(end))=flux1(end); 
        model.var_ub(RxnF(end))=flux1(end); 
        %Minimize the sum of fluxes 
        model.f(1:2*size(model.rxns,2))=-1;
        [flux1]=solveProblem(model,RxnF,RxnR);

         flux(ii,rxn)=flux1(1:end-1);
         flux_bio(ii)=flux1(end)+(flux1(end)<=0)*(-vshrinkage);

    end

    flux=[flux,zeros(Ncells,num_metInert)]; 

end




function [flux]=solveProblem(model,RxnF,RxnR)
    sol=solveTFAmodelCplex(model);

    flux=zeros(1,length(RxnF));
    if isempty(sol.val)==0 && isnan(sol.val)==0
        flux=sol.x(RxnF)-sol.x(RxnR); 
    end
end
