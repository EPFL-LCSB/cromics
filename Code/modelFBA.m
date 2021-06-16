function [flux,flux_bio]=modelFBA(vshrinkage,limflux,Ncells,num_metInert,sp,GEMs,~)

    GEM=eval(strcat('GEMs.sp',num2str(sp)));
    rxn=find(GEM.indexMetF);

    RxnF=[GEM.indexMetF(rxn);GEM.indexBioF];
  
    flux=zeros(Ncells,num_metInert+size(limflux,2));
    flux_bio=zeros(Ncells,1);


    parfor ii=1:Ncells

        model=GEM.model;

        %Modify the upper limits for exchange fluxes
        model.lb(RxnF(1:end-1))=-limflux(ii,rxn); 

        %Max biomass 
        model.f(:)=0;
        model.f(RxnF(end))=1;
        flux1=FBA_cplex(model,RxnF);

        flux(ii,rxn)=flux1(1:end-1);
        flux_bio(ii)=flux1(end)+(flux1(end)<=0)*(-vshrinkage);

    end

    flux=[flux,zeros(Ncells,num_metInert)]; 

end




function [flux]=FBA_cplex(model,RxnF)

    LP.A=model.S;
    LP.b=model.b;
    LP.lb=model.lb;
    LP.ub=model.ub;
    LP.c=model.c;
    LP.osense=-1; %Maximization;
    LP.csense(size(model.mets))='E';
    sol=solveCobraLP(LP);

    flux=zeros(1,length(RxnF));
    if isempty(sol.val)==0 && isnan(sol.val)==0
        flux=sol.x(RxnF); 
    end
end
