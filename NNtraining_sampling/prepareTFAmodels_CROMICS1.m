%CROMICS 1

%E. coli core mutant for: lcts, ac, o2, meth
clc
clear 
add_cplex;
load GSmodel_Ecoli.mat 
model=ttmodel;
model.ub(ismember(model.rxns,'SHSL1'))=0; %Block reaction O-succinylhomoserine lyase (L-cysteine), i.e. cystathionine gamma-synthase
model.lb(ismember(model.rxns,'SHSL1'))=0;
model.lb(ismember(model.rxns,'ATPM'))=3.15; 


lcts='DM_lcts_e';
o2='DM_o2_e';
ac='DM_ac_e';
meth='DM_met-L_e';
gal='DM_gal_e';
bio='Ec_biomass_iJO1366_core_53p95M';
glc='DM_glc_e';

model.lb(ismember(model.rxns,glc))=0; %glc exchange
model.ub(ismember(model.rxns,glc))=0; %glc exchange
model.lb(ismember(model.rxns,bio))=0;
model.ub(ismember(model.rxns,bio))=1000;
model.lb(ismember(model.rxns,'Ec_biomass_iJO1366_WT_53p95M'))=0;
model.ub(ismember(model.rxns,'Ec_biomass_iJO1366_WT_53p95M'))=0;
model.lb(ismember(model.rxns,ac))=-10;
model.ub(ismember(model.rxns,ac))=1000;
model.lb(ismember(model.rxns,meth))=-10;
model.ub(ismember(model.rxns,meth))=0;
model.lb(ismember(model.rxns,o2))=-10;
model.ub(ismember(model.rxns,o2))=0;
model.lb(ismember(model.rxns,gal))=-10;
model.ub(ismember(model.rxns,gal))=1000;
model.lb(ismember(model.rxns,lcts))=-10;
model.ub(ismember(model.rxns,lcts))=0;

model.c(:)=0;
model.c(ismember(model.rxns,bio))=1;

load DB_AlbertyUpdate_keng.mat

model.CompartmentData.compMinConc(:)=1e-5;
model.CompartmentData.compMinConc(5)=1e-8;
model.CompartmentData.compMaxConc(:)=0.02;
model.CompartmentData.compMaxConc(5)=0.1;


% lcts + gal + o2 --> ac + biomass
ld=10;
ad=0;
od=10; 
md=1;
gl=10;

model.lb(ismember(model.rxns,lcts))=-ld;
model.ub(ismember(model.rxns,lcts))=0;
model.lb(ismember(model.rxns,ac))=-ad;
model.ub(ismember(model.rxns,ac))=1000;
model.lb(ismember(model.rxns,o2))=-od;
model.ub(ismember(model.rxns,o2))=0;
model.lb(ismember(model.rxns,meth))=-md;
model.ub(ismember(model.rxns,meth))=0;
model.lb(ismember(model.rxns,gal))=-gl;
model.ub(ismember(model.rxns,gal))=1000;

tmodel = prepModelforTFBA(model,DB_AlbertyUpdate,model.CompartmentData);
[ttmodel] = convToTFA(tmodel, DB_AlbertyUpdate, {}, 'DGo', [], 0.07, 1, 1); 
lb1=ttmodel.var_lb;
ub1=ttmodel.var_ub;


% lcts + gal --> ac + biomass
ld=10;
ad=0;
od=0; 
md=1;
gl=10;

model.lb(ismember(model.rxns,lcts))=-ld;
model.ub(ismember(model.rxns,lcts))=0;
model.lb(ismember(model.rxns,ac))=-ad;
model.ub(ismember(model.rxns,ac))=1000;
model.lb(ismember(model.rxns,o2))=-od;
model.ub(ismember(model.rxns,o2))=0;
model.lb(ismember(model.rxns,meth))=-md;
model.ub(ismember(model.rxns,meth))=0;
model.lb(ismember(model.rxns,gal))=-gl;
model.ub(ismember(model.rxns,gal))=1000;

tmodel = prepModelforTFBA(model,DB_AlbertyUpdate,model.CompartmentData);
[ttmodel] = convToTFA(tmodel, DB_AlbertyUpdate, {}, 'DGo', [], 0.07, 1, 1); 
lb2=ttmodel.var_lb;
ub2=ttmodel.var_ub;


% ac + gal + o2 --> biomass
ld=0;
ad=10;
od=10;  
md=1;
gl=10;

model.lb(ismember(model.rxns,lcts))=-ld;
model.ub(ismember(model.rxns,lcts))=0;
model.lb(ismember(model.rxns,ac))=-ad;
model.ub(ismember(model.rxns,ac))=1000;
model.lb(ismember(model.rxns,o2))=-od;
model.ub(ismember(model.rxns,o2))=0;
model.lb(ismember(model.rxns,meth))=-md;
model.ub(ismember(model.rxns,meth))=0;
model.lb(ismember(model.rxns,gal))=-gl;
model.ub(ismember(model.rxns,gal))=1000;

tmodel = prepModelforTFBA(model,DB_AlbertyUpdate,model.CompartmentData);
[ttmodel] = convToTFA(tmodel, DB_AlbertyUpdate, {}, 'DGo', [], 0.07, 1, 1); 
lb3=ttmodel.var_lb;
ub3=ttmodel.var_ub;

%compare the limits:
ttmodel.var_lb=min([lb1,lb2,lb3],[],2);
ttmodel.var_ub=max([ub1,ub2,ub3],[],2);

save ('Ecoli_core_lcts_TFA.mat','ttmodel')


%%
%S. enterica mutant for: lcts, ac, o2, meth, gal
clc
clear 

%****** change "EX_" for "DM_"
load GEM_Senterica.mat

nrxn=length(ttmodel.rxns); 
for i=1:nrxn
    aa=strfind(ttmodel.rxns(i),'EX_');
    if isempty(aa{1})==0
        ttmodel.rxns(i)=strrep(ttmodel.rxns(i),'EX_','DM_');
        ttmodel.rxns(i)=strrep(ttmodel.rxns(i),'(e)','_e');
    end
end
%************************
model=ttmodel;

gal='DM_gal_e';
o2='DM_o2_e';
ac='DM_ac_e';
meth='DM_met-L_e';
bio='biomass_iRR1083_metals';
glc='DM_glc-D_e';
lcts='DM_lcts_e';
model.lb(ismember(model.rxns,glc))=0; %glc exchange
model.ub(ismember(model.rxns,glc))=0; 
model.lb(ismember(model.rxns,lcts))=0; %lcts exchange
model.ub(ismember(model.rxns,lcts))=0; 
model.lb(ismember(model.rxns,bio))=0;
model.ub(ismember(model.rxns,bio))=1000;
model.lb(ismember(model.rxns,ac))=-10;
model.lb(ismember(model.rxns,gal))=-10;

load db.mat 
load Ecoli_CompartmentData.mat
model = mapMetIDs(model,DB);
model.metSEEDID=model.metSEEDID';
model.CompartmentData=Ecoli_CompartmentData;
model = getMetCompartment(model,Ecoli_CompartmentData,'[]');
load DB_AlbertyUpdate_keng.mat

model.CompartmentData.compMinConc(:)=1e-5;
model.CompartmentData.compMinConc(5)=1e-8;
model.CompartmentData.compMaxConc(:)=0.02;
model.CompartmentData.compMaxConc(5)=0.1;

% ac + o2 --> meth + biomass
ad=10;
od=10;  

model.lb(ismember(model.rxns,ac))=-ad;
model.ub(ismember(model.rxns,ac))=0;
model.lb(ismember(model.rxns,o2))=-od;
model.ub(ismember(model.rxns,o2))=0;
model.lb(ismember(model.rxns,meth))=0;
model.ub(ismember(model.rxns,meth))=1000;

tmodel = prepModelforTFBA(model,DB_AlbertyUpdate,model.CompartmentData);
[ttmodel] = convToTFA(tmodel, DB_AlbertyUpdate, {}, 'DGo', [], 0.07, 1, 1); 
lb1=ttmodel.var_lb;
ub1=ttmodel.var_ub;


% gal + o2 --> ac + biomass
ld=10;
ad=0;
od=10;  
md=1;

model.lb(ismember(model.rxns,gal))=-ld;
model.ub(ismember(model.rxns,gal))=0;
model.lb(ismember(model.rxns,ac))=-ad;
model.ub(ismember(model.rxns,ac))=1000;
model.lb(ismember(model.rxns,o2))=-od;
model.ub(ismember(model.rxns,o2))=0;
model.lb(ismember(model.rxns,meth))=0;
model.ub(ismember(model.rxns,meth))=1000;

tmodel = prepModelforTFBA(model,DB_AlbertyUpdate,model.CompartmentData);
[ttmodel] = convToTFA(tmodel, DB_AlbertyUpdate, {}, 'DGo', [], 0.07, 1, 1); 
lb2=ttmodel.var_lb;
ub2=ttmodel.var_ub;


% gal --> ac + biomass
ld=10;
ad=0;
od=0; 
md=1;

model.lb(ismember(model.rxns,gal))=-ld;
model.ub(ismember(model.rxns,gal))=0;
model.lb(ismember(model.rxns,ac))=-ad;
model.ub(ismember(model.rxns,ac))=1000;
model.lb(ismember(model.rxns,o2))=-od;
model.ub(ismember(model.rxns,o2))=0;
model.lb(ismember(model.rxns,meth))=0;
model.ub(ismember(model.rxns,meth))=1000;

tmodel = prepModelforTFBA(model,DB_AlbertyUpdate,model.CompartmentData);
[ttmodel] = convToTFA(tmodel, DB_AlbertyUpdate, {}, 'DGo', [], 0.07, 1, 1); 
lb3=ttmodel.var_lb;
ub3=ttmodel.var_ub;


%compare the limits:
ttmodel.var_lb=min([lb1,lb2,lb3],[],2);
ttmodel.var_ub=max([ub1,ub2,ub3],[],2);

save ('Senterica_meth_TFA.mat','ttmodel')


%% E. coli WT: glc, ac, o2

clc
clear 
load GSmodel_Ecoli.mat 
model=ttmodel;

lcts='DM_lcts_e';
glc='DM_glc_e';
o2='DM_o2_e';
ac='DM_ac_e';
bio='Ec_biomass_iJO1366_WT_53p95M';
fe='DM_fe3_e';
model.lb(ismember(model.rxns,'ATPM'))=3.15; 
model.lb(ismember(model.rxns,fe))=0; %Fe3+ exchange
model.lb(ismember(model.rxns,bio))=0;
model.ub(ismember(model.rxns,bio))=1000;

load DB_AlbertyUpdate_keng.mat

% glc + o2 --> ac + biomass
gd=10;
ad=0;
od=15; 

model.lb(ismember(model.rxns,glc))=-gd;
model.ub(ismember(model.rxns,glc))=0;
model.lb(ismember(model.rxns,ac))=-ad;
model.ub(ismember(model.rxns,ac))=1000;
model.lb(ismember(model.rxns,o2))=-od;
model.ub(ismember(model.rxns,o2))=0;

tmodel = prepModelforTFBA(model,DB_AlbertyUpdate,model.CompartmentData);
[ttmodel] = convToTFA(tmodel, DB_AlbertyUpdate, {}, 'DGo', [], 0.07, 1, 1); 
lb1=ttmodel.var_lb;
ub1=ttmodel.var_ub;


% glc --> ac + biomass
gd=10;
ad=0;
od=0; 

model.lb(ismember(model.rxns,glc))=-gd;
model.ub(ismember(model.rxns,glc))=0;
model.lb(ismember(model.rxns,ac))=-ad;
model.ub(ismember(model.rxns,ac))=1000;
model.lb(ismember(model.rxns,o2))=-od;
model.ub(ismember(model.rxns,o2))=0;

tmodel = prepModelforTFBA(model,DB_AlbertyUpdate,model.CompartmentData);
[ttmodel] = convToTFA(tmodel, DB_AlbertyUpdate, {}, 'DGo', [], 0.07, 1, 1); 
lb2=ttmodel.var_lb;
ub2=ttmodel.var_ub;


% ac + o2 --> biomass
gd=0;
ad=17;
od=15; 

model.lb(ismember(model.rxns,glc))=-gd;
model.ub(ismember(model.rxns,glc))=0;
model.lb(ismember(model.rxns,ac))=-ad;
model.ub(ismember(model.rxns,ac))=1000;
model.lb(ismember(model.rxns,o2))=-od;
model.ub(ismember(model.rxns,o2))=0;

tmodel = prepModelforTFBA(model,DB_AlbertyUpdate,model.CompartmentData);
[ttmodel] = convToTFA(tmodel, DB_AlbertyUpdate, {}, 'DGo', [], 0.07, 1, 1); 
lb3=ttmodel.var_lb;
ub3=ttmodel.var_ub;


%compare the limits:
ttmodel.var_lb=min([lb1,lb2,lb3],[],2);
ttmodel.var_ub=max([ub1,ub2,ub3],[],2);

model=ttmodel;

save ('Ecoli_glc_ac_o2_TFA.mat','model')



