%CROMICS 2

%E. coli for: glc, ac, o2. 

clc
clear 
add_cplex;
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


%*******  aTFA  *******
%Add the proteome allocation constraint based on CAFBA (Mori et. al, 2016).

%Name of the reactions whose metabolic flux will be modified
glc='DM_glc_e';
model=setUpModel_aTFA(model,'glc'); %If the weight assigned for the protein group C (transport) is equal to zero then the substrate "subst" selected is not important  

%*********************


save ('Ecoli_glc_ac_o2_aTFA.mat','model')