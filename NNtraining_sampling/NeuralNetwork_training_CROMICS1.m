% Training one NN with data transformation


%% E. coli on lactose
clc
clear 

load TestPoints_Ecoli_core_lcts_pTFA_45K.mat

x = (abs(input_tab(:,1:5)))';
sqrt_x = sqrt(x);

temp_t = target_solTable'; %[lcts, o2, meth, ac, gal, bio]
% we separate the galactose consumed or produced: 
t = [-temp_t(1,:); -temp_t(2,:); -temp_t(3,:); temp_t(4,:); -(temp_t(5,:)<0).*temp_t(5,:); (temp_t(5,:)>0).*temp_t(5,:); temp_t(6,:)];


sqrt_t = sqrt(t);

% Choose a Training Function
trainFcn = 'trainbr';
hiddenLayerSize = [15,15]; %For 2 layer of 15
nb_epochs = 10000;
functionName = 'NN15x15_Ecoli_core_lcts_pTFA';

transferFct = 'poslin';

[testPerformance, trainPerformance, effective_param, mus, net] = NN_train(sqrt_x,sqrt_t, trainFcn, hiddenLayerSize, nb_epochs, functionName, transferFct);




%% S. enterica

clc
clear 

load TestPoints_Senterica_pTFA_40K.mat

x = ([abs(input_tab(:,1)),abs(input_tab(:,3:4)),rmeth])'; %[o2, ac, gal]
sqrt_x = sqrt(x);

temp_t = target_solTable'; %[o2, meth, ac, gal, bio]
t = [-temp_t(1,:); temp_t(2,:); -(temp_t(3,:)<0).*temp_t(3,:); (temp_t(3,:)>0).*temp_t(3,:); -temp_t(4,:); temp_t(5,:)];
sqrt_t = sqrt(t);


trainFcn = 'trainbr';
hiddenLayerSize = [15 15];
nb_epochs = 10000;

functionName = 'NN15x15_Senterica_pTFA';
transferFct = 'poslin';


[testPerformance, trainPerformance, effective_param, mus, net] = NN_train(sqrt_x,sqrt_t, trainFcn, hiddenLayerSize, nb_epochs, functionName, transferFct);



%% E. coli on glucose
clc
clear 

load TestPoints_Ecoli_glc_o2_ac_pTFA_30K.mat

x = (abs(input_tab(:,1:5)))';
sqrt_x = sqrt(x);

temp_t = target_solTable'; %[lcts, o2, meth, ac, gal, bio]
% we separate the galactose consumed or produced: 
t = [-temp_t(1,:); -temp_t(2,:); -temp_t(3,:); temp_t(4,:); -(temp_t(5,:)<0).*temp_t(5,:); (temp_t(5,:)>0).*temp_t(5,:); temp_t(6,:)];


sqrt_t = sqrt(t);

% Choose a Training Function
trainFcn = 'trainbr';
hiddenLayerSize = [15,15]; %For 2 layer of 15
nb_epochs = 10000;
functionName = 'NN15x15_Ecoli_glc_pTFA';

transferFct = 'poslin';

[testPerformance, trainPerformance, effective_param, mus, net] = NN_train(sqrt_x,sqrt_t, trainFcn, hiddenLayerSize, nb_epochs, functionName, transferFct);



