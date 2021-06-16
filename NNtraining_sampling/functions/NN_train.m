function [testPerformance, trainPerformance, trainTime, effective_param, mus, net] = NN_train(x,t, trainFcn, hiddenLayerSize, nb_epochs, functionName, transferFct)
    
% This function trains a neural network with MATLAB neural network toolbox,
% based on given samples of input-outputs (x and t).
%
%  
% INPUTS :
%
% x (3 x nb_samples) : contains the input fluxes (glc, o2, ac)
%
% t (4 x nb_samples) : contains the targets (glc, o2, ac, bio), i.e.
% their true values given by FBA
%
% trainFcn : training function
% For a list of all training functions type: help nntrain
% 'trainlm' is usually fastest.
% 'trainbr' takes longer but may be better for challenging problems.
% 'trainscg' uses less memory. NFTOOL falls back to this in low memory situations.
% NB : trainbr seems to be the best training function for our task.
%
% hiddenLayerSize : architecture of the trained neural network
% e.g. 30 for 30 neurons, [10,10] for two layers of 10 neurons, etc.
%
% nb_epochs : maximal number of epochs (iterations of the training
% algorithm)
%
% functionName : desired name of the MATLAB function that will be
% created to re-use your neural network later for prediction.
%
% transferFct : transfer function of the last layer
%   'purelin' : linear function (y=x). Used by default.
%   'tansig' : hyperbolic tangent (tanh) function. To use when targets can
%   be negative. 
%   'poslin' : ReLU function. To use when targets are only positive.
%   'logsig' : logistic function. To use when targets are bounded between 0
%   and 1.
%
%
% OUTPUTS :
%
% testPerformance : standardized MSE obtained on the test set
% trainPerformance : standardized MSE obtained on the training set
% trainTime : time taken to train the neural network
%
% NB : this function automatically creates another function in the same
% folder, containing the trained neural network.
    
    

    % Create a Fitting Network
    net = fitnet(hiddenLayerSize,trainFcn);
    
    % Changing the transfer function of the last (output) layer
    net.layers{length(net.layers)}.transferFcn = transferFct;

    % Choose Input and Output Pre/Post-Processing Functions
    % For a list of all processing functions type: help nnprocess
    % (no pre/post-processing here)
    net.input.processFcns = {};
    net.output.processFcns = {};

    % Setup Division of Data for Training, Validation, Testing
    % For a list of all data division functions type: help nndivide
    % NB : with 'trainbr' as training function, the algorithm cannot accept
    % any testing set (i.e. the parameter 'testRatio' will be ignored).
    net.divideFcn = 'dividerand';  % Divide data randomly
    net.divideMode = 'sample';  % Divide up every sample
    net.divideParam.trainRatio = 70/100;
    net.divideParam.valRatio = 15/100;
    net.divideParam.testRatio = 15/100;
    
    % Choose a Performance Function
    % For a list of all performance functions type: help nnperformance
    net.performFcn = 'mse';  % Mean squared error
    % normalization of outputs between -1 and 1 
    net.performParam.normalization = 'standard'; 

    % Minimal value of the gradient. If reached, the training will be
    % stopped.
    net.trainParam.min_grad = 10^(-9);
    % Maximal number of epochs (iterations). If reached,the training will
    % be stopped.
    net.trainParam.epochs = nb_epochs;

    % Choose Plot Functions
    % For a list of all plot functions type: help nnplot
    net.plotFcns = {'plotperform','plottrainstate','ploterrhist', ...
      'plotregression', 'plotfit'};

    % Train the Network
    [net,tr] = train(net,x,t);
    
    % Saving the training output parameters
    if(strcmp(trainFcn,'trainbr'))
        effective_param = tr.gamk;
    else
        effective_param = [];
    end
    mus = tr.mu; 

    % Test the Network on the whole dataset
    y = net(x);

    % Computing the performance
    e = gsubtract(t,y);
    performance = perform(net,t,y);

    % Recalculate Training, Validation and Test Performance
    trainTargets = t .* tr.trainMask{1};
    valTargets = t  .* tr.valMask{1};
    testTargets = t  .* tr.testMask{1};
    trainPerformance = perform(net,trainTargets,y);
    valPerformance = perform(net,valTargets,y);
    testPerformance = perform(net,testTargets,y);
    trainTime = tr.time(end);


    % Deployment
    % Change the (false) values to (true) to enable the following code blocks.
    if (false)
      % Generate MATLAB function for neural network for application deployment
      % in MATLAB scripts or with MATLAB Compiler and Builder tools, or simply
      % to examine the calculations your trained neural network performs.
      genFunction(net,functionName);
    end
    if (true)
      % Generate a matrix-only MATLAB function for neural network code
      % generation with MATLAB Coder tools.
      genFunction(net,functionName,'MatrixOnly','yes');
    end
    if (false)
      % Generate a Simulink diagram for simulation or deployment with.
      % Simulink Coder tools.
      gensim(net);
    end


end