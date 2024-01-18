%% Options

System=FPS;

optionsGP.OptimizeHyperparameters = 'auto';

optionsGP.Optimizer = 'quasinewton';
optionsGP.KernelFunction = 'ardsquaredexponential';
optionsGP.Standardize = true;
optionsGP.HyperparameterOptimizationOptions = ...
    struct('UseParallel',true, 'Verbose', 0);

options.GP.optionsGP = optionsGP;
options.GP.slopeInput = {'HYST', 'per', 'strength', 'hard'};
options.GP.sigmaInput = {'HYST', 'per', 'strength', 'hard'};

options.General.minPointsPSDM = 70;
AccRatioLim = 5; % Used as limit for the analyses

options.PSDM.fit = 'bi-linear'; %Either 'power law' or 'bi-linear'
options.PSDM.IMused = 'SA';

kFold = 10;

allData = LRB; % Database with equivalent SDoF cloud-NLTHAs results 
%% Training
tic

% load('LRB.mat')

% Express acceleration in terms of g
    for k=length(allData.per):-1:1 
        allData.edpTHacc{k}=allData.edpTHacc{k,1}./9.81;
    end

fullFitAcc = surrogatedPSDMacc_ISO(allData, options);

fullFitAcc = fullFitAcc.getTrainingPSDM(AccRatioLim);

fullFitAcc = fullFitAcc.fitGPregression;

save('fullFitAcc.mat', 'fullFitAcc')
toc
%% Validation 

% test within the training set
fullFitAcc = fullFitAcc.getErrorGP;

% k-fold cross validation
[summaryKfold, kFoldGPs] = fullFitAcc.kFoldGPs(kFold);

save('fullFitKfoldAcc.mat', 'summaryKfold', 'kFoldGPs')
