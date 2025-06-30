%% Options

optionsGP.OptimizeHyperparameters = 'auto';

optionsGP.Optimizer = 'quasinewton';
optionsGP.KernelFunction = 'ardsquaredexponential';
optionsGP.Standardize = true;
optionsGP.HyperparameterOptimizationOptions = ...
    struct('UseParallel',true, 'Verbose', 0);

options.GP.optionsGP = optionsGP;
options.GP.variablesInput = {'HYST', 'per', 'strength', 'hard'};

options.PSDM.minPoints = 70;
options.PSDM.fx = 'powerlaw'; % 'bilinear', 'powerlaw'
MUds4 = 7; % ductility at DS4 (used as limit for the analyses)

kFold = 10;

%% Training

load('allDataTrain.mat')

fullFit = surrogatedPSDM(allData, options);

fullFit = fullFit.getTrainingPSDM(MUds4);

fullFit = fullFit.fitGPregression;

%% Validation 

% test within the training set
fullFit = fullFit.getErrorGP;

% k-fold cross validation
[summaryKfold, kFoldGPs] = fullFit.kFoldGPs(kFold);

%% Save

if strcmp(options.PSDM.fx, 'bilinear')
    save('fullFit.mat', 'fullFit')
    save('fullFitKfold.mat', 'summaryKfold', 'kFoldGPs')
else
    save('fullFitPowerLaw.mat', 'fullFit', '-v7.3')
    save('fullFitKfoldPowerLaw.mat', 'summaryKfold', 'kFoldGPs', '-v7.3')
end
