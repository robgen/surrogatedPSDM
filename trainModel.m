%% Options

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
MUds4 = 7; % ductility at DS4 (used as limit for the analyses)

kFold = 10;

%% Training

load('allDataTrain.mat')

fullFit = surrogatedPSDM(allData, options);

fullFit = fullFit.getTrainingPSDM(MUds4);

fullFit = fullFit.fitGPregression;

save('fullFit.mat', 'fullFit')

%% Validation 

% test within the training set
fullFit = fullFit.getErrorGP;

% k-fold cross validation
[summaryKfold, kFoldGPs] = fullFit.kFoldGPs(kFold);

save('fullFitKfold.mat', 'summaryKfold', 'kFoldGPs')
