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

options.GP.aInput = {'HYST', 'per', 'strength', 'hard'};
options.GP.bInput = {'HYST', 'per', 'strength', 'hard'};
options.GP.sigmaInput = {'HYST', 'per', 'strength', 'hard'};

options.General.minPointsPSDM = 70;
MUds4 = 250; % used as limit for the analyses

options.PSDM.fit = 'bi-linear'; %Either 'power law' or 'bi-linear'
options.PSDM.IMused = 'SA';

kFold = 10;

allData = LRB; % Database with equivalent SDoF cloud-NLTHAs results 
%% Training
tic

fullFit = surrogatedPSDMdisp_ISO(allData, options);

fullFit = fullFit.getTrainingPSDM(MUds4);

fullFit = fullFit.fitGPregression;

save('fullFitDisp.mat', 'fullFit')
toc
%% Validation 

% test within the training set
fullFit = fullFit.getErrorGP;

% k-fold cross validation
[summaryKfold, kFoldGPs] = fullFit.kFoldGPs(kFold);

save('fullFitKfoldDISP.mat', 'summaryKfold', 'kFoldGPs')
