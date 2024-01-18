classdef surrogatedPSDMacc
    %surrogatedPSDM for acc on top of isolation system
    
    properties
        parameters
        
        allData
        PSDMtrain
        isAcceptable
        GPregs
        PSDMtrainGP % predicted w/ the GP using training input
        PSDMformula
        fragFormula
        
        
    end
    
    
    properties(Access = 'private')

        nameRegs = {'slope', 'sigma'};
    end
    
    
    methods
        
        function self = surrogatedPSDMacc(allData, options)
            % surrogatePSDM is the construct of this class
            
            if nargin == 1, options = struct; end
            self = setAllParameters(self, options);
            
            self.allData = allData;
            
        end
        
        
        % Modified EDP as Acc/Accyieldiso
        function self = getTrainingPSDM(self, AccRatioLim)
                            
            if nargin > 1
                if numel(AccRatioLim)==1
                    AccRatioLim = AccRatioLim * ones(size(self.allData,1),1);
                end
            else
                for m = size(self.allData,1) : -1 : 1
                    AccRatioLim(m,1) = self.allData.PUSH{1}(4,1) / ...
                        self.allData.PUSH{1}(3,1);
                end
            end
            
            for m = size(self.allData,1) : -1 : 1
                [self.PSDMtrain(m,1), self.PSDMtrain(m,2), ...
                    self.PSDMtrain(m,3)] = self.fitPSDM(...
                    self.allData.SA{m} / self.allData.strength(m), ...
                    self.allData.edpTHacc{m,1}'/self.allData.strength(m), ...
                    AccRatioLim(m));
            end
            
            [~,~,~,self.PSDMformula, self.fragFormula] = self.fitPSDM(...
                    self.allData.SA{m} / self.allData.strength(m), ...
                    self.allData.edpTHacc{m,1}'/self.allData.strength(m), ...
                    AccRatioLim(m));
            
            self.isAcceptable = self.PSDMtrain(:,3) >= ...
                self.parameters.General.minPointsPSDM;
        end
        
        
        function self = fitGPregression(self)
            
            optionsNameValue = self.struct2pairs(self.parameters.GP.optionsGP);
            
            for r = 1 : numel(self.nameRegs)
                inputInd = getInputIndices(self, self.allData, self.nameRegs{r});
                
                self.GPregs.(self.nameRegs{r}) = fitrgp(...
                    self.allData(self.isAcceptable,inputInd), ...
                    self.PSDMtrain(self.isAcceptable,r), optionsNameValue{:});
            end
        end
        
        
        function [summaryKfold, kFoldGPs] = kFoldGPs(self, k)
            
            if nargin < 2;  k = 10; end
            
            for r = 1 : numel(self.nameRegs)
                kFoldGPs.(self.nameRegs{r}) = crossval(...
                    self.GPregs.(self.nameRegs{r}), 'kfold', k);
                
                summaryKfold.(self.nameRegs{r}).inFoldLoss = kfoldLoss(...
                    kFoldGPs.(self.nameRegs{r}));
                
                summaryKfold.(self.nameRegs{r}).inFoldPred(:,1) = ...
                    kfoldPredict(kFoldGPs.(self.nameRegs{r}));
                
                summaryKfold.(self.nameRegs{r}).inFoldNRMSE = ( mean((...
                    summaryKfold.(self.nameRegs{r}).inFoldPred(:,1) - ...
                    self.GPregs.(self.nameRegs{r}).Y ).^2) ) .^0.5 ./ ...
                    mean(self.GPregs.(self.nameRegs{r}).Y);
            end
            
            
            figure;F=figure('Position', get(0, 'Screensize')); hold on
%             title('inFold Norm. Root Mean Square Deviation  - Acc PSDM')
            for r = 1 : numel(self.nameRegs)
                scatter(summaryKfold.(self.nameRegs{r}).inFoldPred, ...
                    self.GPregs.(self.nameRegs{r}).Y, ...
                    50, 'filled','MarkerEdgeColor',[0.3 .3 .3],...
                    'LineWidth',1.0,'DisplayName', ...
                    sprintf('%s - inFold NRMSE=%1.3f', self.nameRegs{r}, ...
                    summaryKfold.(self.nameRegs{r}).inFoldNRMSE))
            end
            legend('Location', 'SouthEast')
            plot([0 0.6],[0 0.6],'color',[0.5 0.5 0.5])
            xlabel('Modelled')
            ylabel('Surrogated')
            set(gca, 'FontSize', 16)
            saveas(F,'Figures\NRMSEinFold_Acc.png','png')
            
            
            figure;F=figure(); hold on
%             title('Surrogated vs Modelled inFold Distribution - Acc PSDM')
            for r = 1 : numel(self.nameRegs)
                histogram(summaryKfold.(self.nameRegs{r}).inFoldPred ./ ...
                    self.GPregs.(self.nameRegs{r}).Y, ...
                    'normalization', 'pdf')
            end
            legend(self.nameRegs, 'Location', 'NorthEast')
            xlabel('Surrogated/Modelled')
            ylabel('Frequency (PDF)')
            set(gca, 'FontSize', 16)
            saveas(F,'Figures\SvsMinFoldDist_Acc.png','png')
            
        end
        
        
        function self = getErrorGP(self)
            
            nReg = numel(self.nameRegs);
            for r = 1 : nReg
                inputInd = getInputIndices(self, ...
                    self.allData, self.nameRegs{r});
                
                self.PSDMtrainGP(:,r) = predict(...
                    self.GPregs.(self.nameRegs{r}), ...
                    self.allData(:,inputInd) );
            end
            
            NRMSE = ( mean((self.PSDMtrainGP(self.isAcceptable,:) - ...
                self.PSDMtrain((self.isAcceptable),1:nReg)).^2) ) .^0.5 ./ ...
                mean(self.PSDMtrain((self.isAcceptable),1:nReg));
            
            maxPlot = max(max(self.PSDMtrainGP));
            
            figure;F=figure(); hold on
%             title('Norm. Root Mean Square Deviation  - Acc PSDM')
            plot(maxPlot*[0 1.2], maxPlot*[0 1.2], 'k')
            for r = 1 : nReg
                scat(r) = scatter(self.PSDMtrain(self.isAcceptable,r), ...
                    self.PSDMtrainGP(self.isAcceptable,r), ...
                    50, 'filled', 'DisplayName', ...
                    sprintf('%s - NRMSE=%1.3f', self.nameRegs{r}, NRMSE(r)));
            end
            cols = get(gca,'ColorOrder');
            for r = 1 : nReg
                scatter(self.PSDMtrain(~self.isAcceptable,r), ...
                    self.PSDMtrainGP(~self.isAcceptable,r), ...
                    50, 'MarkerEdgeColor', cols(r,:))
            end
            legend(scat, 'Location', 'SouthEast')
            xlabel('Modelled')
            ylabel('Surrogated')
            set(gca', 'FontSize', 16)
            saveas(F,'Figures\NRMSE_Acc.png','png')
            
            figure;F=figure(); hold on
%             title('Surrogated vs Modelled Distribution - Acc PSDM')
            for r = 1 : nReg
                histogram(self.PSDMtrainGP(self.isAcceptable,r) ./ ...
                    self.PSDMtrain(self.isAcceptable,r), ...
                    'normalization', 'pdf')
            end
            legend(self.nameRegs, 'Location', 'NorthEast')
            xlabel('Surrogated/Modelled')
            ylabel('Frequency (PDF)')
            set(gca', 'FontSize', 16)
            saveas(F,'Figures\SvsMDist_Acc.png','png')
            
        end
        
    
        function [powerLaw, isExtrapolated] = getPowerLawAcc(self,...
                hyst, T, fyiso, alphaiso)
                    % related to Acc/fyiso-1=a(R-1)

            % hyst must be a column cell array
            inputTable = table;
            inputTable.HYST = hyst;
            inputTable.per(:) = T;
            inputTable.strength(:) = fyiso;
            inputTable.hard(:) = alphaiso;
             
            % predict power law parameters
            isExtrapolated = zeros(numel(T),1);
            for r = numel(self.nameRegs) : -1 : 1
                inputInd = getInputIndices(self, ...
                    inputTable, self.nameRegs{r});
                powerLaw(:,r) = predict(...
                    self.GPregs.(self.nameRegs{r}), inputTable(:,inputInd) );

                % find extrapolations
                namePredictors = self.GPregs.(self.nameRegs{r}).PredictorNames;
                namePredictors(strcmp(namePredictors, 'HYST')) = [];
                for p = 1 : numel(namePredictors)
                    minP = min(self.GPregs.(self.nameRegs{r}).X.(namePredictors{p}));
                    maxP = max(self.GPregs.(self.nameRegs{r}).X.(namePredictors{p}));
                    isExtrapolated = isExtrapolated | ...
                        inputTable.(namePredictors{p}) < minP | ...
                        inputTable.(namePredictors{p}) > maxP;
                end
            end
        end
        
        
        function [fragMedian, fragStDev] = predictFragGPSS(self, ...
                fyiso,fySS, accRatioThresholdsSS,powerLaw)
            % Returns the fragility curves median and standard dev. for all
            % the SS and isolation systems where their PSDM are predicted 
            % by powerLaw the limit states are defined  in accRatioThresholdsSS 

            for SS=length(fySS):-1:1 
                if nargin == 5
                    for m = numel(fyiso) : -1 : 1
                    
                    elasticCases = accRatioThresholdsSS(m,:,SS) <= 1;
                    
                    fragMedian(m,elasticCases,SS) = ...
                        accRatioThresholdsSS(m,elasticCases,SS) * fyiso(m);
                    fragMedian(m,~elasticCases,SS) = self.fragFormula.eta(...
                        accRatioThresholdsSS(m,~elasticCases,SS), powerLaw(m,1)) * fyiso(m);
                    
                    fragStDev(m,elasticCases,SS) = 0.01;
                    fragStDev(m,~elasticCases,SS) = self.fragFormula.beta(...
                        powerLaw(m,2));
                    
                    end
                else
                    fragMedian = NaN;
                    fragStDev = NaN;
                end
            end
        end
        
        
        function [fragMedian, fragStDev] = predictFragGPNSCD(self, ...
                fyiso,fySS, accRatioThresholdsNSCD,powerLaw)
            
            % Returns the fragility curves median and standard dev. for all
            % the SS, isolation systems and storeys where their PSDM are 
            % predicted by powerLaw. 
            % Computes the fragility curves for the NSCD damage 
            % defined by the limit states in accRatioThresholdsNSCD

            for SS=length(fySS):-1:1 
                if nargin == 5
                   for ST= size(accRatioThresholdsNSCD,3):-1:1
                        for m = numel(fyiso) : -1 : 1

                        elasticCases = accRatioThresholdsNSCD(m,:,ST,SS) <= 1;
                        fragMedian(m,elasticCases,ST,SS) = ...
                            accRatioThresholdsNSCD(m,elasticCases,ST,SS) * fyiso(m);
                        fragMedian(m,~elasticCases,ST,SS) = self.fragFormula.eta(...
                            accRatioThresholdsNSCD(m,~elasticCases,ST,SS), powerLaw(m,1)) * fyiso(m);

                        fragStDev(m,elasticCases,ST,SS) = 0.01;
                        fragStDev(m,~elasticCases,ST,SS) = self.fragFormula.beta(...
                            powerLaw(m,2));
                        end
                    end
                else
                    fragMedian = NaN;
                    fragStDev = NaN;
                end
            end
        end
        
        
        
        function [fragMedian, fragStDev] = predictFragGPNSCA(self, ...
                fyiso,fySS, accRatioThresholdsNSCA,powerLaw)
            
            % Returns the fragility curves median and standard dev. for all
            % the SS, isolation systems and storeys where their PSDM are 
            % predicted by powerLaw. 
            % Computes the fragility curves for the NSCA damage 
            % defined by the limit states in accRatioThresholdsNSCA

            for SS=length(fySS):-1:1 
                if nargin == 5
                   for ST= size(accRatioThresholdsNSCA,3):-1:1
                        for m = numel(fyiso) : -1 : 1

                        elasticCases = accRatioThresholdsNSCA(m,:,ST,SS) <= 1;
                        fragMedian(m,elasticCases,ST,SS) = ...
                            accRatioThresholdsNSCA(m,elasticCases,ST,SS) * fyiso(m);
                        fragMedian(m,~elasticCases,ST,SS) = self.fragFormula.eta(...
                            accRatioThresholdsNSCA(m,~elasticCases,ST,SS), powerLaw(m,1)) * fyiso(m);

                        fragStDev(m,elasticCases,ST,SS) = 0.01;
                        fragStDev(m,~elasticCases,ST,SS) = self.fragFormula.beta(...
                            powerLaw(m,2));
                        end
                    end
                else
                    fragMedian = NaN;
                    fragStDev = NaN;
                end
            end
        end
        
    end

    
    
    
    methods(Access = 'private')
        
        function self = setAllParameters(self, options)
            % setAllParameters deals with the optional parameters
            
            defOptionsGP = self.setDefaultGP;
            
            % build basic parameters
            macroFieldsPar = {'GP', 'Parametric', 'General'};
            
            microFieldsPar{1} = {'optionsGP', 'slopeInput', 'sigmaInput'};
            microFieldsParVals{1} = { defOptionsGP, {'HYST', 'per'}, {'HYST', 'per'} };
            
            microFieldsPar{2} = {'Dummy'};
            microFieldsParVals{2} = {{''}};
            
            microFieldsPar{3} = {'minPointsPSDM' };
            microFieldsParVals{3} = { 10 };
            
            for F = 1 : numel(macroFieldsPar)
                for f = 1 : numel(microFieldsPar{F})
                    self.parameters.(macroFieldsPar{F}).(microFieldsPar{F}{f}) = microFieldsParVals{F}{f};
                end
            end
            
            % overwrite fields if some parameter is specified
            macroFieldsOptional = fieldnames(options);
            for OF = 1 : numel(macroFieldsOptional)
                microFieldsOptional = fieldnames(options.(macroFieldsOptional{OF}));
                for of = 1 : numel(microFieldsOptional)
                    if isfield(self.parameters.(macroFieldsOptional{OF}), microFieldsOptional{of}) == 1
                        self.parameters.(macroFieldsOptional{OF}).(microFieldsOptional{of}) = options.(macroFieldsOptional{OF}).(microFieldsOptional{of});
                    else
                        error('Field %s.%s is not allowed', macroFieldsOptional{OF}, microFieldsOptional{of})
                    end
                end
            end
        end
        
        
        function inputInd = getInputIndices(self, targetTable, nameReg)
            inputNames = self.parameters.GP.([nameReg 'Input']);
            for p = numel(inputNames) : -1 : 1
                inputInd(p) = find(...
                    strcmp(targetTable.Properties.VariableNames, ...
                    inputNames{p}));
            end
        end
        
        
    end
    
    
    methods(Static)
        

        function [a, sigma, NpointsReg, PSDMformula, fragFormula] = ...
                fitPSDM(IM, EDP, AccRatLim) %ds4=max(fySS/fyiso)
            
            % fitPSDM fits the law Acc/fyiso = a.*(R-1) + 1; R = SA/SAy
            % if you want different fitting methods, copy from "comparePSDMs"
            
            % the 1.05 limit (instead of 1) avoids logs too close to zero
            toInclude = EDP > 1.05 & IM > 1 & EDP <= AccRatLim; 
            NpointsReg = sum(toInclude);
            
            b = 1;
            lna = mean( log(EDP(toInclude)-1) - b*log(IM(toInclude)-1) );

            a = exp(lna);

            resid = log(EDP(toInclude)-1) - (lna + b.*log(IM(toInclude)-1));
            sigma = std( resid );

            PSDMformula.AccRatio = @(R,a,sigma) ( a.*exp(sigma).*(R-1) + 1 );
            PSDMformula.R = @(AccRatio,a,sigma) ( (AccRatio-1)./(a.*exp(sigma)) + 1 );
            
            fragFormula.eta = @(AccRatio,a) ((AccRatio-1)/a + 1);
            fragFormula.beta = @(sigma) (sigma);
            
%             Control plot
%             dummyAccRatio = 1:0.01:5;
%             figure; hold on
%             title('Acc-based PSDM')
%             grid(gca,'minor');
%             scatter(IM(toInclude), EDP(toInclude),36,'MarkerEdgeColor',[128 128 128]/265,'MarkerFaceColor',[128 128 128]/285)
%             plot([0 PSDMformula.R(dummyAccRatio, a, 0)], [0 dummyAccRatio], 'k', ...
%                 'LineWidth', 3)
%             plot([0 PSDMformula.R(dummyAccRatio, a, sigma)], ...
%                 [0 dummyAccRatio], '--k', 'LineWidth', 1.5,'Color',[0.5,0.5,0.5])
%             plot([0 PSDMformula.R(dummyAccRatio, a, -sigma)], ...
%                 [0 dummyAccRatio], '--k', 'LineWidth', 1.5,'Color',[0.5,0.5,0.5])
%             set(gca,'YTick',[1 6],'XTick',[1 26])
%             ylim([0,5])
%             xlim([0,30])
%             xlabel('R [-]')
%             ylabel('\alpha [-]')
%             set(gca, 'FontSize', 15)
            
        end
        
        
        function C = struct2pairs(S)
            %Turns a scalar struct S into a cell of string-value pairs C
            
            if iscell(S)
                C=S; return
            elseif length(S)>1
                error 'Input must be a scalar struct or cell';
            end
            C=[fieldnames(S).'; struct2cell(S).'];
            C=C(:).';
        end
        
        
        function defOptionsGP = setDefaultGP
                       
            defOptionsGP.Optimizer = 'quasinewton';
            defOptionsGP.Basis = 'constant';
            defOptionsGP.KernelFunction = 'squaredexponential'; % ARD kernels cannot be optimised.
            defOptionsGP.Sigma = 0.00001;
            defOptionsGP.constantSigma = false;
            defOptionsGP.Standardize = true;
            defOptionsGP.OptimizeHyperparameters = {'KernelScale'};
            defOptionsGP.HyperparameterOptimizationOptions = struct('UseParallel',false, 'Verbose', 0);
        end
        
        
    end
end

