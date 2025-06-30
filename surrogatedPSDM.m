classdef surrogatedPSDM
    %surrogatedPSDM
    
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

        nameRegs
    end
    
    
    methods
        
        function self = surrogatedPSDM(allData, options)
            % surrogatePSDM is the construct of this class
            
            if nargin == 1, options = struct; end
            self = setAllParameters(self, options);
            
            self.allData = allData;

            if strcmp(self.parameters.PSDM.fx, 'bilinear')
                self.nameRegs = {'slope', 'sigma'};
            elseif strcmp(self.parameters.PSDM.fx, 'powerlaw')
                self.nameRegs = {'a', 'b', 'sigma'};
            else
                error('fx can only be bilinear or powerlaw')
            end
            
        end
        
        
        function self = getTrainingPSDM(self, MUds4)
            
            if nargin > 1
                if numel(MUds4)==1
                    MUds4 = MUds4 * ones(size(self.allData,1),1);
                end
            else
                for m = size(self.allData,1) : -1 : 1
                    MUds4(m,1) = self.allData.PUSH{1}(4,1) / ...
                        self.allData.PUSH{1}(3,1);
                end
            end
            
            if strcmp(self.parameters.PSDM.fx, 'bilinear')
                for m = size(self.allData,1) : -1 : 1
                    [self.PSDMtrain(m,1), self.PSDMtrain(m,2), ...
                        self.PSDMtrain(m,3)] = self.fitPSDMbilinear(...
                        self.allData.SA{m} / self.allData.strength(m), ...
                        self.allData.edpTHmax{m} / self.allData.Dyield(m), ...
                        MUds4(m));
                end
                
                [~,~,~,self.PSDMformula, self.fragFormula] = self.fitPSDMbilinear(...
                        self.allData.SA{m} / self.allData.strength(m), ...
                        self.allData.edpTHmax{m} / self.allData.Dyield(m), ...
                        MUds4(m));

                self.isAcceptable = self.PSDMtrain(:,3) >= ...
                    self.parameters.PSDM.minPoints;
            else
                for m = size(self.allData,1) : -1 : 1
                    [self.PSDMtrain(m,1), self.PSDMtrain(m,2), ...
                        self.PSDMtrain(m,3), self.PSDMtrain(m,4)] = self.fitPSDMpowerlaw(...
                        self.allData.SA{m} / self.allData.strength(m), ...
                        self.allData.edpTHmax{m} / self.allData.Dyield(m), ...
                        MUds4(m));
                end
                
                [~,~,~,~, self.PSDMformula, self.fragFormula] = self.fitPSDMpowerlaw(...
                        self.allData.SA{m} / self.allData.strength(m), ...
                        self.allData.edpTHmax{m} / self.allData.Dyield(m), ...
                        MUds4(m));

                self.isAcceptable = self.PSDMtrain(:,4) >= ...
                    self.parameters.PSDM.minPoints;
            end
            
        end
        
        
        function self = fitGPregression(self)
            
            optionsNameValue = self.struct2pairs(self.parameters.GP.optionsGP);
            
            for r = 1 : numel(self.nameRegs)
                inputInd = getInputIndices(self, self.allData);
                
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
            
            
            figure; hold on
            for r = 1 : numel(self.nameRegs)
                scatter(summaryKfold.(self.nameRegs{r}).inFoldPred, ...
                    self.GPregs.(self.nameRegs{r}).Y, ...
                    50, 'filled', 'DisplayName', ...
                    sprintf('%s - inFold NRMSE=%1.3f', self.nameRegs{r}, ...
                    summaryKfold.(self.nameRegs{r}).inFoldNRMSE))
            end
            legend('Location', 'SouthEast')
            xlabel('Modelled')
            ylabel('Surrogated')
            set(gca, 'FontSize', 18)
            
            
            figure; hold on
            for r = 1 : numel(self.nameRegs)
                histogram(summaryKfold.(self.nameRegs{r}).inFoldPred ./ ...
                    self.GPregs.(self.nameRegs{r}).Y, ...
                    'normalization', 'pdf')
            end
            legend(self.nameRegs, 'Location', 'NorthEast')
            xlabel('Surrogated/Modelled')
            ylabel('Frequency (PDF)')
            set(gca, 'FontSize', 18)
        end
        
        
        function self = getErrorGP(self)
            
            nReg = numel(self.nameRegs);
            for r = 1 : nReg
                inputInd = getInputIndices(self, self.allData);
                
                self.PSDMtrainGP(:,r) = predict(...
                    self.GPregs.(self.nameRegs{r}), ...
                    self.allData(:,inputInd) );
            end
            
            NRMSE = ( mean((self.PSDMtrainGP(self.isAcceptable,:) - ...
                self.PSDMtrain((self.isAcceptable),1:nReg)).^2) ) .^0.5 ./ ...
                mean(self.PSDMtrain((self.isAcceptable),1:nReg));
            
            maxPlot = max(max(self.PSDMtrainGP));
            
            figure; hold on
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
            set(gca', 'FontSize', 18)
            
            
            figure; hold on
            for r = 1 : nReg
                histogram(self.PSDMtrainGP(self.isAcceptable,r) ./ ...
                    self.PSDMtrain(self.isAcceptable,r), ...
                    'normalization', 'pdf')
            end
            legend(self.nameRegs, 'Location', 'NorthEast')
            xlabel('Surrogated/Modelled')
            ylabel('Frequency (PDF)')
            set(gca', 'FontSize', 18)
            
        end
        
        
        function [fragMedian, fragStDev, powerLaw, ...
                isExtrapolated, maxAllowedIM] = predictFragGP(self, ...
                hyst, T, fy, hard, ductThresholds)
            % related to mu-1=a(R-1)
            
            % hyst must be a column cell array
            inputTable = table;
            inputTable.HYST = hyst;
            inputTable.per(:) = T;
            inputTable.strength(:) = fy;
            inputTable.hard(:) = hard;
            
            % predict power law parameters
            isExtrapolated = zeros(numel(T),1);
            for r = numel(self.nameRegs) : -1 : 1
                inputInd = getInputIndices(self, inputTable);
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

            if nargin == 6
                for m = numel(T) : -1 : 1
                    if strcmp(self.parameters.PSDM.fx, 'bilinear')
                        % it will fail if there are no elastic cases
                        elasticCases = ductThresholds(m,:) <= 1;

                        fragMedian(m,elasticCases) = ...
                            ductThresholds(m,elasticCases) * fy(m);
                        fragMedian(m,~elasticCases) = self.fragFormula.eta(...
                            ductThresholds(m,~elasticCases), powerLaw(m,1)) * fy(m);
                        
                        fragStDev(m,elasticCases) = 0.01;
                        fragStDev(m,~elasticCases) = self.fragFormula.beta(...
                            powerLaw(m,2));
                    else 
                        fragMedian(m,:) = self.fragFormula.eta(...
                            ductThresholds(m,~elasticCases), ...
                            powerLaw(m,1), powerLaw(m,2)) * fy(m);
                        
                        fragStDev(m,:) = self.fragFormula.beta(...
                            powerLaw(m,2), powerLaw(m,3));
                    end

                end
            else
                fragMedian = NaN;
                fragStDev = NaN;
            end

            if nargout == 5
                maxAllowedIM = self.getMaxAllowedIM(inputTable);
            else
                maxAllowedIM = NaN;
            end
        end
        
        
        function maxAllowedIM = getMaxAllowedIM(self, inputTable)
            
            maxAllowedIM = zeros(size(inputTable,1),1);

            hyst = unique(inputTable.HYST);
            for h = 1 : numel(hyst)
                oneHystInputTable = inputTable(strcmp(inputTable.HYST, hyst(h)),:);
                oneHystAllData = self.allData(strcmp(self.allData.HYST, hyst(h)),:);

                gridPoints = [oneHystAllData.per, ...
                              oneHystAllData.strength, ...
                              oneHystAllData.hard];
                kdtree = KDTreeSearcher(gridPoints);

                for p = size(oneHystInputTable,1) : -1 : 1
                    targetPoint = [oneHystInputTable.per(p), ...
                        oneHystInputTable.strength(p), ...
                        oneHystInputTable.hard(p)];
                    [index, ~] = knnsearch(kdtree, targetPoint);
                    oneHystMaxAllowedIM(p,1) = max(self.allData.SA{index});
                end
                maxAllowedIM(strcmp(inputTable.HYST, hyst(h))) = oneHystMaxAllowedIM;
            end
        end

    end
    
    
    methods(Access = 'private')
        
        function self = setAllParameters(self, options)
            % setAllParameters deals with the optional parameters
            
            defOptionsGP = self.setDefaultGP;
            
            % build basic parameters
            macroFieldsPar = {'GP', 'Parametric', 'PSDM'};
            
            microFieldsPar{1} = {'optionsGP', 'variablesInput'};
            microFieldsParVals{1} = { defOptionsGP, {'HYST', 'per', 'strength', 'hard'}};
            
            microFieldsPar{2} = {'Dummy'};
            microFieldsParVals{2} = {{''}};
            
            microFieldsPar{3} = {'minPoints', 'fx' };
            microFieldsParVals{3} = { 10, 'bilinear' };
            
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
        
        
        function inputInd = getInputIndices(self, targetTable)
            inputNames = self.parameters.GP.variablesInput;
            for p = numel(inputNames) : -1 : 1
                inputInd(p) = find(...
                    strcmp(targetTable.Properties.VariableNames, ...
                    inputNames{p}));
            end
        end
        

    end
    
    
    methods(Static)
        
        function [a, sigma, NpointsReg, PSDMformula, fragFormula] = ...
                fitPSDMbilinear(IM, EDP, EDPds4)
            %fitPSDM fits the law MU = a.*(R-1) + 1; R = SA/SAy
            % if you want different fitting methods, copy from "comparePSDMs"
            
            % the 1.05 limit (instead of 1) avoids logs too close to zero
            toInclude = EDP > 1.05 & IM > 1 & EDP <= EDPds4;
            NpointsReg = sum(toInclude);
            
            b = 1;
            lna = mean( log(EDP(toInclude)-1) - b*log(IM(toInclude)-1) );

            a = exp(lna);

            resid = log(EDP(toInclude)-1) - (lna + b.*log(IM(toInclude)-1));
            sigma = std( resid );

            PSDMformula.mu = @(R,a,sigma) ( a.*exp(sigma).*(R-1) + 1 );
            PSDMformula.R = @(mu,a,sigma) ( (mu-1)./(a.*exp(sigma)) + 1 );
            
            fragFormula.eta = @(muCap,a) ((muCap-1)/a + 1);
            fragFormula.beta = @(sigma) (sigma);
            
%             %%% Control plot
%             dummyMU = 1:0.01:5;
%             figure; hold on
%             scatter(IM(toInclude), EDP(toInclude))
%             plot([0 PSDMformula.R(dummyMU, a, 0)], [0 dummyMU], 'k', ...
%                 'LineWidth', 2)
%             plot([0 PSDMformula.R(dummyMU, a, sigma)], ...
%                 [0 dummyMU], '--k', 'LineWidth', 1)
%             plot([0 PSDMformula.R(dummyMU, a, -sigma)], ...
%                 [0 dummyMU], '--k', 'LineWidth', 1)
%             xlabel('R [-]')
%             ylabel('\mu [-]')
%             set(gca, 'FontSize', 18)
            
        end
        
        
        function [a, b, sigma, NpointsReg, PSDMformula, fragFormula] = ...
                fitPSDMpowerlaw(IM, EDP, EDPds4)
            %fitPSDM fits the law MU = a.*R^b; R = SA/SAy
            
            % the 1.05 limit (instead of 1) avoids logs too close to zero
            toInclude = EDP > 1.00 & IM > 1 & EDP <= EDPds4;
            NpointsReg = sum(toInclude);
            
            [xData, yData] = prepareCurveData( ...
                log(IM(toInclude)), log(EDP(toInclude)) );
            
            lm = fitlm(xData, yData);
            a = exp(lm.Coefficients.Estimate(1));
            b = lm.Coefficients.Estimate(2);
            
            resid = log(EDP(toInclude)) - (log(a) + b.*log(IM(toInclude)));
            sigma = std( resid );
            
            % TODO: check formulas
            PSDMformula.mu = @(R,a,b,sigma) a.*R.^b + 0*sigma;
            PSDMformula.R = @(mu,a,b,sigma) exp((log(mu)-log(a))/b ) + 0*sigma;

            % TODO: check formulas
            fragFormula.eta = @(muCap,a,b) exp((log(mu)-log(a))/b );
            fragFormula.beta = @(sigma,b) (sigma/b);
            
%             %%% Control plot
%             dummyMU = 0:0.01:6;
%             figure; hold on
%             scatter(IM(toInclude), EDP(toInclude))
%             plot(PSDMformula.R(dummyMU, a, b, 0), dummyMU, 'k', ...
%                 'LineWidth', 2)
%             plot([0 PSDMformula.R(dummyMU, a, b, sigma)], ...
%                 [0 dummyMU], '--k', 'LineWidth', 1)
%             plot([0 PSDMformula.R(dummyMU, a, b, -sigma)], ...
%                 [0 dummyMU], '--k', 'LineWidth', 1)
%             xlabel('R [-]')
%             ylabel('\mu [-]')
%             set(gca, 'FontSize', 18)
            
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

