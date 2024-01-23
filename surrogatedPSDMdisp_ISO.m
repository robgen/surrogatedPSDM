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
        
        nameRegs;

    end
    
    
    methods
        
        function self = surrogatedPSDM(allData, options)
            % surrogatePSDM is the construct of this class
            
            if nargin == 1, options = struct; end
            self = setAllParameters(self, options);
            
            self.allData = allData;

            switch self.parameters.PSDM.fit
                case 'bi-linear'
                    self.nameRegs = {'slope', 'sigma'};

                case 'power law'
                    self.nameRegs = {'a', 'b', 'sigma'};
            end
            
        end
        
        
        function self = getTrainingPSDM(self, Mulim)
            
            if nargin > 1
                if numel(Mulim)==1
                    Mulim = Mulim * ones(size(self.allData,1),1);
                end
            else
                for m = size(self.allData,1) : -1 : 1
                    Mulim(m,1) = self.allData.PUSH{1}(4,1) / ...
                        self.allData.PUSH{1}(3,1);
                end
            end

            IMused = self.parameters.PSDM.IMused;

            switch self.parameters.PSDM.fit

            case 'bi-linear'

                for m = size(self.allData,1) : -1 : 1
                    [self.PSDMtrain(m,1), self.PSDMtrain(m,2), ...
                        self.PSDMtrain(m,3),self.PSDMtrain(:,4),~] = ...
                        self.fitPSDMbilinear(...
                        self.allData.(IMused){m}', self.allData.R{m}', ...
                        self.allData.edpTHmax{m} / self.allData.Dyield(m), ...
                        Mulim(m));
                end
                
                [~,~,~,~,self.PSDMformula, self.fragFormula] = self.fitPSDMbilinear(...
                        self.allData.(IMused){m}', self.allData.R{m}', ...
                        self.allData.edpTHmax{m} / self.allData.Dyield(m), ...
                        Mulim(m));
                
                self.isAcceptable = self.PSDMtrain(:,4) >= ...
                    self.parameters.General.minPointsPSDM;

            case 'power law'

                for m = size(self.allData,1) : -1 : 1
                    [self.PSDMtrain(m,1), self.PSDMtrain(m,2), ...
                        self.PSDMtrain(m,3),self.PSDMtrain(m,4),~,~] = ...
                        self.fitPSDMPowerLaw(...
                        self.allData.(IMused){m}',  ...
                        self.allData.edpTHmax{m} / self.allData.Dyield(m), ...
                        0, Mulim(m));
                end
                
                [~,~,~,~,self.PSDMformula, self.fragFormula] = ...
                    self.fitPSDMPowerLaw(...
                        self.allData.(IMused){m}', ...
                        self.allData.edpTHmax{m} / self.allData.Dyield(m), ...
                        0, Mulim(m));
                
                self.isAcceptable = self.PSDMtrain(:,4) >= ...
                    self.parameters.General.minPointsPSDM;
            end
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
            
            
            figure;F=figure(); hold on
%             title('inFold Norm. Root Mean Square Deviation  - Ductility PSDM')
            for r = 1 : numel(self.nameRegs)
                scatter(summaryKfold.(self.nameRegs{r}).inFoldPred, ...
                    self.GPregs.(self.nameRegs{r}).Y, ...
                    50, 'filled','MarkerEdgeColor',[0.3 .3 .3],...
                    'LineWidth',1.0,'DisplayName', ...
                    sprintf('%s - inFold NRMSE=%1.3f', self.nameRegs{r}, ...
                    summaryKfold.(self.nameRegs{r}).inFoldNRMSE))
            end
            legend('Location', 'SouthEast')
            xlabel('Modelled')
            ylabel('Surrogated')
            set(gca, 'FontSize', 16)
            plot([0 80],[0 80],'color',[0.5 0.5 0.5])
            saveas(F,'Figures\NRMSEinFold_Disp.png','png')
            
            figure;F=figure(); hold on
%             title('Surrogated vs Modelled inFold Distribution - Ductility PSDM')
            for r = 1 : numel(self.nameRegs)
                histogram(summaryKfold.(self.nameRegs{r}).inFoldPred ./ ...
                    self.GPregs.(self.nameRegs{r}).Y, ...
                    'normalization', 'pdf')
            end
            legend(self.nameRegs, 'Location', 'NorthEast')
            xlim([0.8,1.2])
            xlabel('Surrogated/Modelled')
            ylabel('Frequency (PDF)')
            set(gca, 'FontSize', 16)
            saveas(F,'Figures\SvsMinFoldDist_Disp.png','png')
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
%             title('Norm. Root Mean Square Deviation  - Ductility PSDM')
            plot(maxPlot*[0 1.2], maxPlot*[0 1.2], 'k')
            for r = 1 : nReg
                scat(r) = scatter(self.PSDMtrain(self.isAcceptable,r), ...
                    self.PSDMtrainGP(self.isAcceptable,r), ...
                    50, 'filled','MarkerEdgeColor',[0.3 .3 .3],...
                    'LineWidth',1.0, 'DisplayName', ...
                    sprintf('%s - NRMSE=%1.3f', self.nameRegs{r}, NRMSE(r)));
            end
            cols = get(gca,'ColorOrder');
%             for r = 1 : nReg
%                 scatter(self.PSDMtrain(~self.isAcceptable,r), ...
%                     self.PSDMtrainGP(~self.isAcceptable,r), ...
%                     50, 'MarkerEdgeColor', cols(r,:))
%             end
            legend(scat, 'Location', 'SouthEast')
            xlabel('Modelled')
            ylabel('Surrogated')
            set(gca', 'FontSize', 16)
            saveas(F,'Figures\NRMSE_Disp.png','png')
            
            
            figure;F=figure(); hold on
%             title('Surrogated vs Modelled Distribution - Ductility PSDM')
            for r = 1 : nReg
                histogram(self.PSDMtrainGP(self.isAcceptable,r) ./ ...
                    self.PSDMtrain(self.isAcceptable,r), ...
                    'normalization', 'pdf')
            end
            legend(self.nameRegs, 'Location', 'NorthEast')
            xlabel('Surrogated/Modelled')
            ylabel('Frequency (PDF)')
            set(gca', 'FontSize', 16)
            saveas(F,'Figures\SvsMDist_Disp.png','png')
        end
        
        
        function [powerLaw, isExtrapolated] = getPowerLaw(self, hyst, T, fyiso, alphaiso)
            
                        % related to mu-1=a(R-1)
            
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
        
        
        function [fragMedian, fragStDev] = predictFragGPIso(self, ...
                 fyiso, ductThresholdsIso, powerLaw)
           
            if nargin == 4
                for m = numel(fyiso) : -1 : 1
                    elasticCases = ductThresholdsIso(m,:) <= 1;
                    
                    fragMedian(m,elasticCases) = ...
                        ductThresholdsIso(m,elasticCases) * fyiso(m);
                    fragMedian(m,~elasticCases) = self.fragFormula.eta(...
                        ductThresholdsIso(m,~elasticCases), powerLaw(m,1)) * fyiso(m);
                    
                    fragStDev(m,elasticCases) = 0.01;
                    fragStDev(m,~elasticCases) = self.fragFormula.beta(...
                        powerLaw(m,2));
                end
            else
                fragMedian = NaN;
                fragStDev = NaN;
            end
        end

    end
    
    
    
    methods(Access = 'private')
        
        function self = setAllParameters(self, options)
            % setAllParameters deals with the optional parameters
            
            defOptionsGP = self.setDefaultGP;
            
            % build basic parameters
            macroFieldsPar = {'GP', 'Parametric', 'General','PSDM'};
            
            microFieldsPar{1} = {'optionsGP', 'slopeInput', 'sigmaInput', 'aInput', 'bInput', 'sigmaInput'};
            microFieldsParVals{1} = { defOptionsGP, {'HYST', 'per'}, {'HYST', 'per'}, {'HYST', 'per'}, {'HYST', 'per'}, {'HYST', 'per'} };
            
            microFieldsPar{2} = {'Dummy'};
            microFieldsParVals{2} = {{''}};
            
            microFieldsPar{3} = {'minPointsPSDM' };
            microFieldsParVals{3} = { 10 };

            microFieldsPar{4} = {'fit', 'IMused' };
            microFieldsParVals{4} = { 'bi-linear', 'PGV' };
            
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
        
        function [a, b, sigma, NpointsReg, PSDMformula, fragFormula] =...
                fitPSDMPowerLaw(IM, EDP, yieldEDP,EDPds4)
        
        toInclude =  EDP <= EDPds4;
        NpointsReg = sum(toInclude);

        if yieldEDP == 0
            % single regression with all the points
            [xData, yData] = prepareCurveData( log(IM), log(EDP) );
        
            ft = fittype( 'poly1' );
            [fitresult, gof] = fit( xData, yData, ft );
            a = exp(fitresult.p2);
            b = fitresult.p1;
        else
            % line to interpolate the points up to yielding (get a)
            preYield = EDP <= yieldEDP;
            
            [xData, yData] = prepareCurveData( ...
                IM(preYield), EDP(preYield) );
            ft = fittype( {'x'} );
            fitresult = fit( xData, yData, ft );
            a = fitresult.a;
            
            % power law of all the points, given a (get b)
            [xData, yData] = prepareCurveData( log(IM), log(EDP)-log(a) );
            ft = fittype( {'x'} );
            [fitresult, gof]  = fit( xData, yData, ft );
            b = fitresult.a;
        end
        
        resid = log(EDP) - (log(a) + b.*log(IM));
        sigma = std( resid );
        
        powerLaw = [a b];

        PSDMformula.mu = @(IM,a,b,sigma) ( a.*exp(sigma).*IM.^b);
        PSDMformula.IM = @(mu,a,b,sigma) ( (mu./(a*exp(sigma))).^(1/b) );
        
        fragFormula.eta = @(muCap,a) (a.*IM.^b);
        fragFormula.beta = @(sigma) (sigma);

        

%             % Control plot
%              y = linspace(0,max(EDP),100);
%              x = linspace(0,max(IM),100);
% 
%             figure; hold on
%             scatter(IM, EDP)
%             plot(x,PSDMformula.mu(x, a, b, 0), 'k', ...
%                 'LineWidth', 3)
%             plot(x,PSDMformula.mu(x, a, b, sigma), '--k', ...
%                 'LineWidth', 1.5,'Color',[0.5,0.5,0.5])
%             plot(x,PSDMformula.mu(x, a, b, -sigma), '--k', ...
%                 'LineWidth', 1.5,'Color',[0.5,0.5,0.5])
% %             set(gca,'YTick',[1 200],'XTick',[1 30])
%             xlabel('IM')
%             ylabel('EDP')       

%             figure; hold on
%             scatter(IM, EDP)
%             plot(PSDMformula.IM(y, a, b, 0), y, 'k', ...
%                 'LineWidth', 3)
%             plot(PSDMformula.IM(y, a, b, sigma), y, '--k', ...
%                 'LineWidth', 1.5,'Color',[0.5,0.5,0.5])
%             plot(PSDMformula.IM(y, a, b, -sigma), y, '--k', ...
%                 'LineWidth', 1.5,'Color',[0.5,0.5,0.5])
%             set(gca,'YTick',[1 200],'XTick',[1 30])
%             xlabel('IM')
%             ylabel('EDP')    
%             xlim([0 3])
%         
        end


        function [a, sigma, b, NpointsReg, PSDMformula, fragFormula] = ...
                fitPSDMbilinear(IM, R, EDP, EDPds4)
            %fitPSDM fits the law MU = a.*(R-1) + 1; R = SA/SAy
            % if you want different fitting methods, copy from "comparePSDMs"
            
            
            % the 1.05 limit (instead of 1) avoids logs too close to zero
            
            toInclude = EDP > 1.05 & R > 1 & EDP <= EDPds4;
            NpointsReg = sum(toInclude);
            
            b = 1;
            lna = mean( log(EDP(toInclude)-1) - b*log(IM(toInclude)-min(IM(toInclude))*0.9) );
            % replaced 1 for min(IM)*1.1, since the IM at yielding is not
            % clearly defined for other IMs other than R. Assumes that the
            % training database starts at a IM level 10% higher than the
            % yield IM

            a = exp(lna);

            resid = log(EDP(toInclude)-1) - (lna + b.*log(IM(toInclude)-min(IM(toInclude))*0.9));
            sigma = std( resid );

            PSDMformula.mu = @(R,a,sigma) ( a.*exp(sigma).*(R-1) + 1 );
            PSDMformula.R = @(mu,a,sigma) ( (mu-1)./(a.*exp(sigma)) + 1 );
            
            fragFormula.eta = @(muCap,a) ((muCap-1)/a + 1);
            fragFormula.beta = @(sigma) (sigma);
            
            % Control plot
%             Mu=300;
%             dummyMU = 1:0.01:Mu;
%             figure; hold on
%             title('Ductility-based PSDM')
%             grid(gca,'minor');
%             scatter(IM(toInclude), EDP(toInclude),36,'MarkerEdgeColor',[128 128 128]/265,'MarkerFaceColor',[128 128 128]/285)
%             plot([0 PSDMformula.R(dummyMU, a, 0)], [0 dummyMU], 'k', ...
%                 'LineWidth', 3)
%             plot([0 PSDMformula.R(dummyMU, a, sigma)], ...
%                 [0 dummyMU], '--k', 'LineWidth', 1.5,'Color',[0.5,0.5,0.5])
%             plot([0 PSDMformula.R(dummyMU, a, -sigma)], ...
%                 [0 dummyMU], '--k', 'LineWidth', 1.5,'Color',[0.5,0.5,0.5])
%             set(gca,'YTick',[1 200],'XTick',[1 30])
%             xlim([0, 30])
%             xlabel('R [-]')
%             ylabel('\mu [-]')
%             set(gca, 'FontSize', 25)
%             
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
