classdef poissyFit< matlab.mixin.Copyable
    % A class to estimate Poisson rates (lambda) from measured Fluorescence.
    %
    % The theory behind this code is described here:
    % Ganmor, E., Krumin, M., Rossi, L. F., Carandini, M. & Simoncelli, E. P.
    % Direct Estimation of Firing Rates from Calcium Imaging Data. 1–34 (2016).
    % http://arxiv.org/abs/1601.00364
    %
    % My implementation started from the code provided here:
    % https://github.com/eladganmor/Imaging_Analysis
    % and added some code checking, bootstrapping, code comments, and
    % additinional functionality to fit double von Mises, do splithalves
    % consistency checks and comparison with deconvolved spike estimates.
    %
    % See README.md , demos/parammetric, demos/simple, demos/twoVonMises
    %
    % BK - May 2023.
    properties (SetAccess=public,GetAccess=public)
        nrWorkers       = 0; % Set to >0 to use parfor for bootstrap resampling
        
        spikeCountDistribution (1,1) string = "POISSON";% Assumed spike count distribution. POISSON or EXPONENTIAL
        
         measurementNoise  = 1;      % Standard deviation of Fluorescence measurement

        tuningFunction  = []  % Tuning function for direct parametric estimation.
        % Optimization parameters.
        options =    optimoptions(@fminunc,'Algorithm','trust-region', ...
            'MaxIter',1e5, ...
            'SpecifyObjectiveGradient',true, ...
            'display','final', ...
            'CheckGradients',true, ... % Set to true to check supplied gradients against finite differences
            'diagnostics','off', ...
            'FiniteDifferenceType','central'); % Central is slower but more accurate

        hasDerivatives         % Integer indicating whether the tuning function has derivatives (1) and or second derivatives (2)
        % This is normally set automatically, but if you use an anonymous function, you may
        % have to overrule the automated setting (See demos/parametric)
    end

    properties (SetAccess=protected)
        % Matrices that are constant for a given problem and computed only
        % once in initialize:
        nrSpikes                % [nrSpikesMax 1] - Spike counts over which the algorithm marginalizes.
        pFGivenSpikes           % [nrSpikesMax nrTimePoints nrTrials] - Probability of the observed F given the spikes
        nrSpikesFactorial       % nrSpikes!
        nrSpikesPerTimepoint    % [nrSpikesMax nrTimePoints nrTrials] - Matrix with nrSpikes in each column
        DM                      % [nrPredictors nrTimePoints nrTrials] - Design matrix
        fTilde                  % [nrTimePoints nrTrials]
        % The raw data
        fluorescence            % [nrTimePoints nrTrials]
        stimulus                % [1 nrTrials]
        binWidth                % Time bins for the fluorescence signal in seconds.

        uStimulus               % Unique stimulus values
        stimulusIx              % Index into uStimulus for each trial
        tuningCurve             % Non-parametric estimate of the tuning curve (F per second)
        tuningCurveErr          % IQR for tuning curve
        bestGuess               % Initial guess of the parameters (from nonparametric tuning)
        blank                   % Response to blank (coded as NaN in the stimulus)
        blankErr                %
        gof                     % Goodness of fit (Spearman correlation of the estimate with the nonparameteric tuning curve)
        residualZ               % Z score of the mean residuals (0 = good; no bias)
        residualNoise           % Stdev of the residuals relative to measurement noise  (1 =good given the noise)
        parms                   % Estimated parameters that capture the fluorescence best
        parmsError              % Errors on the parameter estimates.
        bootstrap               % Number of bootstrap sets used to estimate parms/parmsError.
        % If this is 1, the parms are the result of
        % the ML fit, if it is >1 they are the
        % mean of the bootstrap sets (and the error
        % the standard deviation across sets).
        bootParms               % Estimates for each bootstrap set.
        bootParmsError           % Estimates for each bootstrap set.
        bootGof                   % Estimates for each bootstrap set.

        exitFlag                % Exitflag of the optimization routine.


        % Model parameters  - Set in constructor
        maxRate             %  The maximum number of spikes (per second) that a neuron can reasonable fire under the conditions of the experiment
        tau                 % Fluorescence decay time constant [s]
        fPerSpike            % Fluorescence per spike. This is a true free parameter
     
    end

    properties (Dependent)
        nrTimePoints        % Time points per trial
        nrTrials            % Number of trials
        nrObservations      % Total number of Fluorescence observations (nrTimePoints*nrTrials)
        nrSpikesMax         % Maximum number of spikes in a single bin (derived from maxRate)
        nrParms             % Nunmber of estimated parameters
    end

    %% Get methods for dependent members
    methods
        function v = get.nrTimePoints(o)
            v = size(o.fluorescence,1);
        end
        function v = get.nrObservations(o)
            v = size(o.fluorescence);
        end

        function v = get.nrTrials(o)
            v = size(o.fluorescence,2);
        end

        function v = get.nrSpikesMax(o)
            v = ceil(o.maxRate*o.binWidth);
        end

        function v= get.nrParms(o)
            if isempty(o.tuningFunction)
                v = size(o.DM,1);
            else
                v = numel(o.bestGuess);
            end
        end

    end

    %% Public functions
    methods (Access=public)

        function o = poissyFit(stim,fluor,dt,fun,pv)
            % Construct an object, based on the stimulus,fluorescence and
            % timebins of the data.
            arguments
                stim (1,:) double   % The stimulus presented in each trial
                fluor (:,:) double  % The response in each bin (rows) and each trial (rows)
                dt (1,1) double     % The duration of a single bin (in seconds).
                fun = []            % Tuning function for parametric fits.
                pv.hasDerivatives (1,1) double =NaN % Specify how many derivatives the fun has (needed for anonymous fun)
                pv.scaleToMax (1,1) logical = false % Scale fluorescene to match the o.maxRate
                pv.maxRate     (1,1) double = 100  % The maximum spike rate that is considered.
                pv.fPerSpike (1,1) double = 80 % Each spike produces this much Ca/F
                pv.tau  (1,1) double =1.3       % Exponential decay time of Ca/F
            end
            o.maxRate  = pv.maxRate;
            o.fPerSpike= pv.fPerSpike;
            o.tau= pv.tau;
            o.binWidth = dt;
            o.stimulus = stim;
            
            if pv.scaleToMax
                peakF = max(fluor,[],"all");
                maxModeledF = o.nrSpikesMax*o.fPerSpike*(1+exp(-o.binWidth/o.tau));  
                scale = peakF/maxModeledF;         
                fluor = fluor/scale; % Scale fluorescence         
            end          
            o.fluorescence = fluor;
           
            
            o.tuningFunction = fun;           
            % Check whether the fun has deriviatives. This fails if the fun
            % is an anonymous function - the user then has to set
            % .hasDerivatives manually.
            if isempty(fun)
                o.hasDerivatives  = 2;
            elseif ~isnan(pv.hasDerivatives)
                o.hasDerivatives = pv.hasDerivatives;
            elseif nargout(fun) <0
                warning('Please set o.hasDerivatives manually for this (anonymous?) tuning function')
            else
                o.hasDerivatives  = nargout(fun)-1;
            end
        end

        function v = copyWithNewData(o,stim,fluor,dt,fun,pv)
            % Construct a new object with the same parameters as an
            % existing object, but new stimuli ,fluorescence, and
            % (optionally) new timebin and tuning function. Used by
            % the bootstrap procedure in solve.
            arguments
                o (1,1) poissyFit
                stim (1,:) double   % The stimulus presented in each trial
                fluor (:,:) double  % The response in each bin (rows) and each trial (rows)
                dt (1,1) double =o.timeBin     % The duration of a single bin (in seconds).
                fun = o.tuningFunction  % Tuning function for parametric fits.
                pv.maxRate     (1,1) double = o.maxRate % The maximum spike rate that is considered.
                pv.fPerSpike (1,1) double = o.fPerSpike % Each spike produces this much Ca/F
                pv.tau  (1,1) double = o.tau       % Exponential decay time of Ca/F
            end
            v = copy(o); % Make a shallow copy
            v.maxRate  = pv.maxRate;
            v.fPerSpike= pv.fPerSpike;
            v.tau= pv.tau;
            
            % Replace the 'data'
            v.stimulus = stim;
            v.fluorescence = fluor;
            v.binWidth = dt;
            v.tuningFunction = fun;
            v.initialize
        end

            function [lambda,dLogLambda,d2LogLambda]  = lambdaFromStimulus(o,parms)
            % Given the current parameter estimates, return the underlying lambda
            % (i.e. the poisson rate) for each time bin in the data.
            if isempty(o.tuningFunction)
                lambda = squeeze(exp(tensorprod(parms,o.DM,2,1)));
                dLogLambda = o.DM;
                d2LogLambda = zeros([numel(parms),numel(parms),o.nrObservations]);
            else
                if o.hasDerivatives ==2
                    [logLambda, dLogLambda, d2LogLambda] = o.tuningFunction(o.stimulus,parms);
                    dLogLambda = permute(repmat(dLogLambda,[1 1 o.nrTimePoints]),[1 3 2]);
                    d2LogLambda = permute(repmat(d2LogLambda,[1 1 1 o.nrTimePoints]),[1 2 4 3]);
                elseif o.hasDerivatives ==1
                    [logLambda, dLogLambda] = o.tuningFunction(o.stimulus,parms);
                    dLogLambda = permute(repmat(dLogLambda,[1 1 o.nrTimePoints]),[1 3 2]);
                else
                    logLambda = o.tuningFunction(o.stimulus,parms);
                end
                lambda = exp(logLambda);
                lambda = repmat(lambda,[o.nrTimePoints 1]);
            end
        end

        function plot(o,pv)
            arguments
                o (1,1) poissyFit
                pv.showErrorbars  = false
                pv.showBootstrapSets = false
            end
            % Plot an overview of the data and tuning curve fits/
            if isempty(o.tuningFunction)
                predictedTuningCurve= exp(o.parms(2:end));% + o.parms(1));
                predictedTuningCurveHigh= exp(o.parms(2:end)+o.parmsError(2:end)  + o.parms(1) + o.parmsError(1));
                predictedTuningCurveLow= exp(o.parms(2:end)-o.parmsError(2:end)  + o.parms(1) - o.parmsError(1));
            else
                predictedTuningCurve = exp(o.tuningFunction(o.uStimulus,o.parms));
                predictedTuningCurveHigh = exp(o.tuningFunction(o.uStimulus,o.parms+o.parmsError));
                predictedTuningCurveLow= exp(o.tuningFunction(o.uStimulus,o.parms-o.parmsError));
            end
            pos = predictedTuningCurveHigh- predictedTuningCurve;
            neg = predictedTuningCurve - predictedTuningCurveLow;

            yyaxis left
            % Show the non-parametric tuning curve (already in /s)
            h1= errorbar(o.uStimulus,o.tuningCurve,o.tuningCurveErr/2,'b');
            ylabel 'Fluorescence (/s)'
            yyaxis right
            hold on
            if pv.showErrorbars
                h2= errorbar(o.uStimulus,1./o.binWidth*predictedTuningCurve,1./o.binWidth*pos,1./o.binWidth.*neg,'r');
            else
                h2= plot(o.uStimulus,1./o.binWidth*predictedTuningCurve,'r');
            end
            if pv.showBootstrapSets
                if isempty(o.tuningFunction)
                    bsPredictedTuningCurve= exp(o.bootParms(:,2:end));% + o.parms(1));                
                else
                    bsPredictedTuningCurve = nan(o.bootstrap,numel(predictedTuningCurve));
                    for i=1:o.bootstrap
                        bsPredictedTuningCurve(i,:) = exp(o.tuningFunction(o.uStimulus,o.bootParms(i,:)))';
                    end                
                end
                plot(o.uStimulus,1./o.binWidth*bsPredictedTuningCurve,'Color',[0.8 0.8 0.8])               
            end
            ylabel 'Lambda (spk/s)'
            xlabel 'Stimulus'
            title (['Parms: ' num2str(o.parms,2)])
            legend([h1(1) h2(1)],'Fluorescence with IQR','Rate with SE')

        end



        function [r,bootParms,rSpikes,spikeBootParms,rCross] = splitHalves(o,nrBoot,guess,spikes,pv)
            % Asess performance by comparing parameter estimates based on
            % split halves of the data.
            %
            % If the estimation of the rate is successful then identical
            % repeats of the same condition (e.g., stimulus orientation)
            % should result in identical estimates of the underlying rate.
            % Therefore the correlation between parameters based on one half of
            % trials and comparing it with the other half of trials
            % quantifies how well the algorithm works.
            % Even a perfect estimator would have a correlation less than 1
            % because of Poisson noise on the spike counts.
            % See Pachitariu et al. J Neurosci 2018 38(37):7976 –7985.
            %
            % For direct, non-parametric estimates the correlation is
            % calculated directly between the estimates based on the split
            % halves.
            % For parametric estimates, we estimate the parameters of the
            % tuning function and from those the predicted response at each
            % of the conditions, and correlate these tuning functions based on
            % the split halves.
            %
            % To compare poissyFits performance with a sequential estimate
            % based on spike deconvolution, pass the deconvolved spike
            % count estimate as 'spikes'. This matrix should match the fluorescence
            % matrix [nrTimePoints nrTrials] used by poissyFit.
            %

            arguments 
                o (1,1) poissyFit
                nrBoot  (1,1) double {mustBeInteger,mustBePositive}  =250  % The number of randomly sampled split halves to use
                guess = []  % The initial guess for parameter estimation ([] uses the internal guess).
                spikes (:,:) double = [];
                pv.corrFun = @(x,y)corr(x,y,'Type','Spearman')
            end
%             arguments (Output)
%                 r (1,1) double  % Mean of the split-half correlation acros bootstrap sets
%                 bootParms double  % Parameter estimates in each bootstrap set
%                 rSpikes (1,1) double  % Mean of the split-half correlation in the deconvolved spikes
%                 spikeBootParms double % Parameter estimates based one the spikes in each bootstrap set
%                 rCross (1,1) double % Mean of the correlation between fluorescence and spike estimates
%             end
            initialize(o);
            % Avoid warnings on display during bootstrapping
            o.options.Display = 'off';
            o.options.Diagnostics = 'off';
            o.options.CheckGradients = false;
            bootParmsOne   = nan(o.nrParms,nrBoot);
            bootParmsOther = nan(o.nrParms,nrBoot);
            r =nan(nrBoot,1);
            if ~isempty(spikes)
                % Compare with results based on deconvolved spikes
                % (e.g. obtained from some other package like suite2p or caiman that
                % does the deconvolution)
                % We can use the same fitting methods by creating a
                % poissyFit object in which ca/F decays instantly
                % (tau=0) and each spike produces 1 unit of ca/F.
                % We also assume that the measurement noise on the spikes
                % is a scaled version of the measurement noise on the
                % fluorescence.
                oSpk = o.copyWithNewData(o.stimulus,spikes,o.binWidth,o.tuningFunction);
                oSpk.tau =0;
                oSpk.fPerSpike =1;
                % Scale the measurement noise by the nonparametric estimate
                % of the noise in the tuning curves
                oSpk.measurementNoise = mean(oSpk.tuningCurveErr)./mean(o.tuningCurveErr).*o.measurementNoise;
                % Preallocate
                spikeBootParmsOne = nan(o.nrParms,nrBoot);
                spikeBootParmsOther= nan(o.nrParms,nrBoot);
                rSpikes =nan(nrBoot,1);
                rCross = nan(nrBoot,1);
            else
                oSpk = []; % Define to avoid parfor complaints
            end
            % Bootstrap fitting on split halves
            %parfor (i=1:nrBoot, o.nrWorkers)
            for i=1:nrBoot % Uncomment for debugging
                [oneHalfTrials,otherHalfTrials] =resampleTrials(o,false,0.5) ;
                thisO = o.copyWithNewData(o.stimulus(:,oneHalfTrials),o.fluorescence(:,oneHalfTrials),o.binWidth,o.tuningFunction);
                solve(thisO,1,guess);
                bootParmsOne(:,i) = thisO.parms';

                thisO = o.copyWithNewData(o.stimulus(:,otherHalfTrials),o.fluorescence(:,otherHalfTrials),o.binWidth,o.tuningFunction);
                solve(thisO,1,guess);
                bootParmsOther(:,i) = thisO.parms';


                if ~isempty(spikes)
                    % Use the same trials to estimat spike tuning from the
                    % oSpk object
                    thisO = oSpk.copyWithNewData(oSpk.stimulus(:,oneHalfTrials),oSpk.fluorescence(:,oneHalfTrials),o.binWidth,o.tuningFunction); %#ok<PFBNS>
                    solve(thisO,1,guess);
                    spikeBootParmsOne(:,i) = thisO.parms';

                    thisO = oSpk.copyWithNewData(oSpk.stimulus(:,otherHalfTrials),oSpk.fluorescence(:,otherHalfTrials),o.binWidth,o.tuningFunction);
                    solve(thisO,1,guess);
                    spikeBootParmsOther(:,i) = thisO.parms';
                end
            end


            %% Determime correlations
            
            for i=1:nrBoot
                if isempty(o.tuningFunction)
                    % Ignore the grand mean, correlate the rest
                    r(i)= pv.corrFun(bootParmsOne(2:end,i),bootParmsOther(2:end,i));
                    if ~isempty(spikes)
                        rSpikes(i)= pv.corrFun(spikeBootParmsOne(2:end,i),spikeBootParmsOther(2:end,i));
                        rCross(i) = (pv.corrFun(spikeBootParmsOne(2:end,i),bootParmsOne(2:end,i)) +pv.corrFun(spikeBootParmsOther(2:end,i),bootParmsOther(2:end,i)))/2;
                    end
                else
                    % Predict the tuning and correlate those
                    tc1= o.tuningFunction(o.uStimulus,bootParmsOne(:,i))';
                    tc2 = o.tuningFunction(o.uStimulus,bootParmsOther(:,i))';
                    r(i) = corr(tc1,tc2);
                    if ~isempty(spikes)
                        tcSpk1= o.tuningFunction(o.uStimulus,spikeBootParmsOne(:,i))';
                        tcSpk2 = o.tuningFunction(o.uStimulus,spikeBootParmsOther(:,i))';
                        rSpikes(i)= pv.corrFun(tcSpk1,tcSpk2);
                        rCross(i) = (pv.corrFun(tc1,tcSpk1) +pv.corrFun(tc2,tcSpk2))/2;
                    end
                end
            end
            r = mean(r);
            if ~isempty(spikes)
                rSpikes= mean(rSpikes);
            end
            if nargout>1
                bootParms = cat(3,bootParmsOne,bootParmsOther);
            end
            if nargout >3
                spikeBootParms = cat(3,spikeBootParmsOne,spikeBootParmsOther);
            end
            if nargout >4
                rCross = mean(rCross);
            end

        end

        function v = estimateNoise(o)
            % This estimates the noise in the fluorescence measurements by
            % finding timepoints that are identical and have an identical previous timepoint
            % , then determining the std across those identical repeats ,
            % and then the mean across those std over conditions.
            % This is an estimate of the noise introduced by the Poisson
            % nature of the spiking plus the measurement noise. So it
            % cannot simply be used as the measurement noise, but it may
            % help to determine whether some independent measure for
            % measurement noise is reasonable. 
           if o.tau>0               
            state = cat(1,o.DM(:,2:end,:),o.DM(:,1:end-1,:));
            state = reshape(state,size(o.DM,1)*2,(o.nrTimePoints-1)*o.nrTrials)';
           else
               % No history ; state is a single time bin
               state = o.DM(:,2:end,:);              
               state = reshape(state,size(o.DM,1),(o.nrTimePoints-1)*o.nrTrials)';
           end
            [~,~,ix] = unique(state,'rows');
            f = o.fluorescence(2:end,:);
            nrSpk = f./o.fPerSpike;
            % meanPerState = accumarray(ix,nrSpk(:)',[],@mean);
            stdPerState = accumarray(ix,nrSpk(:)',[],@std);
            v  = mean(stdPerState);

        end


        function solve(o,boot,guess)
            % solve(o,boot,guess)
            % Solve the fitting problem. Optionally do this in bootstrap
            % resampling fashion.
            % boot - Number of boostrap sets
            % guess - Optional initial guess of the parameters of the fit.
            %        If left unspecified, it uses a non-parametric estimate
            %        of the tuning cruve.
            % Each boostrap estimate is stored in bootstrapParms.
            % If bootstrapping is used, the .parm and .parmError are the
            % bootstrap estimates, otherwise they are the outcome of the
            % (single) optimization procedure with errors based on the
            % Hessian at the solution.
            arguments
                o (1,1)
                boot (1,1) double {mustBeInteger,mustBePositive} = 1  % Number of bootstrap sets
                guess = [];
            end

            if o.hasDerivatives<1 && o.options.SpecifyObjectiveGradient
                warning('The tuning function does not provide derivatives; switching to quasi-newton')
                o.options.SpecifyObjectiveGradient = false;
                o.options.Algorithm = 'quasi-newton';
            end
            o.initialize
            if isempty(guess)
                guess  = o.bestGuess;
            end
            % Solve the  minimization problem.
            [o.parms,ll,o.exitFlag,output,derivative,hessian] = fminunc(@o.logLikelihood,guess,o.options); %#ok<ASGLU> 
            % The Hessian is usually not very informative.
             state = warning('query');
             warning('off','MATLAB:nearlySingularMatrix')
             iH = diag(inv(hessian));
             iH(iH<0) = NaN; % Negative element means we're not at a (local) minimum. Probably a flat piece of the LL.
             o.parmsError = sqrt(iH)';           
             warning(state)
            fixup(o,"SOLVE");
            computeQuality(o)
            
            % Perform bootstrap resampling to get a robust estimate of the
            % parms and their errors
            if boot>1
                tmpBootParms = nan(boot,o.nrParms);
                tmpBootParmsError = nan(boot,o.nrParms);
                tmpBootGof = nan(boot,1);
                %parfor (i=1:boot,o.nrWorkers)
                    for i=1:boot % For debugging
                    thisTrials = resampleTrials(o,true,1);
                    thisO = o.copyWithNewData(o.stimulus(:,thisTrials),o.fluorescence(:,thisTrials),o.binWidth,o.tuningFunction);                   
                    solve(thisO,1,guess);
                    tmpBootParms(i,:) = thisO.parms;
                    tmpBootParmsError(i,:) = thisO.parmsError;
                    tmpBootGof(i) = thisO.gof;
                end
                % Store the mean and std across sets as the outcome
                % (This is corrected in fixup for specific tuning
                % functions)
                o.parms      = mean(tmpBootParms,1,'omitnan');
                o.parmsError = std(tmpBootParms,0,1,'omitnan');
                o.gof = mean(tmpBootGof,1,"omitnan");
                o.bootstrap  = boot;
                o.bootParms = tmpBootParms;
                o.bootParmsError = tmpBootParmsError;
                o.bootGof = tmpBootGof;
                fixup(o,"BOOTSTRAP");
                % Recompute residuals based on bs estimate
                computeQuality(o);
            else
                o.bootParms = [];
                o.bootParmsError =[];
                o.bootstrap = 0;
            end

           
        end


    end



    %% Protected, internal functions
    methods  (Access=protected)
        function computeQuality(o)

            if isempty(o.tuningFunction)
                % Grand mean plus tuning
                tc = exp(o.parms(2:end)')+o.parms(1);                                                
            else
                % Predict the tuning 
                tc= exp(o.tuningFunction(o.uStimulus,o.parms)');
            end
            o.gof = corr(tc,o.tuningCurve,'type','Spearman');

            
            residuals = tc(o.stimulusIx)'-o.fluorescence;
            meanR = mean(residuals,"all","omitnan");
            stdR =  std(residuals,0,"all","omitnan");
            o.residualZ = abs(meanR./stdR);
            o.residualNoise = stdR./o.measurementNoise;
        end

        function initialize(o)
            % Pre-compute some quantities that do not depend on the
            % estimated parameters to save time during optimization
            % Becasue the first bin in each trial has no preceding bin,
            % we cannot estimate lambda in that bin hence these matrices
            % correspond to timepoint 2 and later.
            % The logLikelihood computation will combine these
            % matrix with the lambda estimates from time points 2 and
            % later.

            switch (o.spikeCountDistribution)
                case "POISSON"
                    o.nrSpikes             = (0:o.nrSpikesMax)';
                    o.nrSpikesFactorial    = repmat(factorial(o.nrSpikes),[1 o.nrTimePoints-1 o.nrTrials]);
                    o.nrSpikesPerTimepoint = repmat(o.nrSpikes,[1 o.nrTimePoints-1 o.nrTrials]);
                    % The F is determined  by the previous time bin's F (which decays at a rate of tau)
                    % plus the new F generated by the spikes in the current bin
                    % F is  [nrTimePoints nrTrials]. The special case of
                    % o.tau ==0 (and o.fPerSpike =1) can be used to fit
                    % directly to estimate spike counts in independent bins
                    % (i.e. without Ca decay effects from the previous time
                    % bin)
                    F = permute(repmat(o.fluorescence,[1 1 o.nrSpikesMax+1]),[3 1 2]);
                    fExpected = o.fPerSpike.*o.nrSpikesPerTimepoint + F(:,1:end-1,:).*exp(-o.binWidth/o.tau);
                    % The probability of observing an F given the spikes is a
                    % normal distribution with a stdev equal ot the measurement
                    % noise. Adding eps to avoid underflow.
                    o.pFGivenSpikes = eps+normpdf(F(:,2:end,:) - fExpected,0,o.measurementNoise);
                case "EXPONENTIAL"
                    o.fTilde =  1/o.fPerSpike*(o.fluorescence(1:end-1,:).*exp(-o.binWidth/o.tau)-o.fluorescence(2:end,:));

            end

            % Construct the design matrix
            o.DM = designMatrix(o);

            % Determine a fluorescence tuning curve (mean rate per unique stimulus
            % value)
            mResponse = mean(o.fluorescence,1,"omitnan"); % Average over time bins
            isBlank  = isnan(o.stimulus);
            o.blank = median(mResponse(isBlank));
            o.blankErr = diff(prctile(mResponse(isBlank),[25 75]),1);
            [o.uStimulus,~,o.stimulusIx] = unique(o.stimulus(~isBlank));
            o.tuningCurve = accumarray(o.stimulusIx,mResponse(~isBlank)',[],@mean)/o.binWidth;
            o.tuningCurveErr = accumarray(o.stimulusIx,mResponse(~isBlank)',[],@std)/o.binWidth;
            

            % Heuristic sanity checks
            peakF = prctile(o.fluorescence(:),99);
            % The maximum fluorescence that could occur is two successive
            % bins with the maxmimum number of modeled spikes:
            maxModeledF = o.nrSpikesMax*o.fPerSpike*(1+exp(-o.binWidth/o.tau));  
            newFPerSpike = peakF./(o.nrSpikesMax*(1+exp(-o.binWidth/o.tau)));
            if peakF > maxModeledF
                fprintf(2,'The peak Fluorescence (%.3f) is too high for the modeled maximum spike rate (%d). \n \t Scale the input Fluorescence, or increase fPerSpike to at least %.2f\n.',peakF,o.maxRate,newFPerSpike);
            end

            fixup(o,"BESTGUESS");
        end


        function fixup(o,what)
            % Various fixups needed for tuningFunctions. All put in one
            % place to make future additions of additional tuningFunctions
            % a bit easier
            arguments
                o (1,1) poissyFit
                what (1,1) string {mustBeMember(what,["BESTGUESS","SOLVE","BOOTSTRAP"])}
            end

            switch (what)
                case "BESTGUESS"
                    if isempty(o.tuningFunction)
                        % Because the LL is convex for direct estimation, the
                        % starting guess does not matter.
                        % 0 baseline, median(F)/FPerSpike
                        o.bestGuess = zeros(1,size(o.DM,1));
                    elseif contains((func2str(o.tuningFunction)),'logVonMises')
                        % Start at what looks like the preferred
                        % orientation.
                        [~,prefIx] = max(o.tuningCurve);
                        o.bestGuess = [0 o.uStimulus(prefIx) 1];
                    elseif contains((func2str(o.tuningFunction)),'logTwoVonMises')
                        % Start at what looks like the preferred
                        % orientation.
                        [~,prefIx] = max(o.tuningCurve);
                        o.bestGuess = [0 o.uStimulus(prefIx) 0 0 0];
                    else
                        errror('PLease add a %s initialization in poissyFit/fixup',what);
                    end

                case "SOLVE"
                    % Post processing for consistent parameter interpretation
                    if isempty(o.tuningFunction)
                        % Nothing to do
                    elseif contains((func2str(o.tuningFunction)),'logVonMises')
                          % Cannot wrap as this could be 180 or 360
                          % periodicity
                    elseif contains((func2str(o.tuningFunction)),'logTwoVonMises')
                        if o.parms(5)>o.parms(4)
                            % Either amp1 or amp2 could be the preferred direction. For consistency
                            % make sure amp1 is the largest one and flip
                            % preferred/anti-preferred if needed
                            o.parms = o.parms([1 2 3 5 4]);% Swap amplitude
                            o.parms(2) = mod(o.parms(2)+180,360);  % Redefine preferred
                        else
                            o.parms(2) = mod(o.parms(2),360);% Convenience
                        end
                    else
                        errror('PLease add a %s initialization in poissyFit/fixup',what);
                    end
                case 'BOOTSTRAP'
                    % Determine errors from bootstrap parms
                    if isempty(o.tuningFunction)
                        % Nothing to do
                    elseif contains((func2str(o.tuningFunction)),'logVonMises') || ...
                            contains((func2str(o.tuningFunction)),'logTwoVonMises')
                        % Estimate distributions are highly skewed use
                        % median
                        o.parms      = median(o.bootParms,1,'omitnan');
                        iqr = @(x) (diff(prctile(x,[25 75]),1));
                        o.parmsError = iqr(o.bootParms);
                
                        % Determine mean angle and circular standard deviation                        
                        z = exp(pi/180*1i*o.bootParms(:,2));
                        R = mean(z);
                        o.parmsError(2) = 180/pi*sqrt(-2*log(abs(R))); % Circular standard deviation
                        o.parms(2) = 180/pi*angle(R);
                    else
                        errror('PLease add a %s initialization in poissyFit/fixup',what);
                    end
            end
        end

        function DM = designMatrix(o)
            % Given a row vector of stimulus values (e.g. orientation per
            % trial),  construct a design matrix where each row corresponds
            % to a stimulus (orientation) , each column an observation (time bin: nrTimePoints*nrTrials)
            % and the value 1 indicates that that stimulus (orientation) was on the screen
            % during that time bin.
            [uX,~,xIx] = unique(o.stimulus);
            nrX = numel(uX);
            xIx = repmat(xIx(:)',[o.nrTimePoints 1]);

            DM = zeros([nrX,o.nrObservations]);
            timeIx = repmat((1:o.nrTimePoints)',[o.nrTrials 1]);
            trialIx = repmat((1:o.nrTrials),[o.nrTimePoints 1]);
            ix = sub2ind(size(DM),xIx(:),timeIx(:), trialIx(:));
            DM(ix) = 1;
            if isempty(o.tuningFunction)
                % Add a constant term.
                DM = [ones([1,o.nrObservations]) ;DM]; % Constant term
            end

        end


        function [inSet,notInSet] = resampleTrials(o,withReplacement,frac)
            % Resample trials to use in bootstrapping. This resampling
            % makes sure to resample trials from each of the unique
            % stimulus conditions so that the resampled trials have all of
            % the conditions.
            %
            % Input
            % withReplacement - set to true to sample with replacement
            % (used by bootstrapping).
            % frac  - The fraction of trials to resample. 1 means resample
            %           all, 0.5 with replacement=false means split halves.
            % OUTPUT
            % inSet - The set of selected trials
            % outSet - The trials not in the set.
            %
            arguments
                o (1,1) poissyFit
                withReplacement (1,1) logical = false
                frac (1,1) double {mustBeInRange(frac,0,1)} =  1
            end
            stimCntr= 1;
            trialsPerStim = cell(1,numel(o.uStimulus));
            for u= o.uStimulus
                trialsPerStim{stimCntr} = find(o.stimulus==u);
                stimCntr = stimCntr+1;
            end
            if withReplacement
                % Sampling with replacement
                inSet = cellfun(@(x) (x(randi(numel(x),[1 ceil(frac*numel(x))]))),trialsPerStim,'uni',false);
            else
                % Sampling without replacement (e.g. to split 80/20 or
                % split halves)
                inSet = cellfun(@(x) (x(randperm(numel(x),ceil(frac*numel(x))))),trialsPerStim,'uni',false);
            end
            notInSet = cellfun(@(x,y) setxor(x,y),trialsPerStim,inSet,'uni',false);
            inSet = cat(2,inSet{:});
            notInSet =cat(2,notInSet{:});
        end


        function [ll, derivative, hessian] = logLikelihood(o,parms)
            % Calculates the log likelihood given the data and parameters
            % returns negative of log likelihood, plus the analytic
            % derivatives, Hessian (if available and requested)
            switch o.hasDerivatives
                case 2
                    [lambda,dLogLambda,d2LogLambda] = lambdaFromStimulus(o,parms);
                case 1
                    [lambda,dLogLambda] = lambdaFromStimulus(o,parms);
                otherwise
                    lambda = lambdaFromStimulus(o,parms);
            end
            switch  o.spikeCountDistribution
                case "POISSON"
                    % Assume the spike count distribution is poisson (with
                    % a cut off of o.nrSpikesMax per bin).
                    lambda =reshape(lambda(2:end,:),[1 o.nrTimePoints-1 o.nrTrials]);
                    LAMBDA= repmat(lambda,[o.nrSpikesMax+1 1 1]);
                    pX = eps+ exp(-LAMBDA).*LAMBDA.^o.nrSpikesPerTimepoint./o.nrSpikesFactorial;
                    pXY = pX.*o.pFGivenSpikes;
                    sum_PXY = sum(pXY,1); % Sum over nrSpikes
                    ll = -sum(log(sum_PXY),[2 3]); % Sum over nrTimePoints and nrTrials

                    % Derivatives, if requested and available.
                    if o.hasDerivatives>0 && nargout>1
                        sumX_PXY = tensorprod(o.nrSpikes,pXY,1,1);
                        sumX_PXY_over_sum_PXY = sumX_PXY./sum_PXY;
                        gamma = sumX_PXY_over_sum_PXY - lambda;
                        derivative = -tensorprod(dLogLambda(:,2:end,:),gamma,[2 3],[2 3]);
                    end
                    if o.hasDerivatives>1 && nargout>2
                        % Hessian
                        c = lambda- tensorprod(o.nrSpikes.^2,pXY./sum_PXY,1,1) + (sumX_PXY_over_sum_PXY).^2;
                        part1  = tensorprod((c.*dLogLambda(:,2:end,:)),dLogLambda(:,2:end,:),[2 3],[2 3]);
                        part2 = -sum(d2LogLambda(:,:,2:end,:).*reshape(gamma,1,1,o.nrTimePoints-1,o.nrTrials),[3,4]);
                        hessian = part1 + part2;
                    end
                case "EXPONENTIAL"
                    % Assume the spike "count" distribution is
                    % exponential. With that the sum over spike counts can
                    % be done as an integral worked out analytically. (See
                    % appendix in Gammor et al.
                    lambda = lambda(2:end,:); % Only timebins 2 and later can be used.
                    % Variable naming follows the notation in Gammor et al.
                    % Need to add eps to avoid infs in the logs.
                    lfTilde = lambda.*o.fTilde ; % Used a few times
                    k1  = sqrt(0.5*pi).*exp(0.5*(lambda*o.measurementNoise).^2).*erfc((lambda*o.measurementNoise-o.fTilde/o.measurementNoise)/sqrt(2));
                    z   = exp(-0.5*(o.fTilde/o.measurementNoise).^2+lfTilde);

                    ll = -(sum(log(lambda+eps)- lfTilde+ log(k1+eps),[1 2],"omitnan"));
                    if o.hasDerivatives > 0 && nargout >1
                        gamma =  (1-lfTilde-o.measurementNoise.*z.*lambda./k1 +(lambda*o.measurementNoise).^2);
                        derivative = -tensorprod(dLogLambda(:,2:end,:),gamma,[2 3],[1 2]);
                    end
                    if o.hasDerivatives>1 && nargout>2
                        % Hessian
                        A = -lfTilde+ 2*o.measurementNoise*lambda.^2-o.measurementNoise*z.*lambda./k1.*(lfTilde+o.measurementNoise*lambda.*(z./k1 -o.measurementNoise*lambda)+1);
                        A  = repmat(reshape(A,[1 o.nrTimePoints-1 o.nrTrials]),[o.nrParms 1 1]);
                        B= (1-lfTilde- o.measurementNoise*z.*lambda./k1 + (o.measurementNoise*lambda).^2);
                        part1 = tensorprod(A.*dLogLambda(:,2:end,:),dLogLambda(:,2:end,:),[2 3],[2 3]);
                        part2 = tensorprod(d2LogLambda(:,:,2:end,:),B,[3 4],[1 2]);
                        hessian  = (part1 +  part2);
                    end
                    
                case 'EXPONENTIALLOGARITHMIC'
                        n = o.nrTimePoints*o.nrTrials;
                        beta = parms(1); 
                        p = exp(parms(2));
                        
                        ll = -n.*log(-log(p)+eps) + n.*log(beta) + n.*log(1-p+eps) - beta.*sum(o.fluorescence,"all") - sum(log(1-(1-p).*exp(-beta.*o.fluorescence)),"all");
                        if o.hasDerivatives > 0 && nargout >1

                            %derivative = 
                        end
            end
        end

    end

    %% Static
    methods (Static)
        % Pre-defined tuning functions for parametric fits. If you want to add
        % additional tuning functions, note how these fit the log of the rate
        % and logs of the parameters to enforce positive rates without putting any
        % constraints on the underlying parameters that the optimization
        % routine varies. That way fminunc can be used.

        function [y, firstDerivative, secondDerivative] = logVonMises(x,parms,domain)
            % Log Von Mises tuning function y = k*cos((x - mu)) + b;
            % INPUT
            % x - Stimulus angles in degrees. Use NaN for blanks (no stimulus)
            % parms - parameter vector - [offset, preferred , kappa]
            % domain - Domain of angles ([360] or 180)
            %
            % OUTPUT
            %  y - Value per x
            %  firstDerivative  - first derivative (analytic)
            %  secondDerivative - second derivative (analytic)
            %
            % SEE demos/parametric
            arguments
                x (:,:) double
                parms (1,3) double
                domain (1,1) double = 180
            end
            offset = parms(1); preferred = parms(2); kappa = exp(parms(3));
            domainMultiplier = 360/domain*pi/180; % If only considering orientation, map angles to [0 180]

            delta = (x-preferred);
            y = kappa*cos(domainMultiplier*delta) + offset;

            % The analytic derivatives (checked against the numeric finite
            % differences).

            % [dy/doffset  dy/dpreferred dy/dkappa ]
            firstDerivative = [ones(1,length(x));
                domainMultiplier*kappa*sin(domainMultiplier*delta);
                kappa*cos(domainMultiplier*delta)]; %Yes kappa, because it is exp(parms(3))

            secondDerivative = zeros(3,3,length(x));
            %d/dpreferred,dpreferred
            secondDerivative(2,2,:) = -domainMultiplier^2*kappa*cos(domainMultiplier*delta);
            %d/dpreferred,dkappa
            secondDerivative(2,3,:) = domainMultiplier*kappa*sin(domainMultiplier*delta);
            % d/dkappa,dpreferred (same)
            secondDerivative(3,2,:) = secondDerivative(2,3,:);
            % d/dkappa,dkappa
            secondDerivative(3,3,:) = kappa*cos(domainMultiplier*delta);

            % Handle blank (no stimulus)
            blankInds = isnan(x);
            y(blankInds) = offset;
            firstDerivative(:,blankInds) = 0;
            secondDerivative(:,:,blankInds) = 0;
        end

        function [y,firstDerivative] = logTwoVonMises(x,parms)
            % Use this to fit a direction selective tuning function with
            % one bump at the preferred, and a (potentially smaller) bump
            % at the anti-preferred (=preferred +180).
            %
            % This uses the log of a the sum of two Von Mises functions:
            % y = amp1*exp(kappa*cos((x - preferred))) +amp2*exp(kappa*cos((x - preferred-180)))+offset;
            %
            % INPUT
            % x - Stimulus angles in degrees. Use NaN for blanks (no stimulus)
            % parms - parameter vector - [offset, preferred , kappa, amp1 ,amp2]
            %
            % OUTPUT
            %  y - Value per x
            % derivative - Analytic derivatice of log(y) with respect to the
            %           five parameters. This helps the optimization
            %           procedure converge.
            %
            % SEE demos/twoVonMises for an example usage.

            arguments 
                x (:,:) double
                parms (1,5) double
            end
%             arguments (Output)
%                 y (1,:) double
%                 firstDerivative (5,:) double
%             end
            offset = exp(parms(1)); preferred = parms(2); kappa = exp(parms(3));amp1=exp(parms(4)); amp2=exp(parms(5));
            deg2rad =pi/180;
            term1 = amp1*exp(kappa*cos(deg2rad*(x-preferred)));
            term2 = amp2*exp(kappa*cos(deg2rad*(x-preferred-180)));
            twoVonM = term1 + term2 + offset;
            y = log(twoVonM);
            % The analytic derivatives (checked against the numeric finite
            % differences).
            % [dy/doffset  dy/dpreferred dy/dkappa dy/damp1 dy/damp2]
            firstDerivative = repmat(1./twoVonM,[5 1]).*cat(1,offset*ones(1,numel(x)),...
                term1.*kappa.*deg2rad.*sin(deg2rad*(x-preferred)) + term2.*deg2rad.*kappa.*sin(deg2rad*(x-preferred-180)),...
                kappa.*cos(deg2rad.*(x-preferred)).*term1 +kappa.*cos(deg2rad.*(x-preferred-180)).*term2,...
                amp1.*exp(kappa*cos(deg2rad.*(x-preferred))),...
                amp2.*exp(kappa*cos(deg2rad.*(x-preferred-180))));

            % Handle blank (no stimulus)
            blankInds = isnan(x);
            y(blankInds) = offset;
            firstDerivative(:,blankInds) = 0;

        end

    end

end

