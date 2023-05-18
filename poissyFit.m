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
    % additinional functionality to fit double von Mises. 
    %
    % See README.md , demos/parammetric, demos/simple, demos/twoVonMises
    %
    % BK - May 2023.
    properties (SetAccess=public,GetAccess=public)
        % Meta parameters
        maxRate         = 100 %  The maximum number of spikes (per second) that a neuron can reasonable fire under the conditions of the experiment

        % Model parameters
        stimHistory     = 0.1 % in seconds. was stimHistoryLength
        tau             = 1.2 % Fluorescence decay time constant [s]
        fPerSpike      = 80  % Fluorescence per spike. This is a true free parameter
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
        pFGivenSpikes           % [nrSpikesMax nrObservations] - Probability of the observed F given the spikes
        nrSpikesFactorial       % nrSpikes!
        nrSpikesPerTimepoint    % [nrSpikesMax nrObservations] - Matrix with nrSpikes in each column
        DM                      % [nrPredictors nrObservations] - Design matrix
        % The raw data
        fluorescence            % [nrTimePoints nrTrials]
        stimulus                % [1 nrTrials]
        binWidth                % Time bins for the fluorescence signal in seconds.
     
        uStimulus               % Unique stimulus values
        stimulusIx              % Index into uStimulus for each trial
        tuningCurve             % Non-parametric estimate of the tuning curve
        tuningCurveErr          % IQR for tuning curve
        bestGuess               % Initial guess of the parameters (from nonparametric tuning)

        parms                   % Estimated parameters that capture the fluorescence best
        parmsError              % Errors on the parameter estimates.
        bootstrap               % Number of bootstrap sets used to estimate parms/parmsError.
        % If this is 1, the parms are the result of
        % the ML fit, if it is >1 they are the
        % mean of the bootstrap sets (and the error
        % the standard deviation across sets).
        exitFlag                % Exitflag of the optimization routine.
    end

    properties (Dependent)
        nrTimePoints        % Time points per trial
        nrTrials            % Number of trials
        nrObservations      % Total number of Fluorescence observations (nrTimePoints*nrTrials)
        nrSpikesMax         % Maximum number of spikes in a single bin (derived from maxRate)
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


    end

    %% Public functions
    methods (Access=public)

        function o = poissyFit(stim,fluor,dt,fun)
            % Construct an object, based on the stimulus,fluorescence and
            % timebins of the data.
            arguments
                stim (1,:) double   % The stimulus presented in each trial
                fluor (:,:) double  % The response in each bin (rows) and each trial (rows)
                dt (1,1) double     % The duration of a single bin (in seconds).
                fun = []            % Tuning function for parametric fits.
            end
            o.stimulus = stim;
            o.fluorescence = fluor;
            o.binWidth = dt;
            o.tuningFunction = fun;
            % Check whether the fun has deriviatives. This fails if the fun
            % is an anonymous function - the user then has to set
            % .hasDeriviatives manually.
            if isempty(fun)
                o.hasDerivatives  = 2;
            elseif nargout(fun) <0
                warning('Please set o.hasDerivatives manually for this anonymous tuning function')
            else
                o.hasDerivatives  = nargout(fun)-1;
            end
        end

        function v = copyWithNewData(o,stim,fluor,dt,fun)
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
            end
            v = copy(o); % Make a shallow copy
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
                lambda = exp(tensorprod(parms,o.DM,2,1));
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
                lambda=reshape(repmat(lambda,[o.nrTimePoints 1]),[1 o.nrTimePoints o.nrTrials]);

            end
        end

        function plot(o)
            % Plot an overview of the data and tuning curve fits/
            if isempty(o.tuningFunction)

                predicted = reshape(lambdaFromStimulus(o,o.parms),[o.nrTimePoints o.nrTrials]);
                meanPredicted = mean(predicted,1,"omitnan"); % Average over time points

                predictedTuningCurve= accumarray(o.stimulusIx,meanPredicted',[],@median);

                predictedHigh = reshape(lambdaFromStimulus(o,o.parms+o.parmsError),[o.nrTimePoints o.nrTrials]);
                meanPredictedHigh = mean(predictedHigh,1,"omitnan"); % Average over time points
                predictedTuningCurveHigh= accumarray(o.stimulusIx,meanPredictedHigh',[],@median);
                predictedLow = reshape(lambdaFromStimulus(o,o.parms-o.parmsError),[o.nrTimePoints o.nrTrials]);
                meanPredictedLow = mean(predictedLow,1,"omitnan"); % Average over time points
                predictedTuningCurveLow= accumarray(o.stimulusIx,meanPredictedLow',[],@median);
            else
                predictedTuningCurve = exp(o.tuningFunction(o.uStimulus,o.parms));
                predictedTuningCurveHigh = exp(o.tuningFunction(o.uStimulus,o.parms+o.parmsError));
                predictedTuningCurveLow= exp(o.tuningFunction(o.uStimulus,o.parms-o.parmsError));
            end
            yyaxis left
            h1= errorbar(o.uStimulus,o.tuningCurve,o.tuningCurveErr/2,'b');
            ylabel 'Fluorescence (a.u.)'
            yyaxis right
            hold on
            
                pos = predictedTuningCurveHigh- predictedTuningCurve;
                neg = predictedTuningCurve - predictedTuningCurveLow;

                h2= errorbar(o.uStimulus,1./o.binWidth*predictedTuningCurve,1./o.binWidth*pos,1./o.binWidth.*neg,'r');
            ylabel 'Lambda (spk/s)'
            xlabel 'Stimulus'
            title (['Parms: ' num2str(o.parms)])
           legend([h1(1) h2(1)],'Fluorescence with IQR','Rate with SE')

        end



        function [bootParms,bootParmsError] = solve(o,boot,guess)
            %[bootParms,bootParmsError] = solve(o,boot,guess)
            % Solve the fitting problem. Optionally do this in bootstrap
            % resampling fashion
            % boot - Number of boostrap sets
            % guess - Optional initial guess of the parameters of the fit.
            %        If left unspecified, it uses a non-parametric estimate
            %        of the tuning cruve.
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
            % Solve the minimization problem.
            [o.parms,o.exitFlag,~,~,~,hessian] = fminunc(@o.logLikelihood,guess,o.options);
            iH = diag(inv(hessian));
            iH(iH<0) = NaN; % Negative element means we're not at a (local) minimum. Probably a flat piece of the LL.
            o.parmsError = sqrt(iH)';

            % Perform bootstrap resampling to get a robust estimate of the
            % parms and their errors
            if boot>1
                h = waitbar(0,'Bootstrapping');
                bootParms       = nan(boot,numel(o.parms));
                bootParmsError  = nan(boot,numel(o.parms));
                for i=2:boot
                    h = waitbar(i/boot,h,'Bootstrapping');
                    thisTrials = randi(o.nrTrials,[1 o.nrTrials]);
                    thisO = o.copyWithNewData(o.stimulus(:,thisTrials),o.fluorescence(:,thisTrials),o.binWidth,o.tuningFunction);
                    solve(thisO,1,guess);
                    bootParms(i,:) = thisO.parms;
                    bootParmsError(i,:) = thisO.parmsError;

                end
                close(h)
                % Store the mean and std across sets as the outcome
                o.parms      = mean(bootParms,1,'omitnan');
                o.parmsError = std(bootParms,0,1,'omitnan');
                o.bootstrap  = boot;
            else
                o.bootstrap = 1;
            end
        end


    end



    %% Protected, internal functions
    methods  (Access=protected)
        function initialize(o)
            o.nrSpikes             = (0:o.nrSpikesMax)';
            o.nrSpikesFactorial    = repmat(factorial(o.nrSpikes),[1 o.nrObservations]);
            o.nrSpikesPerTimepoint = repmat(o.nrSpikes,[1 o.nrObservations]);
            % The F is determined  by the previous time bin's F (which decays at a rate of tau)
            % plus the new F generated by the spikes in the current bin
            % F is  [nrTimePoints nrTrials] so the first bin cannot really
            % be used as its preceding bin may not be in the matrix (e.g.
            % if only a subset of trials is analyzed). 
            F = permute(repmat(o.fluorescence,[1 1 o.nrSpikesMax+1]),[3 1 2]);
            fExpected = o.fPerSpike.*o.nrSpikesPerTimepoint(:,2:end,:) + F(:,1:end-1,:).*exp(-o.binWidth/o.tau);
            % The probability of observing an F given the spikes is a
            % normal distribution with a stdev equal ot the measuremnt
            % noise.Adding eps to avoid underflow.
            
            o.pFGivenSpikes = eps+normpdf(F(:,2:end,:) - fExpected,0,o.measurementNoise);
            % The first bin has no preceding time bin; set it to NaN.
            o.pFGivenSpikes = cat(2, nan(o.nrSpikesMax+1, 1,o.nrTrials),o.pFGivenSpikes);
            % Construct the design matrix
            o.DM = designMatrix(o);

            % Determine a tuning curve (median rate per unique stimulus
            % value)
            [o.uStimulus,~,o.stimulusIx] = unique(o.stimulus);
            mResponse = mean(o.fluorescence,1,"omitnan"); % Average over time bins
            o.tuningCurve = accumarray(o.stimulusIx,mResponse',[],@median);
            o.tuningCurveErr = accumarray(o.stimulusIx,mResponse',[],@(x) (diff(prctile(x,[25 75]),1)));

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
                o.bestGuess = [];
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


        function [ll, derivative, hessian] = logLikelihood(o,parms)
            % Calculates the log likelihood given the data and parameters
            % returns negative of log likelihood, plus the analytic
            % derivatives, Hessian (if available and requested)
            switch o.hasDerivatives
                case 2
                    [lambda,dLogLambda,d2LogLambda] = lambdaFromStimulus(o,parms);
                    dLogLambda = dLogLambda(:,2:end,:);
                    d2LogLambda= d2LogLambda(:,:,2:end,:);
                case 1
                    [lambda,dLogLambda] = lambdaFromStimulus(o,parms);
                    dLogLambda = dLogLambda(:,2:end,:);
                otherwise
                    lambda = lambdaFromStimulus(o,parms);
            end
            lambda = lambda(:,2:end,:);
            LAMBDA= repmat(lambda,[o.nrSpikesMax+1 1 1]);

            pX = eps+ exp(-LAMBDA).*LAMBDA.^o.nrSpikesPerTimepoint(:,2:end,:)./o.nrSpikesFactorial(:,2:end,:);
            pXY = pX.*o.pFGivenSpikes(:,2:end,:);
            sum_PXY = sum(pXY,1); % Sum over nrSpikes
            ll = -sum(log(sum_PXY),[2 3]); % Sum over nrTimePoints and nrTrials

            % Derivatives, if requested and available.
            if o.hasDerivatives>0 && nargout>1
                sumX_PXY = tensorprod(o.nrSpikes,pXY,1,1);
                sumX_PXY_over_sum_PXY = sumX_PXY./sum_PXY;
                gamma = sumX_PXY_over_sum_PXY - lambda;
                derivative = -tensorprod(dLogLambda,gamma,[2 3],[2 3]);
            end
            if o.hasDerivatives>1 && nargout>2
                % Hessian
                c = lambda - tensorprod(o.nrSpikes.^2,pXY./sum_PXY,1,1) + (sumX_PXY_over_sum_PXY).^2;
                part1  = tensorprod((c.*dLogLambda),dLogLambda,[2 3],[2 3]);
                part2 = - sum(d2LogLambda.*reshape(gamma,1,1,o.nrTimePoints-1,o.nrTrials),[3,4]);
                hessian = part1 + part2;
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

            arguments (Input)
                x (:,:) double
                parms (1,5) double
            end
            arguments (Output)
                y (1,:) double
                firstDerivative (5,:) double
            end
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

