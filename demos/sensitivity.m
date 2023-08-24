%% Simulate data
% A simple test of the poissyFit routines.
%
% Simulate a single neuron with a PO of 90 degree and a fixed von Mises
% kappa and a variable peak firing rate. 
% Estimate the tuning curve for each, show the results, and compare the
% estimated with the ground truth firing rate.

% rng default

% Define "experiment" and simulated neural responses
% Varying these parameters will give some insight into the needed number of
% trials, binsizes, etc.
oriPerTrial      = 0:30:330;  % These are the orientations in the experiment
nrRepeats        =  12;     % Each orientation is shown this many times
nrTimePoints     = 10;      % 10 "bins" in each trial
dt               = 0.1;     % One bin is 100 ms.
tau              = 1.5;     % Fluorescence indicator decay
fPerSpike        = 25;      % Ca/Fluorescence per spike
measurementNoise = 10;  % Stdev of the noise
kappa = .1;
nrBootstraps =100;


close all;
spikeRates = [1 5 10 20 40];
amplitude = log(spikeRates*dt/exp(exp(kappa)));
nrAmplitudes =numel(amplitude);
est = nan(nrAmplitudes,1);
err = nan(nrAmplitudes,1);
for i=1:numel(amplitude)

    tuningParms      = [amplitude(i) 90 kappa]; % Offset Preferred log(Kappa)
    % Use the built-in von Mises function, periodic over 180 (i.e. orientation not direction).
    tc               = @(x,parms)  poissyFit.logVonMises(x,parms,180);

    % Generate data
    ori         = repmat(oriPerTrial,[1 nrRepeats]);
    nrTrials    = numel(oriPerTrial)*nrRepeats;
    lambda      =  exp(tc(ori,tuningParms)); % Lambda, the poisson rate, per trial
    LAMBDA     = repmat(lambda,[nrTimePoints 1]); % Same lambda each time point
    nrSpikes    = exprnd(LAMBDA); % The spike counts
    if tau>0
        decayFun    = fPerSpike*exp(-(0:100)*dt/tau); % F/Ca decay
        pad         = zeros(numel(decayFun)-1,1); % Start with blank
        fluorescence = nan(size(nrSpikes));
        % Simulate that trials are widely spaced, so no Ca spillover from trial to
        % trial
        for tr=1:nrTrials
            fluorescence(:,tr) = conv([pad;nrSpikes(:,tr)],decayFun','valid'); % Generate the ca signal
        end
    else
        fluorescence = nrSpikes;
    end
    fluorescence = fluorescence + normrnd(0,measurementNoise,size(fluorescence)); % Add noise

    %% Now use the poissyFit object to estimate the tuning curve parameters.
    % Setup the object by passing the orientation shown in each trial, the
    % corresponding fluorescence, the time bin and the tuning function we wish
    % to estimate.
    %
    o = poissyFit(ori(1,:),fluorescence,dt,tc,hasDerivatives = 2,tau=tau,fPerSpike=fPerSpike);
    % Make sure the object's assumptions match those of the experiment
    o.measurementNoise = measurementNoise;
    % Set options of the optimization (fminunc)
    o.options =    optimoptions(@fminunc,'Algorithm','trust-region', ...
        'MaxIter',1e8, ...
        'SpecifyObjectiveGradient',true, ...
        'HessianFcn','objective',...  % The logVonMises returns the Hessian; use it
        'FiniteDifferenceType','central', ...
        'display','off');
    % Solve
    solve(o,nrBootstraps)
    % Collect results
    err(i,:) = o.parmsError(1);
    est(i,:) = o.parms(1);
    % Show results for this amplitude
    figure;
    yyaxis right
    plot(ori,lambda/o.binWidth,'g')
    hold on
    plot(o)
    title(sprintf('Spike Rate: %.0f',spikeRates(i)));
end

%% Compare estimated rates with simulated rates across amplitudes
rateEst = exp(est)/(dt/exp(exp(kappa)));
upper  = exp(est+norminv(0.99)*err)/(dt/exp(exp(kappa)));
lower  = exp(est-norminv(0.99)*err)/(dt/exp(exp(kappa)));
figure;
plot(spikeRates,rateEst)
hold
errorbar(spikeRates,rateEst,rateEst-lower,upper-rateEst,'.')
xlabel 'Simulated Rate (spk/s)';
ylabel 'Estimated Rate (spk/s)';
