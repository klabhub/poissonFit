%% Simulate Data

% Define "experiment" and simulated neural responses
oriPerTrial      = 0:30:330;  % These are the orientations in the experiment
nrRepeats        =  12;    % Each orientation is shown this many times
nrTimePoints     = 10;      % 10 "bins" in each trial
dt               = 0.1;     % One bin is 100 ms.
tau              = .5;      % Fluorescence indicator decay
fPerSpike        = 50;      % Ca/Fluorescence per spike
measurementNoise = 5;  % Stdev of the noise
tuningParms      = [0 90 0]; % Offset Preferred log(Kappa)
% Use the built-in von Mises function, periodic over 180 (i.e. orientation not direction).
tc               = @(x,parms)  poissyFit.logVonMises(x,parms,180);
groundTruthTuningCurve = exp(tc(oriPerTrial,tuningParms));
% Generate data
ori         = repmat(oriPerTrial,[1 nrRepeats]);
nrTrials    = numel(oriPerTrial)*nrRepeats;
lambda      =  exp(tc(ori,tuningParms)); % Lambda, the poisson rate, per trial
LAMBDA     = repmat(lambda,[nrTimePoints 1]); % Same lambda each time point
nrSpikes    = poissrnd(LAMBDA); % The spike counts
decayFun    = fPerSpike*exp(-(0:100)*dt/tau); % F/Ca decay
pad         = zeros(numel(decayFun)-1,1); % Start with blank
fluorescence = nan(size(nrSpikes));
% Simulate that trials are widely spaced, so no Ca spillover from trial to
% trial
for tr=1:nrTrials
    fluorescence(:,tr) = conv([pad;nrSpikes(:,tr)],decayFun','valid'); % Generate the ca signal
end
fluorescence = fluorescence + normrnd(0,measurementNoise,size(fluorescence)); % Add noise