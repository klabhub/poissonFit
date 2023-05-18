%% Parametric Demo 
% This demonstrates how to fit a von Mises with two bumps 

%% Simulate data
% A single neuron responds to a set of gratings with an
% orientation tuned (von Mises) response.

rng default
% Define "experiment" and simulated neural responses
oriPerTrial      = 0:30:330;  % These are the orientations in the experiment
nrRepeats        =  12;    % Each orientation is shown this many times
nrTimePoints     = 10;      % 10 "bins" in each trial
dt               = 0.1;     % One bin is 100 ms.
tau              = 1.2;      % Fluorescence indicator decay
fPerSpike        = 50;      % Fluorescence per spike
measurementNoise = 1;  % Stdev of the noise
tuningParms      = [0 90 .01 .5 .1]; % log(Offset)Preferred log(Kappa) log(amp1) log(amp2)
% Use the built-in logTwoVonMises function
tc               = @poissyFit.logTwoVonMises;

% Generate data
ori         = repmat(oriPerTrial,[1 nrRepeats]);
nrTrials    = numel(oriPerTrial)*nrRepeats;
lambda      =  exp(tc(ori,tuningParms)); % Lambda, the poisson rate, per trial
LAMBDA      = repmat(lambda,[nrTimePoints 1]); % Same lambda each time point
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

%% Now use the poissonFit object to estimate the tuning curve parameters.
% Setup the object by passing the orientation shown in each trial, the
% corresponding fluorescence, the time bin and the tuning function we wish
% to estimate.
% 
o = poissyFit(ori(1,:),fluorescence,dt,tc);


% Make sure the object's assumptions match those of the experiment
o.tau =tau;
o.fPerSpike = fPerSpike;
o.measurementNoise = measurementNoise;
% Set options of the optimization (fminunc)
o.options =    optimoptions(@fminunc,'Algorithm','trust-region', ...
    'SpecifyObjectiveGradient',true, ...
    'HessianFcn',[],... % Approximate using finite differences
    'display','final-detailed', ...
    'diagnostics','off', ...   
    'MaxFunctionEvaluations',5e3);

% Boostrap solutions
% No initial value (guess) is provided; poissyFit will estimate the
% preferred from the rates (non-parametrically) and start from that point.
% That seems to be key to finding a good solution.
solve(o,100);
figure;
yyaxis right
plot(ori,lambda/o.binWidth,'g')

plot(o); % Show the result.
hold on

