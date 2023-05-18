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
tau              = .5;      % Fluorescence indicator decay
fPerSpike        = 50;      % Fluorescence per spike
measurementNoise = 1;  % Stdev of the noise
tuningParms      = [0 90 .01 1 .5]; % log(Offset)Preferred log(Kappa) log(amp1) log(amp2)
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
fluorescence = conv([pad;nrSpikes(:)],decayFun','valid'); % Generate the ca signal
fluorescence  =reshape(fluorescence,[nrTimePoints nrTrials]);
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
% Set options of the optimization (fminunc)
o.options =    optimoptions(@fminunc,'Algorithm','trust-region', ...
    'MaxIter',1e8, ...
    'SpecifyObjectiveGradient',true, ...
    'display','final-detailed', ...
    'CheckGradients',true, ... 
    'diagnostics','off', ...   
    'FiniteDifferenceType','central',...
    'MaxFunctionEvaluations',5e3); % Central is slower but more accurate

% Solve 
guess = [0 0 0 0 0];
solve(o,1,guess);
figure;
plot(o); % Show the result.
hold on
yyaxis left
scatter(ori,lambda./max(lambda)*0.95*max(ylim),'k*')

