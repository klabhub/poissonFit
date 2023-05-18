%% Parametric Demo 
% This demonstrates how to fit a specific tuning curve (here a von Mises)
% to a calcium imaging fluorescnce trace.

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
measurementNoise = 0.1;  % Stdev of the noise
tuningParms      = [0 90 0]; % Offset Preferred log(Kappa)
% Use the built-in von Mises function, periodic over 180 (i.e. orientation not direction).
tc               = @(x,parms)  poissonFit.logVonMises(x,parms,180);

% Generate data
ori         = repmat(oriPerTrial,[1 nrRepeats]);
nrTrials    = numel(oriPerTrial)*nrRepeats;
lambda      =  exp(tc(ori,tuningParms)); % Lambda, the poisson rate, per trial
lambda      = repmat(lambda,[nrTimePoints 1]); % Same lambda each time point
nrSpikes    = poissrnd(lambda); % The spike counts
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
o = poissonFit(ori(1,:),fluorescence,dt,tc);
% Make sure the object's assumptions match those of the experiment
o.tau =tau;
o.fPerSpike = fPerSpike;
% Set options of the optimization (fminunc)
o.options =    optimoptions(@fminunc,'Algorithm','trust-region', ...
    'MaxIter',1e8, ...
    'SpecifyObjectiveGradient',true, ...
    'display','final-detailed', ...
    'CheckGradients',false, ... 
    'diagnostics','on', ...
    'FiniteDifferenceType','central'); % Central is slower but more accurate


% Solve 
solve(o);
figure;
plot(o); % Show the result.
hold on
plot(ori,lambda,'g')

%% Bootstrap
% To asses the robustness of the fit,the tool uses bootstrapping: trials
% are sampled with replacement, parameters are estimated for each
% bootstrap set, and the final answer is the average over the bootstrap
% sets (o.parms) and the error is the standard deviation over bootstrap
% stes (o.parmsError). Compainr o.parmsError to the o,parms will reveal
% whether the fit is reliable.
nrBoostrapSets = 1000;
o.options =    optimoptions(@fminunc,'Algorithm','trust-region', ...
    'MaxIter',1e8, ...
    'SpecifyObjectiveGradient',true, ...    
    'display','off', ...
    'CheckGradients',false, ... 
    'diagnostics','off'); 

solve(o,nrBoostrapSets); 
figure;
plot(o)
hold on
plot(ori,lambda,'g')


