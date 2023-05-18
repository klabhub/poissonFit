%% Simple Demo 
% This demonstrates how to determine the Poisson rates from fluorescence
% data with a trial structure.

%% Simulate data
% A single neuron responds to a set of gratings with an
% orientation tuned (von Mises) response.

rng default
% Define "experiment" and simulated neural responses
oriPerTrial      = 0:30:330;  % These are the orientations in the experiment
nrRepeats        =  12;    % Each orientation is shown this many times
nrTimePoints     = 10;  % 10 "bins" in each trial
dt               = 0.25;      % One bin is 100 ms.
tau              = .5;          % Fluorescence indicator decay
fPerSpike        = 50;      % Fluorescence per spike
measurementNoise = 0.1;  % Stdev of the noise
tuningParms      = [0 90 0]; % Offset Preferred log(Kappa)
% Use the built-in von Mises function, to gerneate a direction specific response peaking at 90
tc               = @(x,parms)  poissyFit.logVonMises(x,parms,360);

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

%% Now use the poissonFit object to estimate the rate in each condition.
% Setup the object by passing the orientation shown in each trial, the
% corresponding fluorescence, and the time bin. 
% % Compared with the parametric demo the only difference is that we do not
% specify a tuning function here. The object will then fit one rate per
% unique stimulus value (here a direction) plus a constant for a total of
% 13 parameters. 

o = poissyFit(ori(1,:),fluorescence,dt); 
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

