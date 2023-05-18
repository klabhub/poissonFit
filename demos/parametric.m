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
fPerSpike        = 50;      % Ca/Fluorescence per spike
measurementNoise = 0;  % Stdev of the noise
tuningParms      = [0 90 0]; % Offset Preferred log(Kappa)
% Use the built-in von Mises function, periodic over 180 (i.e. orientation not direction).
tc               = @(x,parms)  poissyFit.logVonMises(x,parms,180);

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

%% Now use the poissyFit object to estimate the tuning curve parameters.
% Setup the object by passing the orientation shown in each trial, the
% corresponding fluorescence, the time bin and the tuning function we wish
% to estimate.
% 
o = poissyFit(ori(1,:),fluorescence,dt,tc);
o.hasDerivatives = 2; % The logVonMises has both derivatives and hessian output
% Make sure the object's assumptions match those of the experiment
o.tau =tau;
o.fPerSpike = fPerSpike;
% Set options of the optimization (fminunc)
o.options =    optimoptions(@fminunc,'Algorithm','trust-region', ...
    'MaxIter',1e8, ...
    'SpecifyObjectiveGradient',true, ...
    'HessianFcn','objective',...  % The logVonMises returns the Hessian; use it
    'display','final-detailed', ...
    'CheckGradients',false, ... 
    'diagnostics','on', ...
    'FiniteDifferenceType','central'); % Central is slower but more accurate


% Solve 
solve(o);
figure;
yyaxis right
plot(ori,lambda/o.binWidth,'g')
hold on
plot(o); % Show the result.

%% Bootstrap
% To asses the robustness of the fit,the tool uses bootstrapping: trials
% are sampled with replacement, parameters are estimated for each
% bootstrap set, and the final answer is the average over the bootstrap
% sets (o.parms) and the error is the standard deviation over bootstrap
% stes (o.parmsError). Compainr o.parmsError to the o,parms will reveal
% whether the fit is reliable.
nrBoostrapSets = 100;
o.options =    optimoptions(@fminunc,'Algorithm','trust-region', ...
    'SpecifyObjectiveGradient',true, ...    
    'HessianFcn','objective',... % The logVonMises returns the Hessian.
    'display','off', ...
    'diagnostics','off'); 

solve(o,nrBoostrapSets); 
figure;
yyaxis right
plot(ori,lambda/o.binWidth,'g')
plot(o)


