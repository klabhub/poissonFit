%% Parametric Demo 
% This demonstrates how to fit a specific tuning curve (here a von Mises)
% to a calcium imaging fluorescnce trace.

%% Simulate data
% A single neuron responds to a set of gratings with an
% orientation tuned (von Mises) response.
simulateData;

%% Now use the poissyFit object to estimate the tuning curve parameters.
% Setup the object by passing the orientation shown in each trial, the
% corresponding fluorescence, the time bin and the tuning function we wish
% to estimate.
% 
o = poissyFit(ori(1,:),fluorescence,dt,tc,tau =tau,fPerSpike=fPerSpike);
o.hasDerivatives = 2; % The logVonMises has both derivatives and hessian output
% Make sure the object's assumptions match those of the experiment
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
plot(oriPerTrial,groundTruthTuningCurve/o.binWidth,'g')
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
plot(oriPerTrial,groundTruthTuningCurve/o.binWidth,'g')
plot(o)

%% Compare F with Spikes fit
[r,~,rSpk,~,~] = splitHalves(o,nrBoostrapSets,[],nrSpikes)

