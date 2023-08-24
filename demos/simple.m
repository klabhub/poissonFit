%% Simple Demo 
% This demonstrates how to determine the Poisson rates from fluorescence
% data with a trial structure.

simulateData
%% Now use the poissonFit object to estimate the rate in each condition.
% Setup the object by passing the orientation shown in each trial, the
% corresponding fluorescence, and the time bin. 
% % Compared with the parametric demo the only difference is that we do not
% specify a tuning function here. The object will then fit one rate per
% unique stimulus value (here a direction) plus a constant for a total of
% 13 parameters. 
nrBoostrapSets =100;
o = poissyFit(ori(1,:),fluorescence,dt,tau=tau,fPerSpike=fPerSpike); 
% Make sure the object's assumptions match those of the experiment

o.measurementNoise = measurementNoise;
% Set options of the optimization (fminunc)
o.options =    optimoptions(@fminunc,'Algorithm','trust-region', ...
    'SpecifyObjectiveGradient',true, ...
    'HessianFcn','objective',...
    'display','final-detailed', ...
    'CheckGradients',false, ... 
    'diagnostics','on', ...
    'FiniteDifferenceType','central'); % Central is slower but more accurate

% Solve 
solve(o);
figure;
yyaxis right
plot(oriPerTrial,groundTruthTuningCurve/o.binWidth,'g')
plot(o); % Show the result.

%% Compare F with Spikes fit
[r,~,rSpk,~,~] = splitHalves(o,nrBoostrapSets,[],nrSpikes);
fprintf("PoissyFit: %.3f vs. Spikes %.3f (split-halves correlation; higher is better, 1 is noise-free maximum)\n",r,rSpk)
