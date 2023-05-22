%% Direction Tuning
% An example application; using poissyFit to determine direction tuning of
% 1000 neurons.
% The data set has :
% F - Neuropil corrected Fluorescence [nrTimePoints nrTrials nrROI]
% stdF - Standard deviation of F
% spk - Deconvolved spiking activity  
% binWidth 
% direction - Direction of motion in each trial [1 nrTrials]
%
% These data were recorded with Scanbox/Neurolabware and preprocessed using
% default parameters of Suite2p.
%

load ../data/exampleData

%% Initialize
nrRois  = size(F,3);
r       = nan(nrRois,1);
rSpk    = nan(nrRois,1);
rCross  = nan(nrRois,1);
parms  = nan(nrRois,5);
parmsError=nan(nrRois,5);

%% FIT
% For each ROI, fit a logTwoVonMises, bootstrap the parameter estimates and
% determined the splitHalves correlation.
nrBoot = 100;
nrWorkers = 8 ; % Parfor for bootstrapping
for roi =299:nrRois 
    fprintf('ROI #%d (%s)\n',roi,datetime('now'))
    o = poissyFit(direction,F(:,:,roi),binWidth,@poissyFit.logTwoVonMises);
    o.spikeCountDistribution = "POISSON";
    o.hasDerivatives = 1;
    o.options =    optimoptions(@fminunc,'Algorithm','trust-region', ...
            'SpecifyObjectiveGradient',true, ...
            'display','none', ...
            'CheckGradients',false, ... % Set to true to check supplied gradients against finite differences
            'diagnostics','off');          
    o.measurementNoise =stdF(roi);
    o.nrWorkers = nrWorkers;    
    solve(o,nrBoot);
    parms(roi,:) =o.parms;
    parmsError(roi,:)= o.parmsError;
    [r(roi),~,rSpk(roi),~,rCross(roi)] = splitHalves(o,nrBoot,[],spk(:,:,roi));
  
end
%% Show Results
figure(1);
clf
scatter(r,rSpk);
axis square;
hold on
plot(xlim,ylim,'k')
xlabel '\sigma_F'
ylabel '\sigma_{spk}'
title 'Tuning curve reliability'
