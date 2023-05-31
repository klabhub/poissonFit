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
%%
%{
% Code to generate the exampleData file - needs access to the klabData and
DataJoint server

roi = sbx.Roi & 'session_date = ''2023-03-23''' & 'pcell>0.95';
expt = ns.Experiment & 'session_date = ''2023-03-23''' & 'paradigm like ''Ori%''' & 'starttime=''12:56:11''';
stepSize = 1/15.5;
trialsToKeep =2:120; % The first trial has some nan in it at the start.
nrRois = 1000;
fetchOptions = ['ORDER BY pCell Desc LIMIT ' num2str(nrRois)] ; % Take the top nrRois neurons
[f] = get(roi ,expt, trial= trialsToKeep, modality = 'fluorescence',start=0.5,stop=1.5,step=stepSize,interpolation ='nearest',fetchOptions = fetchOptions);
[np] = get(roi ,expt,trial= trialsToKeep,  modality = 'neuropil',start=0.5,stop=1.5,step=stepSize,interpolation ='nearest',fetchOptions = fetchOptions);
stdF = [fetch(roi,'stdfluorescence',fetchOptions).stdfluorescence];
[spk] = get(roi,expt,trial= trialsToKeep, modality = 'spikes',start=0.5,stop=1.5,step=stepSize,interpolation ='nearest',fetchOptions = fetchOptions);
direction = get(expt ,'gbr','prm','orientation','atTrialTime',0);
direction =direction(trialsToKeep);

save exampleData f np stdF spk direction nrRois stepSize

%}
load ../data/exampleData

thisTimes =seconds(0.5:stepSize:1.5);
f = retime(f,thisTimes,'linear');
[nrTimePoints, nrTrials] = size(f);
f = permute(double(reshape(f.Variables,[nrTimePoints nrRois nrTrials])),[1 3 2]);
np = retime(np,thisTimes,'linear');
np = permute(double(reshape(np.Variables,[nrTimePoints nrRois nrTrials])),[1 3 2]);
F  = f-0.7*np;
spk = retime(spk,thisTimes,'linear');
spk= permute(double(reshape(spk.Variables,[nrTimePoints nrRois nrTrials])),[1 3 2]);


%% Initialize
nrRois  = size(F,3);
nrRois = 10;
r       = nan(nrRois,1);
rSpk    = nan(nrRois,1);
rCross  = nan(nrRois,1);
parms  = nan(nrRois,5);
parmsError=nan(nrRois,5);

%% FIT
% For each ROI, fit a logTwoVonMises, bootstrap the parameter estimates and
% determined the splitHalves correlation.
nrBoot = 100;
nrWorkers = gcp('nocreate').NumWorkers ; % Parfor for bootstrapping
spikeCountDist = "POISSON";
for roi =1:nrRois
    fprintf('ROI #%d (%s)\n',roi,datetime('now'))
    o = poissyFit(direction,F(:,:,roi),stepSize,@poissyFit.logTwoVonMises);
    o.spikeCountDistribution = spikeCountDist;
    switch spikeCountDist
        case "POISSON"
            o.hasDerivatives = 1;
            o.options =    optimoptions(@fminunc,'Algorithm','trust-region', ...
                'SpecifyObjectiveGradient',true, ...
                'display','none', ...
                'CheckGradients',false, ... % Set to true to check supplied gradients against finite differences
                'diagnostics','off');
        case "EXPONENTIAL"
            o.hasDerivatives = 0;
            o.options =    optimoptions(@fminunc,'Algorithm','quasi-newton', ...
                'SpecifyObjectiveGradient',false, ...
                'display','none', ...
                'diagnostics','off');
    end
    o.measurementNoise =stdF(roi);
    o.nrWorkers = nrWorkers;
    try
        solve(o,nrBoot);
        % Store
        parms(roi,:) =o.parms;
        parmsError(roi,:)= o.parmsError;
        [r(roi),~,rSpk(roi),~,rCross(roi)] = splitHalves(o,nrBoot,[],spk(:,:,roi));
    catch me
    end
end
%% Show Results
figure(1);
clf
subplot(2,2,1)
scatter(r,rSpk);
axis equal
ylim([-1 1])
xlim([-1 1])

hold on
plot(xlim,xlim,'k')
xlabel '\sigma_F'
ylabel '\sigma_{spk}'
title (sprintf('Tuning curve reliability %.2f (p=%.3g)',mean(r-rSpk),ranksum(r,rSpk)))

subplot(2,2,2)
scatter(r,parmsError(:,2),'.')
xlabel '\sigma_F'
ylabel 'stdev (deg)'
hold on
plot(xlim,[20 20])

subplot(2,2,3)
scatter(rSpk,parmsError(:,2),'.')
xlabel '\sigma_{spk}'
ylabel 'stdev (deg)'
hold on
plot(xlim,[20 20])



%% Save
save ("../data/directionTuning" + spikeCountDist + ".mat", 'r','rSpk', 'rCross','parms','parmsError')