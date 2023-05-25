%% Measurement Noise dependnce

load ../data/exampleData

%% Initialize
nrRois  = size(F,3);

noiseFactors = [0.1 0.5 1 2 5 10];
nrNoiseFactors= numel(noiseFactors);
r       = nan(nrRois,nrNoiseFactors);
rSpk    = nan(nrRois,nrNoiseFactors);
rCross  = nan(nrRois,nrNoiseFactors);
parms  = nan(nrRois,nrNoiseFactors,5);
parmsError=nan(nrRois,nrNoiseFactors,5);

%% FIT
% For each ROI, fit a logTwoVonMises, bootstrap the parameter estimates and
% determined the splitHalves correlation.
nrBoot = 100;
nrWorkers = 8 ; % Parfor for bootstrapping
spikeCountDist = "POISSON";
for roi =1:nrRois
    fprintf('ROI #%d (%s)\n',roi,datetime('now'))
    o = poissyFit(direction,F(:,:,roi),binWidth,@poissyFit.logTwoVonMises);
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
    for j = 1:numel(noiseFactors)
        o.measurementNoise =noiseFactors(j)*stdF(roi);
    o.nrWorkers = nrWorkers;
    try
    solve(o,nrBoot);
    
    parms(roi,j,:) =o.parms;
    parmsError(roi,j,:)= o.parmsError;
    [r(roi,j),~,rSpk(roi,j),~,rCross(roi,j)] = splitHalves(o,nrBoot,[],spk(:,:,roi));
    catch me
    end
    end
end


%% Show Results
figure(1);
clf
noiseGroups = repmat(1:nrNoiseFactors,[nrRois 1]);
subplot(2,2,1)
plot(noiseFactors, mean(r,1,"omitnan"));
hold on
plot(noiseFactors, mean(rSpk,1,"omitnan"));
xlabel 'Noise factor'
ylabel '<\sigma>'
title (sprintf('Tuning curve reliability  '))
legend('F','Spk')
subplot(2,2,2)
plot(noiseFactors,mean(parmsError(:,:,2),1,"omitnan"))
xlabel 'Noise factor'
ylabel 'stdev (deg)'

