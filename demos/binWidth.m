%% BinWidth dependnce
load ../data/exampleData.mat

%% Initialize


binWidthFactor = [1 2 4 8];

nrBinWidths= numel(binWidthFactor);
r       = nan(nrRois,nrBinWidths);
rSpk    = nan(nrRois,nrBinWidths);
rCross  = nan(nrRois,nrBinWidths);
parms  = nan(nrRois,nrBinWidths,5);
parmsError=nan(nrRois,nrBinWidths,5);

%% FIT
% For each ROI, fit a logTwoVonMises, bootstrap the parameter estimates and
% determined the splitHalves correlation.
nrBoot = 100;
nrWorkers = gcp('nocreate').NumWorkers ; % Parfor for bootstrapping
spikeCountDist = "POISSON";
for j = 1:nrBinWidths
    thisTimes =seconds(0.5:stepSize*binWidthFactor(j):1.5);
    thisF = retime(f,thisTimes,'linear');
    [nrTimePoints, nrTrials] = size(thisF);
    thisF = permute(double(reshape(thisF.Variables,[nrTimePoints nrRois nrTrials])),[1 3 2]);
    thisNP = retime(f,thisTimes,'linear');
    thisNP = permute(double(reshape(thisNP.Variables,[nrTimePoints nrRois nrTrials])),[1 3 2]);
    thisF = thisF - 0.7*thisNP;
    thisSpk = retime(spk,thisTimes,'linear');
    thisSpk= permute(double(reshape(thisSpk.Variables,[nrTimePoints nrRois nrTrials])),[1 3 2]);
    for roi =1:nrRois 
        fprintf('BinWidth #%d ROI #%d (%s)\n',j,roi,datetime('now'))

        o = poissyFit(direction,thisF(:,:,roi),stepSize*binWidthFactor(j),@poissyFit.logTwoVonMises);
        o.spikeCountDistribution = spikeCountDist;        
        
        switch spikeCountDist
            case "POISSON"
                o.hasDerivatives = 1;
                o.options =    optimoptions(@fminunc,'Algorithm','trust-region', ...
                    'SpecifyObjectiveGradient',true, ...
                    'HessianFcn','objective', ...
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
        o.measurementNoise =stdF(roi)/sqrt(binWidthFactor(j));
        o.nrWorkers = nrWorkers;
        try
            solve(o,nrBoot);
            % Store the results
            parms(roi,j,:) =o.parms;
            parmsError(roi,j,:)= o.parmsError;
            [r(roi,j),~,rSpk(roi,j),~,rCross(roi,j)] = splitHalves(o,nrBoot,[],thisSpk(:,:,roi));
        catch me
            fprintf(2,'Failed on %d\n',roi)
        end
    end
end


%% Show Results
figure(1);
clf
noiseGroups = repmat(1:nrBinWidths,[nrRois 1]);
subplot(2,2,1)
plot(binWidthFactor*stepSize, mean(r,1,"omitnan"));
hold on
plot(binWidthFactor*stepSize, mean(rSpk,1,"omitnan"));
xlabel 'BinWidth (s)'
ylabel '<\sigma>'
title (sprintf('Tuning curve reliability  '))
legend('F','Spk')
subplot(2,2,2)
plot(binWidthFactor*stepSize,mean(parmsError(:,:,2),1,"omitnan"))
xlabel 'BinWidth (s)'
ylabel 'stdev (deg)'


subplot(2,2,3)
scatter(r(:,1),rSpk(:,1));
axis square;
hold on
plot(xlim,ylim,'k')
xlabel '\sigma_F'
ylabel '\sigma_{spk}'
title(sprintf('Tuning curve reliability %.2f (p=%.3g)',mean(r(:,1)-rSpk(:,1)),ranksum(r(:,1),rSpk(:,1))))



