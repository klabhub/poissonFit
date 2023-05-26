%% Parameter exploration
%
%  Fits tuning curves to experimental data using a range of parameter
%  choices and determines a goodness of fit (based on split halves and circular
% standard deviation). This allows a comparison of the spike and Fluorescence fits 
% across a range of parameter choices (to assess the fitting procedures
% sentisitivity).
%

%% Load a data set with 1000 V1 rois
load ../data/exampleData.mat

%% Initialize

nrBoot = 10;  % Nr Bootstrap sets used to estimate circular standard deviation as well as split halves.
nrWorkers = gcp('nocreate').NumWorkers; % Parfor for bootstrapping
spikeCountDist = "POISSON"; % POISSON or EXPONENTIAL
nrRoisToRun = 3;    % Of the 1000 roi in the file, run the first n.
parm = 'noise';     % Which parameter should be varied (others get the default value)
defaultTau = 1.3;   % Default calcium decay tau (Gcamp6s)    
defaultBin = 1/15.5; % Default bin size (framerate of acquisition)
defaultNoiseFactor = 1; % Default noise factor (the suite2p estimate of the noise is multiplied by this).
fun = @poissyFit.logTwoVonMises; 

% Determine the range of the varied parameter
switch (parm)
    case 'tau'
        prms = [0.1 0.5 1 2  5 10]*defaultTau;
        prmLabel = 'tau (s)';
    case 'binwidth'
        prms =  [1 2 4 8]*defaultBin;
        prmLabel = 'bin width (s)';
    case 'noise'
        prms = [0.1 0.5 1 2 5 10];
        prmLabel =  'noise factor ';
    otherwise
        error('Unknown parm')
end

% Preallocate
nrPrms= numel(prms);
r       = nan(nrRois,nrPrms);
rSpk    = nan(nrRois,nrPrms);
rCross  = nan(nrRois,nrPrms);
parms  = nan(nrRois,nrPrms,5);
parmsError=nan(nrRois,nrPrms,5);

%% FIT
% For each ROI, fit a logTwoVonMises, bootstrap the parameter estimates and
% determined the splitHalves correlation.

for j = 1:nrPrms
    thisTau =defaultTau;
    thisBin = defaultBin;
    thisNoise = defaultNoiseFactor;
    switch (parm)
    case 'tau'
        thisTau = prms(j);
    case 'binwidth'
        thisBin =prms(j);
    case 'noise'
        thisNoiseFactor =prms(j);
    otherwise
        error('Unknown parm')
    end

    % Extrac the neuropil corrected fluorescence and deconvolved spikes at
    % the requested binwidth
    thisTimes =seconds(0.5:thisBin:1.5);
    thisF = retime(f,thisTimes,'linear');
    [nrTimePoints, nrTrials] = size(thisF);
    thisF = permute(double(reshape(thisF.Variables,[nrTimePoints nrRois nrTrials])),[1 3 2]);
    thisNP = retime(f,thisTimes,'linear');
    thisNP = permute(double(reshape(thisNP.Variables,[nrTimePoints nrRois nrTrials])),[1 3 2]);
    thisF = thisF - 0.7*thisNP;
    thisSpk = retime(spk,thisTimes,'linear');
    thisSpk= permute(double(reshape(thisSpk.Variables,[nrTimePoints nrRois nrTrials])),[1 3 2]);
    for roi =1:nrRoisToRun
        fprintf('Parm #%d ROI #%d (%s)\n',j,roi,datetime('now'))        
        o = poissyFit(direction,thisF(:,:,roi),thisBin,fun);
        o.spikeCountDistribution = spikeCountDist;
        o.tau = thisTau;
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
        o.measurementNoise =stdF(roi)*thisNoiseFactor;
        o.nrWorkers = nrWorkers;
        try
            solve(o,nrBoot); % Bootstrap the results
            % Store 
            parms(roi,j,:) =o.parms;
            parmsError(roi,j,:)= o.parmsError;
            % Bootstrap the splitHalves correlation
            [r(roi,j),~,rSpk(roi,j),~,rCross(roi,j)] = splitHalves(o,nrBoot,[],thisSpk(:,:,roi));
        catch me
            fprintf(2,'Failed on %d\n',roi)
        end
    end
end


%% Show Results
figure(1);
clf

subplot(2,2,1)
plot(prms, mean(r,1,"omitnan"));
hold on
plot(prms, mean(rSpk,1,"omitnan"));
xlabel(prmLabel)
ylabel '<\sigma>'
title (sprintf('Tuning curve reliability  '))
legend('F','Spk')

subplot(2,2,2)
plot(prms,mean(parmsError(:,:,2),1,"omitnan"))
xlabel(prmLabel)
ylabel 'stdev (deg)'


subplot(2,2,3)
scatter(r(:,1),rSpk(:,1));
axis square;
hold on
plot(xlim,ylim,'k')
xlabel '\sigma_F'
ylabel '\sigma_{spk}'
title(sprintf('Tuning curve reliability %.2f (p=%.3g)',mean(r(:,1)-rSpk(:,1)),ranksum(r(:,1),rSpk(:,1))))


%% Save the results

