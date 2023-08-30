% Compare poissyFit with bayesFit
% A sketch only. Due to the problems with the BF in bayesFit, I gave up on
% using bayesFit. 

nrBoot = 10;        % Bootstrap iterations for splitHalves


% Load the data and rearrange. (see directionTuning demo for details)
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

% Pick an ROI
roi =100;

nrWorkers = gcp('nocreate').NumWorkers ; % Parfor for bootstrapping
    
fprintf('ROI #%d (%s)\n',roi,datetime('now'))
o = poissyFit(direction,F(:,:,roi),stepSize,@poissyFit.logTwoVonMises,fPerSpike=500,scaleToMax=true);
o.hasDerivatives = 1;
o.options =    optimoptions(@fminunc,'Algorithm','trust-region', ...
                    'SpecifyObjectiveGradient',true, ...
                    'display','none', ...
                    'CheckGradients',false, ... % Set to true to check supplied gradients against finite differences
                    'diagnostics','off');
o.measurementNoise =stdF(roi);
o.nrWorkers = nrWorkers;
solve(o,nrBoot);


oSpk = o.copyWithNewData(o.stimulus,spk(:,:,roi),o.binWidth,o.tuningFunction,tau=0,fPerSpike=1, scaleToMax=true);
% Scale the measurement noise by the nonparametric estimate
% of the noise in the tuning curves
oSpk.measurementNoise = mean(oSpk.tuningCurveErr)./mean(o.tuningCurveErr).*o.measurementNoise;
solve(oSpk);

[r,~,rSpk] = splitHalves(o,nrBoot,[],spk(:,:,roi));

figure(1);
clf
subplot(2,2,1) 
plot(o)

subplot(2,2,2)
plot(oSpk)
%% Compare with bayesFit

nrParms = 4; % Circular gaussian 360 has 4 parms
x= repmat(direction,[nrTimePoints 1]);
x=x(:);
fprintf('BayesFit ROI #%d (%s)\n',roi,datetime('now'))
y=spk(:,:,roi);
y = y(:);

subplot(2,2,3)
bayesResults = bayesFit(x,y,"fun","direction_selective_circular_gaussian","graphics",true);
    

subplot(2,2,4)
bayesResults = bayesFit(direction,mean(spk(:,:,roi),1),"fun","direction_selective_circular_gaussian","graphics",true);



 [rBf,bf] = splitHalves(direction',mean(spk(:,:,roi),1)',"fun","direction_selective_circular_gaussian","nrBoot",nrBoot,"nrWorkers",nrWorkers);

