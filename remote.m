% Use kSlurm from a desktop computer
% See the tutorial.mlx on github.com/klabhub/kSlurm how to set this up for
% your machine. 

% Matlab search path
additionalPaths = { '/home/bart/poissyFit','/home/bart/bayesFit'};

% Folders I want to git pull origin before running anything on the cluster
gitFoldersOnAmarel = {'/home/bart/poissyFit','/home/bart/bayesFit'};

% Construc the kSlurm object. Assuming that workers never need more than 
% 2 hours for the work this script sends them.
c = kSlurm('Hours',12,...  
            'NumThreads',8,...   % Speeds up remote processes                        
            'AdditionalPaths',additionalPaths, ...
            'StartupFolder','/home/bart/Documents/MATLAB',...
            'CurrentFolder','/home/bart/poissyFit/demos');
c.gitFolders = gitFoldersOnAmarel;
c.gpo;  % Update remote repo

%%  Run a script
% Us a pool for the parfor inside the poissyFit.
% With 32 on Amarel, the poissyFit analysis ran at 1 ROI per 2-3 seconds. 
expression = "directionTuning";
job = script(c,expression,'Pool',32); 
% Retrieve the results 
fprintf('Job sumitted. Waiting...')
tic;
wait(job)
toc
load(job, 'r','rSpk', 'rCross','parms','parmsError','gof','bf','rBf','parmsBf','errorBf')

%% Show Results
   figure(1);
    clf
    % Scatter of split halves correlation for df/F vs spk
    subplot(2,2,1)
    scatter(r,rSpk);
    axis equal
    ylim([-1 1])
    xlim([-1 1])
    hold on
    plot(xlim,xlim,'k')
    xlabel 'r_F'
    ylabel 'r_{spk}'
    title (sprintf('Split halves correlation: Delta = %.2f (p=%.3g)',mean(rSpk-r),ranksum(r,rSpk)))

    subplot(2,2,2)
    % Scatter of standard devaition of PO across bootstraps against the split halves correlation for df/F and spk
    scatter(r,parmsError(:,2),'.')
    xlabel 'r'
    ylabel 'bootstrap stdev PO (deg)'
    hold on
    plot(xlim,[20 20])
    scatter(rSpk,parmsError(:,2),'.')
    plot(xlim,[20 20])
    legend('dF/F','spk')


    % Bayesfit results.
      % Scatter of split halves correlation for BF vs spk
    subplot(2,2,3)
    notNan = ~isnan(rBf);
    scatter(rSpk(notNan),rBf(notNan))
    hold on
    axis equal
    ylim([-1 1])
    xlim([-1 1])
    plot(xlim,xlim,'k')
    ylabel 'r_{bf}'
    xlabel 'r_{spk}'
    title (sprintf('Split halves correlation: Delta = %.2f (p=%.3g)',mean(rBf(notNan)-rSpk(notNan)),ranksum(rBf(notNan),rSpk(notNan))))

    % Scatter of standard devaition of PO across bootstraps against the split
    % halves correlation for  BF
    subplot(2,2,4)
    scatter(rBf(notNan),errorBf(3,notNan),'.')
    xlabel 'r_{bf}'
    ylabel 'stdev (deg)'
    hold on
    plot(xlim,[20 20])





