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
            'NumThreads',1,...   % Speeds up remote processes                        
            'AdditionalPaths',additionalPaths, ...
            'StartupFolder','/home/bart/Documents/MATLAB',...
            'CurrentFolder','/home/bart/poissyFit/demos');
c.gitFolders = gitFoldersOnAmarel;
c.gpo;  % Update remote repo

%%  Run a script
% Us a pool for the parfor inside the poissyFit.
% With 32 on Amarel, the analysis ran at 1 ROI per 2-3 seconds. 
expression = "directionTuning";
job = script(c,expression,'Pool',11); 


%% Retrieve the results 
wait(job)
load(job, 'r','rSpk', 'rCross','parms','parmsError','gof','bf','bfR','bfParms','bfError')

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
xlabel 'r_F'
ylabel 'r_{spk}'
title (sprintf('Split halves correlation: Delta = %.2f (p=%.3g)',mean(rSpk-r),ranksum(r,rSpk)))

subplot(2,2,2)
scatter(r,parmsError(:,2),'.')
xlabel 'r'
ylabel 'stdev (deg)'
hold on
plot(xlim,[20 20])

hold on
scatter(rSpk,parmsError(:,2),'.')
xlabel 'r'
ylabel 'stdev (deg)'
hold on
plot(xlim,[20 20])
title 'Bootstrap StdDev'

legend('dF/F','spk')

subplot(2,2,3)
notNan = ~isnan(bfR);
scatter(rSpk(notNan),bfR(notNan))
hold on
axis equal
ylim([-1 1])
xlim([-1 1])
plot(xlim,xlim,'k')
ylabel 'r_{bf}'
xlabel 'r_{spk}'
title (sprintf('Split halves correlation: Delta = %.2f (p=%.3g)',mean(bfR(notNan)-rSpk(notNan)),ranksum(bfR(notNan),rSpk(notNan))))

subplot(2,2,4)
scatter(bfR(notNan),bfError(3,notNan),'.')
xlabel 'r_{bf}'
ylabel 'stdev (deg)'
hold on
plot(xlim,[20 20])
title 'Bootstrap StdDev'




