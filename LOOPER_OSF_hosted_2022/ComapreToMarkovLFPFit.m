%% Bayes opt

SAVE_FILE = 'BayesoptResultsLFP_NoDM.mat';
% SAVE_FILE = 'BayesoptResultsLFP.mat';

numVariables = 2;

lowerBounds = zeros(numVariables,1);
upperBounds = ones(numVariables,1);

MAX_ITERATIONS = 100;

ranges = {};
ranges.pcCountMin = 4;
ranges.pcCountMax = 40;
ranges.epsilonMin = 10;
ranges.epsilonMax = 50;
ranges.tMin = 5;
ranges.tMax = 30;
ranges.kMin = 3;
ranges.kMax = 20;
ranges.stateCountsMin = 30;
ranges.stateCountsMax = 90;

pcCount = optimizableVariable('pcCount',[ranges.pcCountMin,ranges.pcCountMax],'Type','integer');
epsilon = optimizableVariable('epsilon',[ranges.epsilonMin,ranges.epsilonMax],'Type','integer');
t = optimizableVariable('t',[ranges.tMin,ranges.tMax],'Type','integer');
k = optimizableVariable('k',[ranges.kMin,ranges.kMax],'Type','integer');
stateCounts = optimizableVariable('stateCounts',[ranges.stateCountsMin,ranges.stateCountsMax],'Type','integer');

fun = @(x)ComapreToMarkovLFPFunctionNoDM(x.pcCount, x.stateCounts);
% fun = @(x)ComapreToMarkovLFPFunction(x.pcCount, x.epsilon, x.t, x.k, x.stateCounts);
% ComapreToMarkovLFPFunctionNoDM(10, 10);
% ComapreToMarkovLFPFunction(10, 10, 10, 10, 10);

if exist(SAVE_FILE,'file')
    load(SAVE_FILE);
    results = resume(BayesoptResults,...
    'SaveFileName', SAVE_FILE, 'MaxObjectiveEvaluations', MAX_ITERATIONS, ...
     'OutputFcn',{@saveToFile}, 'NumCoupledConstraints',1);
else
%     results = bayesopt(fun, [pcCount,stateCounts],'Verbose',1,'AcquisitionFunctionName',...
    results = bayesopt(fun, [pcCount,epsilon,t,k,stateCounts],'Verbose',1,'AcquisitionFunctionName',...
    'expected-improvement-plus', 'MaxObjectiveEvaluations', MAX_ITERATIONS, ...
    'SaveFileName', SAVE_FILE, ...
    'OutputFcn',{@saveToFile}, 'NumCoupledConstraints',1);
end

% results = bayesopt(fun,[pcCount,epsilon,t,k,stateCounts],'Verbose',1,...
%     'UseParallel', false, 'MaxObjectiveEvaluations', 200,...
%     'AcquisitionFunctionName','expected-improvement-plus', 'NumCoupledConstraints',1)

% options = optimoptions('particleswarm','SwarmSize',numVariables*10, 'Display', 'iter', 'PlotFcn', 'pswplotbestf', 'UseParallel', 1, 'MaxStallIterations', numVariables*1000, 'MaxIterations', numVariables*1000);
% bestParams = particleswarm(fun, numVariables, lowerBounds, upperBounds, options);

% options = optimoptions(@patternsearch,'Display','iter', 'UseCompletePoll', true, 'PlotFcn', 'psplotbestf', 'MeshTolerance', 1e-20, 'StepTolerance', 1e-20);
% [bestParams,Fps, ~, output] = patternsearch(fun,ones(5,1)*0.5,[],[],[],[],lowerBounds,upperBounds,options)

% options = optimoptions(@patternsearch,'Display','iter', 'UseCompletePoll', true, 'UseParallel', true, 'PlotFcn', 'psplotbestf', 'MeshTolerance', 1e-20, 'StepTolerance', 1e-20);
% [bestParams,Fps, ~, output] = patternsearch(lossFunction,[],[],[],[],[],lowerBounds,upperBounds,options)


%% Pattern search

numVariables = 5;

lowerBounds = zeros(numVariables,1);
upperBounds = ones(numVariables,1);

ranges = {};
ranges.pcCountMin = 20;
ranges.pcCountMax = 50;
ranges.epsilonMin = 10;
ranges.epsilonMax = 40;
ranges.tMin = 10;
ranges.tMax = 30;
ranges.kMin = 3;
ranges.kMax = 10;
ranges.stateCountsMin = 80;
ranges.stateCountsMax = 120;

fun = @(x)ComapreToMarkovFunction(x(1), x(2), x(3), x(4), x(5), ranges);

% % options = optimoptions('particleswarm','SwarmSize',numVariables*10, 'Display', 'iter', 'PlotFcn', 'pswplotbestf', 'UseParallel', 1, 'MaxStallIterations', numVariables*1000, 'MaxIterations', numVariables*1000);
% % bestParams = particleswarm(fun, numVariables, lowerBounds, upperBounds, options);

options = optimoptions(@patternsearch,'Display','iter', 'UseCompletePoll', true, 'PlotFcn', 'psplotbestf', 'MeshTolerance', 1e-20, 'StepTolerance', 1e-20);
[bestParams,Fps, ~, output] = patternsearch(fun,ones(5,1)*0.5,[],[],[],[],lowerBounds,upperBounds,options)

% options = optimoptions(@patternsearch,'Display','iter', 'UseCompletePoll', true, 'UseParallel', true, 'PlotFcn', 'psplotbestf', 'MeshTolerance', 1e-20, 'StepTolerance', 1e-20);
% [bestParams,Fps, ~, output] = patternsearch(lossFunction,[],[],[],[],[],lowerBounds,upperBounds,options)

%% Pattern search LOOPER

% load('socialMovieSingle.mat');
% 
% numVariables = 3;
% 
% lowerBounds = zeros(numVariables,1);
% upperBounds = ones(numVariables,1);
% 
% ranges = {};
% ranges.nnMin = 2;
% ranges.nnMax = 20;
% ranges.repopulationDensityMin = 0.01;
% ranges.repopulationDensityMax = 0.99;
% ranges.stateCountsMin = 5;
% ranges.stateCountsMax = 400;
% 
% fun = @(x)ComapreToMarkovFMRIFunctionLOOPER(eventData, preloadedData, [], [], saveData, x(1), x(2), x(3), ranges);
% 
% options = optimoptions(@patternsearch,'Display','iter', 'UseCompletePoll', true, 'PlotFcn', 'psplotbestf', 'MeshTolerance', 1e-20, 'StepTolerance', 1e-20);
% [bestParams,Fps, ~, output] = patternsearch(fun,[0.5 0.5 0.5],[],[],[],[],lowerBounds,upperBounds,options)



