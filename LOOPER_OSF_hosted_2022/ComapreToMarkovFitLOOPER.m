%% Bayes opt

SAVE_FILE = 'BayesoptResultsLOOPER.mat';

numVariables = 3;

nn = optimizableVariable('nn',[4, 10],'Type','integer');
repop = optimizableVariable('repop',[0, 1],'Type','integer');
stateCounts = optimizableVariable('stateCounts',[50, 300],'Type','integer');

fun = @(x)ComapreToMarkovFunctionLOOPER(x.nn, x.repop, x.stateCounts);

if exist(SAVE_FILE,'file')
    load(SAVE_FILE);
    results = resume(BayesoptResults,...
    'SaveFileName', SAVE_FILE, 'MaxObjectiveEvaluations', 100, ...
     'OutputFcn',{@saveToFile}, 'NumCoupledConstraints',1);
else
    results = bayesopt(fun, [nn,repop,stateCounts],'Verbose',1,'AcquisitionFunctionName',...
    'expected-improvement-plus', 'MaxObjectiveEvaluations', 100, ...
    'SaveFileName', SAVE_FILE, ...
    'OutputFcn',{@saveToFile}, 'NumCoupledConstraints',1);
end

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



