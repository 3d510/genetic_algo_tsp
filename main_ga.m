%   TSP_GA Traveling Salesman Problem (TSP) Genetic Algorithm (GA).
%   Finds a (near) optimal solution to the TSP by setting up a GA to search
%   for the shortest route (least distance for the salesman to travel to
%   each city exactly once and return to the starting city)
%
% Summary:
%     1. A single salesman travels to each of the cities and completes the
%        route by returning to the city he started from
%     2. Each city is visited by the salesman exactly once
%
% Input:
%     USERCONFIG (structure) with zero or more of the following fields:
%     - xy (float) is an Nx2 matrix of city locations, where N is the number of cities
%     - popSize (scalar integer) is the size of the population (should be divisible by 4)
%     - numGen (scalar integer) is the number of desired iterations for the algorithm to run
%     - mutProb (float) is the probability of mutation per individual
%     - crossProb (float) is the probability of cross-over
%     - eliteFract (float) is used for Elitism. It is the fraction of good members of the current
%     generation which should replace the weak members in the next generation. 

% Input Notes:
%     1. Rather than passing in a structure containing these fields, any/all of
%        these inputs can be passed in as parameter/value pairs in any order instead.
%     2. Field/parameter names are case insensitive but must match exactly otherwise.
%
% Output:
%     RESULTSTRUCT (structure) with the following fields:
%         (in addition to a record of the algorithm configuration)
%     - OPTROUTE (integer array) is the best route found by the algorithm
%     - MINDIST (scalar float) is the cost of the best route

clc;
clear all;

% first, delete any old mat file from the workspace (you do not need to change this code)
if exist('cities.mat', 'file')==2
    delete('cities.mat');
end


%% Then, choose which data file you wish to experiment with.
% Below is the list of various data files with different number of cities:

% loadatt48();         % 48 cities
% loadst70();          % 70 cities
% loadgr96();          % 96 cities

% I have chosen a data file from the above list
 loadatt48();

%% prepare the distance matrix
load('cities.mat');
xy = cities';

% you should update the following code to obtain the average and 95%
% confidence interval for each configuration of numGen
for numGen = 100:100:2000
    for runs = 1:15
        userConfig = struct('xy', xy, 'popSize', 200, 'numGen', numGen, 'crossProb', 0.25, 'mutProb', 0.5, 'eliteFract', 0.02);
        resultStruct = tsp_ga(userConfig);
        
        % the best tour found by GA
        % fprintf('\nBest tour found by GA:\n');
        % resultStruct.optRoute
        
        % the distance of the best tour
        fprintf('\n Number of generations: %d \n Run number: %d \n The distance of the best tour = %d\n',numGen, runs, resultStruct.minDist);
        
    end
end

% Implement your plotting here, using the average and confidence interval results:
% plots ...

