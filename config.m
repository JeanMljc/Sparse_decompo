%% ==========================
%  config.m
%  ==========================
%  Project configuration script
%  Sets up paths, solver
%  Configure before executing the main code.
%  ==========================

%% --- Paths ---

% load path to subfolders
folderPath = '';
% folderPath = '/home/jean/Documents/CODE_matlab/CODE_transmission/';

% load path to CPLEX 12.10 
CPLEXPath = '';
% CPLEXPath = '/home/jean/CPLEX_1210/cplex/matlab/x86-64_linux';

% load code from [Elvander20]
ElvanderPath = '';
% ElvanderPath = '/home/jean/Documents/Papers on subject/Code/multimarginal_omt/multimarginal_omt/tracking';

% Add paths
addpath(genpath(folderPath));
addpath(genpath(CPLEXPath));
addpath(genpath(ElvanderPath));