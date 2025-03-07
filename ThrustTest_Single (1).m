clear; clc; close all;

%% ThrustTest_Single Summary
% First, some of the summary from Thrust:
% This funciton will take in the file location of the two test setups and
% using the file names in those directories, will pull out all of the
% avialable tests, process them, and fit the data into a standard
% formatting for output
%
% ThrustTest_Single is meant to be a testing function where you can test
% the conditioning steps that your group develops. This version of the code
% loads in only one data set at a time so that the workspace is less
% cluttered. It is suggested that students make plots of their thrust data
% often throughout their development in this code section to visually check
% that their conditioning is working as desired. Once conditioning is
% working on one set of data, students are engcouraged to try other single
% data sets. Once groups are satisfied, the full funciton has the same form
% as this, so the conditioning code can simply be dragged and dropped into
% the full "Thrust.m" funciton.
% 
% Note that there is one major step missing from this testing funciton.
% That is averaging over multiple tests from the same setup. This will need
% to be added in the full "Thrust.m" function


%% Outputs:
% ThrustCurves:
%   A table containing 0.5 seconds of thrust data for each of the cases
%   available, this data will have formatting such that there are 501
%   evenly spaced thrust data points (rows) for each test (columns). The
%   ordering of the columns will go from max to min water volume in the 2L
%   bottle and then max to min in the 1.25L bottle
%
% Time:
%   A 1D array corresponding to the times of the thrust data points in the
%   ThrustCurves table
%
% <User defined variables for statistics>

% Set known sampling frequency
f= 1652; % [Hz]

%% Preallocate variables of interest
Time = 0:0.001:0.5; % just go ahead and define this, note that it will be 501 long
ThrustCurves = zeros(length(Time),1);

%% List what configuration is being used
bottleSize = 2000; % [ml]
waterSize = 800;% [ml]
% Be sure that the above matches up with the file name specified below. In
% the full "Thrust" script, all of this data will be loaded automatically
% for you
testName = 'Static Test Stand Data/1250mL Bottle/Test_T02_W0625_B1250'; % Should be a string of the path to the data (including the data file name)

% /////////////////////////////////////////////////////////////////////////
% MODIFY THIS SECTION
% /////////////////////////////////////////////////////////////////////////

%% Load data
% This should not have to be modified
fileName = testName; % again weird indexing is due to string arrays, we have to ask for all the characters in a row
data = readmatrix(fileName); % load the data
data = data(:,3)*4.448; % take only the third column and converting from lbf to N


%% Data Conditioning

%% Data Fitting

%% Sample onto the standard output array format

% /////////////////////////////////////////////////////////////////////////
% END OF SECTION TO MODIFY
% /////////////////////////////////////////////////////////////////////////


