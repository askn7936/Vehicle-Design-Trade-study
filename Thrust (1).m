function [ThrustCurves, Time] = Thrust()
clear; clc; close all;
%% Thrust Summary
% This funciton will take in the file location of the two test setups and
% using the file names in those directories, will pull out all of the
% avialable tests, cendition their data, and fit that data into a standard
% formatting for output. Note that statistics are also requested for the
% student deliverable, but how to pass those out will be left up to the
% students as they are not needed to be passed into any later functions.
% Despite this, the first two outputs of the funciton are not permitted to
% have their form modified.

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
% <User defined variable(s) for statistics>
%

%% Define data locations
% This is hard coded!!!
fileLoc_2L = '/Users/ECABRERA22/Documents/MATLAB/2000ml/'; % path to the data files, be sure to include a trailing slash
fileLoc_1pt25L = '/Users/ECABRERA22/Documents/MATLAB/1250ml/'; % path to the data files, be sure to include a trailing slash

%% Read in all of the avilable data and find what data there is
testInfo_2L = getThrustTestNames(fileLoc_2L);
configs_2L = unique(testInfo_2L.waterVol);
numConfigs_2L = length(configs_2L);

testInfo_1pt25L = getThrustTestNames(fileLoc_1pt25L);
configs_1pt25L = unique(testInfo_1pt25L.waterVol);
numConfigs_1pt25L = length(configs_1pt25L);

numConfigs = numConfigs_2L + numConfigs_1pt25L;

% Set known sampling frequency
f= 1652; % [Hz]

%% Preallocate variables of interest
Time = 0:0.001:0.5; % just go ahead and define this, note that it will be 501 long
ThrustCurves = zeros(length(Time),numConfigs);

ThrustCurvesNames = {};
m=0;
%% Loop over all of the configurations
for N = 1:numConfigs % use upper case N to distiguish that it is counting something different from the aerodynamic modeling loops
    %% Dertemine what configuration to use for this iteration in the loop
    if N <=  numConfigs_2L % determine if we should be reading 2L or 1.25L data
        bottleSize = '2000'; % [ml]
        waterSize = configs_2L(N);
        testIndexes = find(testInfo_2L.waterVol == waterSize); % finds the index of the relavant tests
        numTests = length(testIndexes); % finds the number of tests performed
        testNames = testInfo_2L.fileNames(testIndexes, :); % pulls all of the test names of interest, weird indexing is due to string arrays
    else
        bottleSize = '1250'; % [ml]
        waterSize = configs_1pt25L(N-numConfigs_2L);
        testIndexes = find(testInfo_1pt25L.waterVol == waterSize); % finds the index of the relavant tests
        numTests = length(testIndexes); % finds the number of tests performed
        testNames = testInfo_1pt25L.fileNames(testIndexes, :); % pulls all of the test names of interest, weird indexing is due to string arrays
    end

    % /////////////////////////////////////////////////////////////////////////
    % MODIFY THIS SECTION
    % /////////////////////////////////////////////////////////////////////////
    % Notice that there is little to no guidance in place for this
    % function. This is on purpose as there are many different and equally
    % valid ways to process data (not to say that any way is valid though).
    % The lack of guidance is therefore to encourage you to think about,
    % discuss, and potentially debate as a group the best set of steps to
    % extract just the meaningful part of the thrust profile

    datasum=0;
    for v = 1:numTests
    data=0;
    %% Load data
        % The folloowing three lines will pull all of the files in each
        % test setup for you and give the array "data" which should be
        % conditioned. You should not need to modify any of this section of
        % code
        fileName = testNames(v, :); % again weird indexing is due to string arrays, we have to ask for all the characters in a row
        data = readmatrix(fileName); % load the data
        data = data(:,3)*4.448; % take only the third column and converting from lbf to N

    %% Data Conditioning
[cutoff_value, cutoff_Index]=max(data);

%loop reset
data1=0;
data2full=0;

%declare bounds for 2 sections w/ respect to the max. Max position is lined
%up 50 intervals after start to calc average.
j=0;
for i=(cutoff_Index-50):cutoff_Index
    j=j+1;
    data1(j)=data(i);
    if data1(j)<0
        data1(j)=0;
    end
end
j=0;
for i=cutoff_Index:(length(data))
    j=j+1;
    data2full(j)=data(i);
    if data2full(j)<0
        data2full(j)=0;
    end
end

%combine both arrays for plotting
plot_data=zeros(1,length(Time));

for i=1:length(data1)
    plot_data(i)=data1(i);
end
for i=1:(length(Time)-length(data1))
    plot_data(length(data1)+i)=data2full(i);
end

%Use for plotting all individual thrust data
% m=m+1;
% figure(m)
% plot(Time,plot_data)
% title(fileName)

%% Averaging

%each loop add data sets together
datasum=plot_data+datasum;

    end

    %average calculation
datacut=datasum./numTests;

    % Note that averaging should accour before data fitting. Technically
    % either can be done, but the output will be much more smooth if the
    % fit is applied to an average

%% Data Fitting

[cutoff_valueavg,cutoff_avg]=max(datacut);

%iteration reset, avoids keeping old values
x1=0;
x2=0;
x3=0;
datafit1=0;
datafit2=0;
datafit3=0;
dataval1=0;
dataval2=0;
dataval3=0;

%bound value to end at 0.5 seconds
data2=zeros(1,length(Time)-length(data1)+1);
for i=1:length(data2)
    data2(i)=data2full(i);
end

%create lines of best fit
per=1/f;
x1=per.*(1:length(data1));
x2=per.*(cutoff_avg:length(datacut));

Time1=0;
Time2=0;

[datafit1,S]=polyfit(x1,data1,1);
Time1=Time(1,1:cutoff_avg);
dataval1=polyval(datafit1,Time1, S);

[datafit2,S]=polyfit(x2,data2,2);
Time2=Time(1,1:end-cutoff_avg);
dataval2=polyval(datafit2,Time2,S);

%Create 3rd line of best fit (flat)
%Find when data2 reaches zero
%if data doesnt reach zero
endpoint=Time(dataval2==min(dataval2));

%if data does reach zero
if Time(dataval2<=0)
endpoint=Time(dataval2<=0);
end

%Redeclare data2 and Time2 to end at endpoint
Time3end=length(endpoint(1):0.001:Time(end)-Time1(end));

data2end=0;
Time2end=0;
for i=1:length(Time2)-Time3end+1
    data2end(i)=dataval2(i);
    Time2end(i)=Time2(i)+Time1(end);
end

%y-value of line and best fit for flat section
data3=zeros(1,Time3end)+data2end(end);

x3=per.*(length(datacut):length(datacut)+length(data3)-1);

[datafit3,S]=polyfit(x3,data3,1);
Time3=Time(:,1:length(data3));
dataval3=polyval(datafit3,Time3, S);

%line-up section 1 slope and remove neg values
% if max(dataval1)>cutoff_value
%     dataval1=dataval1-(max(dataval1)-cutoff_value);
% elseif max(dataval1)<cutoff_value
%     dataval1=dataval1+(cutoff_value-max(dataval1));
% end
% for i=1:length(dataval1)
%     if dataval1(i)<0
%         dataval1(i)=0;
%     end
% end

%plot graphs for each volume
figure(N);
plot(Time1,dataval1);
hold on
plot(Time,plot_data);
plot(Time2end,data2end);
plot(Time3+Time2end(end),dataval3);
title(waterSize)
hold off

plot_data(plot_data<0)=0;

thrustOut=plot_data;
% /////////////////////////////////////////////////////////////////////////
% END OF SECTION TO MODIFY
% /////////////////////////////////////////////////////////////////////////
    %% Convert to table for output
    % It is very important that the data is 501 elements long corresponding
    % to 0-0.5 seconds of time at this point!!!
    ThrustCurves(:, N) = thrustOut;
    % Header naming convention of <bottle size (in ml)>_<water volume (in ml)>
    ThrustCurvesNames{N} = [bottleSize, '_', num2str(waterSize)];
end
ThrustCurves = array2table(ThrustCurves);
ThrustCurves.Properties.VariableNames = ThrustCurvesNames; 
end
