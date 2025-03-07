%%Boost-Glide Model Baseline Main Script
% ASEN 2804
% Author: John Mah
% Modifications by: Preston Tee
% Date: Initiated - 15 Dec. 2023
% Last modified - 15 Jan. 2024

%% Clean Workspace and Housekeeping
clear
% clearvars
close all

% removes warnings for table variable names for a cleaner output
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

%% Import and Read Aircraft Design File
Design_Input = readtable("ROCHELLE_final.xlsx",'Sheet','Input','ReadRowNames',true); %Read in Aircraft Geometry File
Count = height(Design_Input); %Number of different aircraft configurations in geometry file

TailAirfoil = [0.0490,0.0239,0.0187,0.0136,0.0106,0.0077, 0.0107,0.0136,0.0159,0.0196,0.0403,0.0595,0.0810,0.1046,0.1276,0.1446,0.1631,0.1801];

% Import Airfoil Data File
Airfoil = readtable("Design Input File.xlsx",'Sheet','Airfoil_Data'); %Read in Airfoil Data

% Import Benchmark Aircraft Truth Data
Benchmark = readtable("Design Input File.xlsx",'Sheet','Benchmark_Truth'); %Read in Benchmark "Truth" Data for model validation only

% Import Material Properties Data
Material_Data = readtable("Design Input File.xlsx",'Sheet','Materials'); %Read in prototyp material densities for weight model

%% Caluations - Conditions and Sizing
% US Standard Atmophere - uses provided MATLAB File Exchange function
[rho,a,T,P,nu,z]= atmos(Design_Input.altitude_o(:,:)); 

ATMOS = table(rho,a,T,P,nu,z); % Reorganize atmopheric conditions into a table for ease of passing into functions
clearvars rho a T P nu z % Clear original variables now that they are in a table

% Call Wing Geometry Calcuation Function
WingGeo_Data = WingGeo(Design_Input,Count); %Calculate specific wing geometry from wing configuration parameters

%% Quick Explainer - Tables
% This code heavily utilizes tables for data organization. You should think
% of a table as a spreadsheet. Tables are also very similar to a 2D array
% of data exept that the columns can be named so that it is clear what data
% is in those columns. There are multiple ways to get the data out of
% columns, through standard indexing, and through dot indexing. 
%
% Standard indexing:
%
% Like when indexing into an array you can get data out of a table by using
% parenthasis, (), which will return another table, which is often a 
% problem for calulations and plotting
%   Example:
%       NewTable = OriginalTable(:,:)
%
% Alternativly, if you index in the same way but with curly braces, {}, a
% standard array will be returned
%   Example:
%       NewArray = OriginalTable{:,:}
%
% Dot indexing:
%
% Finally, if you would like to access just one column, tables support dot
% indexing using the name of the column header. This takes the form of the 
% name of the variable, then a dot, then the name of the column, This will
% return a standard 1D array of data
%   Example: 
%       NewArray = OriginalTable.ColumnName_1
%
% The tables in this code have been purposly organized such that the rows
% will ALWAYS correspond to the different configuration inputs in the input
% file. In other words, the first input row of the input file (1st row not
% including the header row) will match up with row 1 of the tables in this
% code, the second input will be the 2nd row of tables here, etc. Columns
% will always be variables of interest and will be named appropriately.
%
% More MATLAB documentation on getting data out of tables:
% https://www.mathworks.com/help/matlab/matlab_prog/access-data-in-a-table.html
%
% We have provided the necessary code for packaging the data into tables to
% output from and input to functions in order to keep the size of the
% function headers reasonable. It will be your responsibility to unpack and
% use data passed into functions in tables correctly. Please ensure you
% are using the preallocated varaible names and do not modify the code that
% creates the tables. We want to help you with the math, not with general
% coding

%% Calculations - Lift and Drag
% Call Wing Lift & Drag Model Function
[WingLiftModel,AoA,AoA_Count,AirfoilLiftCurve,WingLiftCurve,WingDragCurve] =...
   WingLiftDrag(Design_Input,Airfoil,Count); 

AoA3D_0 = interp1(WingLiftCurve{1,:}, AoA,0); %3D angle of attack at zero lift calculation
Benchmark_AoA0 = interp1(Benchmark{:,2}, Benchmark{:,1},0);
Benchmark_slope = ((Benchmark.CL(2)-Benchmark.CL(1))/(Benchmark.AoA(2)-Benchmark.AoA(1)));
CLerror = (transpose(Benchmark.CL)-AirfoilLiftCurve{1,:}).*100;

% Call Parasite Drag Buildup Model Function
[Parasite_Drag_Data,FF_Table] = ...
   ParasiteDrag(Design_Input,Airfoil,WingGeo_Data,ATMOS,Count);

% Call Induced Drag Model Function
[InducedDrag_Data] = ...
   InducedDrag(Design_Input,WingLiftModel,WingLiftCurve,WingDragCurve,WingGeo_Data,Count);

% Call Complete Drag Polar Function
[DragPolar_mod1,DragPolar_mod2,DragPolar_mod3] = ...
   DragPolar(Parasite_Drag_Data,InducedDrag_Data,Design_Input,AoA_Count,WingLiftCurve,Count);

% Call L/D Analysis Function
[LD_mod1,LD_mod2,LD_mod3,LD_benchmark] = LD(Benchmark,DragPolar_mod1,DragPolar_mod2,DragPolar_mod3,WingLiftCurve,AoA_Count,Count);

AoA0_1 = interp1(LD_mod1{1,:}, AoA,0);
AoA0_2 = interp1(LD_mod2{1,:}, AoA,0);
AoA0_3 = interp1(LD_mod3{1,:}, AoA,0);

% Call Weight Model
[Weight_Data,CG_Data] = ...
    Weight(Design_Input,Count,WingGeo_Data,Airfoil,Material_Data);

%% Calculations - Dynamic Models
%%To Be provided at a later date
% Call Thrust Model
[ThrustCurves, Time] = Thrust();

% Call Boost-Ascent Flight Dynamics Model
[apogee, hApogee, stateStruct] = BoostAscent(Design_Input, ATMOS, Parasite_Drag_Data, Weight_Data, ThrustCurves, Time, Count);

% Call Glide Flight Dynamics Model
[GlideData] = GlideDescent(LD_mod2,apogee,Design_Input, ATMOS, Weight_Data, WingLiftModel, WingLiftCurve, Count);

%% Plotting

%Lift Curves
figure
hold on
plot(AoA,Airfoil{1,(5:22)});
a_wing = polyfit(AoA, Airfoil{1,(5:22)},1);
a_tail = polyfit(AoA, TailAirfoil,1);
plot(Benchmark.AoA(:),Benchmark.CL(:));
plot(AoA,AirfoilLiftCurve{1,:});
plot(AoA,WingLiftCurve{1,:});
xlabel('Angle of Attack (deg)');
ylabel('Coefficient of Lift (CL)');
title('Lift Curve Slope Modeling Comparison');
legend('Airfoil Data','Benchmark Aircraft','Airfoil Model','Wing Model','Location','southeast');
hold off

moment_slope = a_wing(1)*((0.285-0.4024)-((a_tail(1)/a_wing(1))*0.6*(1-0.4)));
epsilon_0 = (1.62*WingLiftCurve{1,"0"})/(pi*Design_Input.AR_w(1));
CMvals = (AoA.*moment_slope)+(0.6*a_tail(1)*epsilon_0);


figure
plot(AoA, CMvals)
yline(0)
xlabel('AoA')
ylabel('Cm')
title('Pitching Moment Diagram')
trimAoA = AoA(CMvals==0);

%Drag Polar Curves
figure
hold on
plot(AirfoilLiftCurve{1,:},Airfoil{1,(24:41)});
plot(WingLiftCurve{1,:},WingDragCurve{1,:});
plot(Benchmark.CL(:),Benchmark.CD(:));
plot(WingLiftCurve{1,:},DragPolar_mod1{1,:});
plot(WingLiftCurve{1,:},DragPolar_mod2{1,:});
plot(WingLiftCurve{1,:},DragPolar_mod3{1,:});
xlabel('Coefficient of Lift (CL)');
ylabel('Coefficient of Drag (CD)');
title('Drag Polar Model Comparison');
legend('Airfoil Drag Polar','Wing Drag Polar','Benchmark Drag Polar','Drag Polar Cavallo','Drag Polar Obert','Drage Polar DATCOM','Location','northwest');
hold off

figure
hold on
for n = 1:Count
    Cdi_1 = WingLiftCurve{1,:}.^2.*InducedDrag_Data.k1_mod1(n);
    Cdi_2 = WingLiftCurve{1,:}.^2.*InducedDrag_Data.k1_mod2(n);
    Cdi_3 = WingLiftCurve{1,:}.^2.*InducedDrag_Data.k1_mod3(n);
    plot(WingLiftCurve{1,:},Cdi_1); % brace indexing for plotting tables
    plot(WingLiftCurve{1,:},Cdi_2);
    plot(WingLiftCurve{1,:},Cdi_3);
end
xlabel('Coefficient of Lift (CL)');
ylabel('Induced Drag (CDi)');
title('Induced Drag vs. CL');
legend('Cavallo','Obert','DATCOM');
hold off

figure
hold on
for n = 1:Count
    plot(WingLiftCurve{1,:},DragPolar_mod1{n,:}); % brace indexing for plotting tables
end
xlabel('Coefficient of Lift (CL)');
ylabel('Coefficient of Drag (CD)');
title('Drag Polar Configuration Comparison - Cavallo Oswald Model 1');
legend(Design_Input.Config(:),'Location','northwest');
hold off

%Lift over Drag Analysis Plots (Configuration 1)
figure
hold on
plot(WingLiftCurve{1,:},LD_mod1{1,:},'--');
plot(WingLiftCurve{1,:},LD_mod2{1,:},'--');
plot(WingLiftCurve{1,:},LD_mod3{1,:},'--');
plot(WingLiftCurve{1,:},LD_benchmark{1,:});
xlabel('Coefficient of Lift (CL)');
ylabel('L/D Ratio (-)');
title('Lift over Drag Comparisons - Configuration 1.01');
legend('L/D Cavallo','L/D Nita-Scholz','L/D Benchmark','Location','southeast');
hold off

%% To Be provided at a later date

%Boost_Ascent Flight Profile Plots
%% Boost_Ascent Flight Profile Plots
% Some setup to make plots more readable in color, look up the
% documentation for 'cmap' for other color map options
cmap = colormap(parula(Count));
set(0,'DefaultAxesColorOrder',cmap)
set(gca(),'ColorOrder',cmap);

%2d plots
fields = fieldnames(stateStruct);
figure(20)
for n = 1:Count
    distBoost = vecnorm([stateStruct.(fields{n}).data(:, 4), stateStruct.(fields{n}).data(:, 5)], 2, 2);
    plot(distBoost,...
            -stateStruct.(fields{n}).data(:, 6), ...
            DisplayName=Design_Input.Properties.RowNames{n}, Color=cmap(n, :))
    if n == 1
        hold on
    end
end
xlabel('Total Distance Traveled [m]');
ylabel('Height Achieved [m]');
title('Boost 2D Total Distance Traveled');
legend();
grid on
hold off

figure(21)
for n = 1:Count
    plot(stateStruct.(fields{n}).data(:, 4),...
            stateStruct.(fields{n}).data(:, 5), ...
            DisplayName=Design_Input.Properties.RowNames{n}, Color=cmap(n, :))
    if n == 1
        hold on
    end
end
xlabel('y [m] - Positive = East');
ylabel('x [m] - Positive = North');
title('Boost Ground Track');
legend();
grid on
hold off

%%3D plots

figure(22)
for n = 1:Count    
    plot3(stateStruct.(fields{n}).data(:, 4),...
            stateStruct.(fields{n}).data(:, 5),...
            stateStruct.(fields{n}).data(:, 6), ...
            DisplayName=Design_Input.Properties.RowNames{n})
    if n == 1
        hold on
    end
end
xlabel('x [m] - Positive = North');
ylabel('y [m] - Positive = East');
zlabel('z [m] - Positive = Down');
title('Boost Trajectory Plots');
set(gca, 'ZDir','reverse')
set(gca, 'YDir','reverse')
legend();
grid on
axis equal
hold off

figure(23)
for n = 1:Count
    Wx = -Design_Input.V_wind(n)*cosd(Design_Input.Wind_Az(n)); 
    Wy = -Design_Input.V_wind(n)*sind(Design_Input.Wind_Az(n));
    quiver3(stateStruct.(fields{n}).data(:, 4), ... % x
        stateStruct.(fields{n}).data(:, 5), ... % y
        stateStruct.(fields{n}).data(:, 6), ... % z
        stateStruct.(fields{n}).data(:, 1)-Wx, ... % Vax
        stateStruct.(fields{n}).data(:, 2)-Wy, ... % Vay
        stateStruct.(fields{n}).data(:, 3), ... % Vz
        DisplayName=Design_Input.Properties.RowNames{n}, ...
        LineWidth=2)
    if n == 1
        hold on
    end
end
xlabel('x [m] - Positive = North');
ylabel('y [m] - Positive = East');
zlabel('z [m] - Positive = Down');
title('Boost Trajectory Plots with Heading Vetors');
set(gca, 'ZDir','reverse')
set(gca, 'YDir','reverse')
legend();
grid on
axis equal
hold off

%% Reset default color order
set(0,'DefaultAxesColorOrder','default')

GlideVel = GlideData.Vsink(1)/sind(GlideData.theta(1));

%Glide Flight Profile Plots

%Merged Boost-Ascent+Glide Flight Profile



