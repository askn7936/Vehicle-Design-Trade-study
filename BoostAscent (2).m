function [apogee, hApogee, stateStruct] = BoostAscent(Design_Input, ATMOS, Parasite_Drag_Data, Weight_Data, ThrustCurves, Time, Count)
%BOOSTASCENT
% The main purpose of this function is to set up for the call to ODE 45 to
% propigate our eqations of motion, and then to format the output of ODE 45
% to be passed back into main for the glide function and for plotting.

%% Outputs:
%
% apogee:
%   A vector of values of the maximum height achieved in [m] for each case 
%   input
%
% hApogee:
%   A table of inertial heading vectors in cartisian coordinates for each
%   case input
%
% stateStruct:
%   A structure containing the full output data from ODE45 for each case
%   input. This structure's second level will correspond to the different
%   input cases, with the data being containied within being the table of
%   states for all time steps. This is meant only for in-depth analysis and
%   is not meant to be munipulated by the students

%% Preallocate variables of interest
apogee = zeros(Count,1);
hApogee = zeros(Count,3);

%% Preallocate ODE varaibles
consts = zeros(12, 1);
S0 = zeros(7, 1);

%% Some useful constants
g = 9.8; % Accelatation due to gravity [m/s^2]
rho_w = 1000; % Density of water [kg/m^3]
mu_k = 0.2; % Launch rail coefficient of dynamic friction []
A_exit = 0.021; % Area of bottle outlet [m^2]

%% Loop through different configurations
for n = 1:Count

    %% Pick the correct thrust curve
    waterSize = Design_Input.Water_Vol(n);
    bottleSize = Design_Input.Bottle_Vol(n)*1000; % Convert to [ml]

    thrustCurveName = [num2str(bottleSize), '_', num2str(waterSize)];

    thrustVec = ThrustCurves.(thrustCurveName); % use parenthasis to index into a table with a variable

    % /////////////////////////////////////////////////////////////////////////
    % MODIFY THIS SECTION
    % /////////////////////////////////////////////////////////////////////////
    %% Make variables for the constants vector
    m_empty = (Weight_Data.Wo(n)-Weight_Data.W_water(n))/g; % [kg] empty weight of the vehicle (not including the water)
    m0 = Weight_Data.Wo(n)/g; % Note that the input water volume should be in ml which approx = grams
    Wx = Design_Input.V_wind(n)*cosd(Design_Input.Wind_Az(n)); % [m/s] Wind velocity in the x-direction
    Wy = Design_Input.V_wind(n)*sind(Design_Input.Wind_Az(n)); % [m/s] Wind velocity in the y-direction

    % /////////////////////////////////////////////////////////////////////////
    % END OF SECTION TO MODIFY
    % /////////////////////////////////////////////////////////////////////////
    %% Pack Constants Vector
    % Basic Properties
    consts(1) = g; % Accelatation due to gravity [m/s^2]
    consts(2) = rho_w; % Density of water [kg/m^3]
    consts(3) = ATMOS.rho(n); % Density of air [kg/m^3]
    consts(4) = mu_k; % Launch rail coefficient of dynamic friction []
    % Vehicle info
    consts(5) = A_exit; % Area of bottle outlet [m^2]
    consts(6) = Parasite_Drag_Data.CDo(n); % C_D of the vehicle (assume zero lift) []
    consts(7) = Design_Input.Sref_w(n); % Wing reference area [m^2]
    consts(8) = m_empty; % Weight of the rocet with no water [kg]
    % Wind
    consts(9) = Wx; % Inertial wind velocity in x [m/s]
    consts(10) = Wy; % Inertial wind velocity in x [m/s]
    % Launch Direction
    consts(11) = Design_Input.Launch_El(n); % Launch elevation [degrees]
    consts(12) = Design_Input.Launch_Az(n); % Launch Azimuth [degrees], measured CW from north when looking down on the map; also known as compass heading

    %% Make an initial condition state vector
    S0(1) = 0; % inertial velocity in x-direction [m/s]
    S0(2) = 0; % inertial velocity in y-direction [m/s]
    S0(3) = 0; % inertial velocity in z-direction [m/s]
    S0(4) = 0; % position in x (inertial) [m]
    S0(5) = 0; % position in y (inertial) [m]
    S0(6) = 0; % position in z (inertial) [m]
    S0(7) = m0; % current total mass [kg]

    %% Final setup
    % set and event to stop the integration 
    opts = odeset('Events', @stoppingPoint, 'RelTol',1e-8);

    % Set a start and max time to integrate over (pick one of the below)
    intSpan = [0, 5]; % [s] - Lets matlab pick its timestep (reccomended)
    % intSpan = 0:0.001:5; % [s] - forces output to have 1ms spacing

    %% Call to ODE 45
    % The @ part is telling ODE45 what our independent vaiable is (t) and
    % what the dependet variables to propigate in time are (S)
    %
    % Then, we pass it the odefun we made that can take in the current
    % time, state, and whatever else it requires and retuns the derivative
    % of the state variables at the current time
    %
    % Finally, we also give ODE45 a span of time to integrate over
    % (intSpan), an initial conditions vector (S0), and any options we set
    % (opts), which in this case holds the function to stop the integration
    [t, S] = ode45(@(t, S) BoostAscent_odefun(t, S, consts, thrustVec, Time), intSpan, S0, opts);

    %% Find outputs of interest
    [apogee(n), iApogee] = max(-S(:, 6));
    hApogee(n, :) = S(iApogee, 4:6)/norm(S(iApogee, 4:6));

    %% Store full output
    configName = ['Config_' num2str(n)];
    stateStruct.(configName).time = t;
    stateStruct.(configName).data = S;

end
%% Convert to tables for output
dirNames = {'x', 'y', 'z'};
hApogee = array2table(hApogee); % Convert to table
hApogee.Properties.VariableNames = dirNames;

   