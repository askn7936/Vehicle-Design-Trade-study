function [WingLiftModel,AoA,AoA_Count,AirfoilLiftCurve,WingLiftCurve,WingDragCurve] =...
   WingLiftDrag(Design_Input,Airfoil,Count)
%% Lift Model Function Summary:
% This function captures and outputs the conversion variables of a 2D
% airfoil's lift curve to a 3D wing.  Specifically, it takes in 2D airfoil
% data from the user inputed "Aifoil" tab in the Design Input spreadsheet,
% and calculates the 2D lift curve slope (a_o), the 3D lift curve slope (a),
% the zero lift angle of attack (AoA_o), and the span efficiency factor (e)
% for every configuration airfoil choice put into teh Design Input
% spreadsheet.  This function also establishes the values of the discrete
% angles of attack that will be evaluated (AoA), and the number (AoA_Count)
% to aid in subsequent functions that will required this information.
% The default is from -5 deg to 12 deg AoA and is sufficient for most/all
% airfoils chosen.
%% Outputs:
%
% WingLiftModel:
%   Table containing lift curve slope, zero lift AoA, and span efficiency
%   info (columns) for each input case (rows)
%
% AoA:
%   1D Array contatining AoA's considered in calculations
%
% AoA_Count:
%   Value containing length(AoA)
%
% AirfoilLiftCurve:
%   Table with data representing the 2D lift curve for each input case
%   (rows) and for each AoA (columns)
%
% WingLiftCurve:
%   Table with data representing the 3D lift curve for each input case
%   (rows) and for each AoA (columns)
%
% WingDragCurve:
%   Table with data representing the 3D drag curve (for the wing only) for
%   each input case (rows) and for each AoA (columns)
%% Preallocate variables of interest
a_0 = zeros(Count, 1); % 2D airfoil lift curve slope
a = zeros(Count, 1); % 3D airfoil lift curve slope
AoA_0 = zeros(Count, 1); % 2D airfoil zero lift AoA [deg]
e = zeros(Count, 1); % Span efficiency factor
%% Initial calculations
AoA = -5:1:12; % Range of AoA considered (-5 to 12 deg increment by 1 deg) - Must match range of AoA's in 2D lift curve slope data
AoA_Count = length(AoA); % Count number of angles considered
%% Preallocate Other data we will need now that we have our AoA Array
AirfoilLiftCurve = zeros(Count,AoA_Count);
WingLiftCurve = zeros(Count,AoA_Count);
WingDragCurve = zeros(Count,AoA_Count);
%% Loop through different configurations
for n = 1:Count
% /////////////////////////////////////////////////////////////////////////
% MODIFY THIS SECTION
% /////////////////////////////////////////////////////////////////////////
   % Determine Span Efficiency Factor (e)
   if Design_Input.QuarterSweep_w(n) == 0 % If no quarter chord sweep:
       f_taper = 0.0524*Design_Input.Taper_w(n)^4-0.15*Design_Input.Taper_w(n)^3+0.1659*Design_Input.Taper_w(n)^2-0.0706*Design_Input.Taper_w(n)+0.0119;
       e(n)= 1/(1+f_taper*Design_Input.AR_w(n));
   else % If there is quarter chord sweep:
       f_taper = 0.0524*Design_Input.Taper_w(n)^4-0.15*Design_Input.Taper_w(n)^3+0.1659*Design_Input.Taper_w(n)^2-0.0706*Design_Input.Taper_w(n)+0.0119; ;
       e(n) = 1/ (1+f_taper*Design_Input.AR_w(n)) .* cosd( Design_Input.QuarterSweep_w(n) );
   end
   % Linear fit to Airfoil lift data
   Cl = polyfit(AoA(1:12),Airfoil{n,(5:16)},1); % Reffer to MATLAB documentation for this function's usage
   a_0(n)=Cl(1); % 2D airfoil lift curve slope
   AoA_0(n)=roots(Cl); % 2D airfoil zero lift AoA (deg)
   AirfoilLiftCurve(n,:) = polyval(Cl,AoA);
   % 3D Wing Lift Curve Slope Model
   a(n)= a_0(n)  ./ ( 1 + (57.3.*a_0(n))/(pi.*e(n).*Design_Input.AR_w(n)))  ; % 3-D lift curve slope
   WingLiftCurve(n,:) = a(n).*(AoA-AoA_0(n)); % Evaluate the 3D coefficient of lift values at AoA's in array (AoA) definied above
   % 3D Wing Drag Coefficient
   Cd = Airfoil{n,24:41};
   CD = Cd + (1/(pi .* e(n).* Design_Input.AR_w(n))) .* WingLiftCurve(n,:).^2;
   WingDragCurve(n,:) = CD(:); % Just the 3D coefficient of drag values at AoA's in array (AoA) definied above
% /////////////////////////////////////////////////////////////////////////
% END OF SECTION TO MODIFY
% /////////////////////////////////////////////////////////////////////////
end
%% Organize into table for output
WingLiftModel = table(a_0, a, AoA_0, e);
%% Convert arrays to tables for clarity
AoA_Names = {'-5', '-4', '-3', '-2', '-1', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'};
AirfoilLiftCurve = array2table(AirfoilLiftCurve);
AirfoilLiftCurve.Properties.VariableNames = AoA_Names; % Name column headers for clarity using vector defined above
WingLiftCurve = array2table(WingLiftCurve);
WingLiftCurve.Properties.VariableNames = AoA_Names;
WingDragCurve = array2table(WingDragCurve);
WingDragCurve.Properties.VariableNames = AoA_Names;
end

