function [Weight_Data,CG_Data] = Weight(Design_Input,Count,WingGeo_Data,Airfoil,Material_Data)
%% Weight Model Summary 
% This function pulls in the material properties from the Data Input file
% (Material_Data table) and creates a weight model estimate for a given
% aircraft configuration based on the Design_Input, WingGeo_Data, and
% Airfoil tables.  Note that you will need to vary this model depending on
% how you choose to fabricate your prototype as the material choices coded
% below may vary for your prototype.

%% Outputs:
%
% Weight_Data:
%   Table containing total weight and a component breakdown of
%   contributions to that total (columns), for each input case (rows)
%
% CG_Data:
%   Table containing overall CG and a component breakdown of contributions
%   to that total (columns), for each input case (rows)s

%% Preallocate variables of interest
Wo = zeros(Count, 1); % Total Weight [N]
W_f = zeros(Count, 1); % Fusalage weight [N]
W_w = zeros(Count, 1); % Wing weight [N]
W_h1 = zeros(Count, 1); % Horizontal tail 1 weight [N]
W_h2 = zeros(Count, 1); % Horizontal tail 2 weight [N]
W_v1 = zeros(Count, 1); % Vertical tail 1 weight [N]
W_v2 = zeros(Count, 1); % Vertical tail 2 weight [N]
W_bottle = zeros(Count, 1); % Bottle weight [N]
W_water = zeros(Count, 1); % Water weight [N]
W_pay = zeros(Count, 1); % Payload weight [N]
W_ballast = zeros(Count, 1); % Ballast Weight [N]

CG_tot = zeros(Count, 1); % Overall CG [m]
CG_f = zeros(Count, 1); % Fusalage CG [m]
CG_w = zeros(Count, 1); % Wing CG [m]
CG_h1 = zeros(Count, 1); % Horizontal tail 1 CG [m]
CG_h2 = zeros(Count, 1); % Horizontal tail 2 CG [m]
CG_v1 = zeros(Count, 1); % Vertical tail 1 CG [m]
CG_v2 = zeros(Count, 1); % Vertical tail 2 CG [m]
CG_bottle = zeros(Count, 1); % Bottle CG [m]
CG_water = zeros(Count, 1); % Water CG [m]
CG_pay = zeros(Count, 1); % Payload CG [m]
CG_ballast = zeros(Count, 1); % Ballast CG [m]


%% Loop through different configurations
for n = 1:Count
% /////////////////////////////////////////////////////////////////////////
% MODIFY THIS SECTION
% /////////////////////////////////////////////////////////////////////////
    %Fuselage Weight & CG Estimate
    W_f(n)=Material_Data.rhoA_foamboard(1)*pi*(Design_Input.Dia_f(n)^2)*Design_Input.Length_f(n)*9.81; %Assumes fuselage created from foamboard; units of Newtons
    CG_f(n) = Design_Input.Length_f(n)/2; %Approximate location of fuselage CG from nose (m)

    %Wing Weight Estimate
    W_w(n)=Material_Data.rho_pinkfoam(1)*Design_Input.Sref_w(n)*Airfoil.Thick_w(n)*WingGeo_Data.MAC_w(n)*9.81; %Assumes pink foam; units of Newtons
    CG_w(n) = 0;

    %Horz Tail 1 Weight Estimate
    W_h1(n)=Material_Data.rhoA_foamboard(1)*Design_Input.Sref_h1(n)*9.81; %Assumes horz tail made from foamboard, units of Newtons
    CG_h1(n) = 0;

    %Horz Tail 2 Weight Estimate
    W_h2(n)=Material_Data.rhoA_foamboard(1)*Design_Input.Sref_h2(n)*9.81; %Assumes horz tail made from foamboard, units of Newtons
    CG_h2(n) = 0;

    %Vert Tail 1 Weight Estimate
    W_v1(n)=Material_Data.rhoA_foamboard(1)*Design_Input.Sref_v1(n)*9.81; %Assumes vert tail made from foamboard, units of Newtons
    CG_v1(n) = 0;

    %Vert Tail 2 Weight Estimate
    W_v2(n)=Material_Data.rhoA_foamboard(1)*Design_Input.Sref_v2(n)*9.81; %Assumes vert tail made from foamboard, units of Newtons
    CG_v2(n) = 0;

    %Pressure vessel (Coke Bottle) Weight
    if Design_Input.Bottle_Vol(n)==1.25
        W_bottle(n) = (0.048/1000)*9.81; %1.25 liter bottle direct measurement 50 grams to Newtons
        CG_bottle(n) = 0; %1.25 liter bottle direct measurement (Newtons)
    else
        W_bottle(n) = (0.063/1000)*9.81; %2 liter bottle direct measurement 68 grams to Newtons
        CG_bottle(n) = 0; %2 liter bottle direct measurement (Newtons)
    end
    
    %Water Propellant Weight
    W_water(n) = (Design_Input.Water_Vol(n)/1000)*9.81; %Water propellant weight direct measurement (Newtons)
    CG_water(n) = 0;

    %Payload Weight (as required)
    W_pay(n) = 0.05*9.81; %Payload weight as required (Netwtons)
    CG_pay(n) = 0;

    %Ballast Weight (as required)
    W_ballast(n) = 0; %Ballast Weight as required (Newtons)
    CG_ballast(n) = 0;

    %Total Weight & CG Location
    Wo(n)=W_f(n)+W_w(n)+W_h1(n)+W_h2(n)+W_v1(n)+W_v2(n)+W_bottle(n)+W_water(n)+W_pay(n)+W_ballast(n);
    CG_tot(n)=(W_f(n)*CG_f(n)+W_w(n)*CG_w(n)+W_h1(n)*CG_h1(n)+W_h2(n)*CG_h2(n)+W_v1(n)*CG_v1(n)+W_v2(n)*CG_v2(n)+W_bottle(n)*CG_bottle(n)+W_water(n)*CG_water(n)+W_pay(n)*CG_pay(n)+W_ballast(n)*CG_ballast(n))/Wo(n);
% /////////////////////////////////////////////////////////////////////////
% END OF SECTION TO MODIFY
% /////////////////////////////////////////////////////////////////////////
end

%% Oraganize into tables for output
Weight_Data= table(Wo, W_f, W_w, W_h1, W_h2, W_v1, W_v2, W_bottle, W_water, W_pay, W_ballast);
CG_Data = table(CG_tot, CG_f, CG_w, CG_h1, CG_h2, CG_v1, CG_v2, CG_bottle, CG_water, CG_pay, CG_ballast);

end
