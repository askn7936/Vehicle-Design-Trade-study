function [InducedDrag_Data] =...
    InducedDrag(Design_Input,WingLiftModel,WingLiftCurve,WingDragCurve,WingGeo_Data,Count)
%% Induced Drag Model Function Summary
% This function evaluates different Oswalds Efficiency Factor models for 
% use in your drag polar model.  It compiles and outputs a variables table
% for three different Oswalds models (mod1, mod2, mod3).  For each Oswalds
% model, the Oswalds Efficiencty Factor (eo_mod1,2,3), and the k1 values
% (k1_mod1,2,3) are outputted in the InducedDrag_Data table.  Note that
% depending on the Oswalds models chosen to evaluate, you may or may not
% need information from teh WingGeo_Data table from the WingGeo function.
% Additionally, this code supports the calculation of the k2 values for
% evaluating non-symmetric airfoil design, but is not required.

%% Outputs:
%
% InducedDrag_Data:
%   Table containing lift curve slope, zero lift AoA, and oswalds info
%   (columns) for each input (rows)


%% Preallocate variables of interest
eo_mod1 = zeros(Count, 1); % Oswalds for Model #1
eo_mod2 = zeros(Count, 1); % Oswalds for Model #2
eo_mod3 = zeros(Count, 1); % Oswalds for Model #1
k1_mod1 = zeros(Count, 1); % k1 for Model #1
k1_mod2 = zeros(Count, 1); % k1 for Model #2
k1_mod3 = zeros(Count, 1); % k1 for Model #3
k2_mod1 = zeros(Count, 1); % k2 for Model #1
k2_mod2 = zeros(Count, 1); % k2 for Model #2
k2_mod3 = zeros(Count, 1); % k2 for Model #3

% NOTE: k2 values not required if only symmetric airfoils used; however,
% this version of the code includes it as an option

%% Loop through different configurations
for n = 1:Count 
% /////////////////////////////////////////////////////////////////////////
% MODIFY THIS SECTION
% /////////////////////////////////////////////////////////////////////////
    %Find CL min Drag value of wing drag polar to estimate k2
    [CD_min,CD_min_index] = min(WingDragCurve{n,:});
    CL_minD = WingLiftCurve{n,CD_min_index};

    %Cavallo Oswalds Model (Baseline required)
    eo_mod1(n) = 1.78*(1-0.045*Design_Input.AR_w(n)^0.68)-0.64;
    k1_mod1(n) =  1/(pi*eo_mod1(n)*Design_Input.AR_w(n));
    k2_mod1(n) = -2*k1_mod1(n)*CL_minD;

    %Student Option 1 Oswalds Model (Obert)

    eo_mod2(n) = 1 / (1.05 + (0.007 * pi * Design_Input.AR_w(n))); %Oswalds Estimate
    k1_mod2(n) = 1/(pi * eo_mod2(n) * Design_Input.AR_w(n));
    %k2_mod2(n) = ; % Optional
   
    %Student Option 2 Oswalds Model (DATCOM)
    e_invisid = 0.83; % Value is a range from 0.83 to 0.95
    P = 0.25; % Range from 0.25-0.45
    eo_mod3(n) = 1 / ((1/e_invisid) + P); %Oswalds Estimate
    k1_mod3(n) = 1/(pi*eo_mod3(n)*Design_Input.AR_w(n));
   % k2_mod3(n) = ; % Optional
% /////////////////////////////////////////////////////////////////////////
% END OF SECTION TO MODIFY
% /////////////////////////////////////////////////////////////////////////   
end

%% Organize into table for output
InducedDrag_Data = table(eo_mod1, eo_mod2, eo_mod3, k1_mod1, k1_mod2, k1_mod3, k2_mod1, k2_mod2, k2_mod3);

end
