function [DragPolar_mod1,DragPolar_mod2,DragPolar_mod3] =...
    DragPolar(Parasite_Drag_Data,InducedDrag_Data,Design_Input,AoA_Count,WingLiftCurve,Count)
%% Drag Polar Summary
% Creates an array for each drag polar model with total CD value. 
% Columns are each configuration tested, rows are variation with angle of
% attack (-5 to 12 deg).
%
% Allows comparison of different configuration's drag polars (per model).
% Once a drag polar model is chosen, other models can be commented out if
% desired.

%% Outputs:
%
% DragPolar_mod1/2/3:
%   Table containing total drag data (parasite and induced) for each
%   induced drag model (1/2/3), each table has columns of AoA and rows of
%   case inputs

%% Preallocate variables of interest
% NOTE: These are being stored in a structure where the second level 
% variables are the different models. The arrays within this second level
% are the arrays discussed above
DragPolar_mod1 = zeros(Count,AoA_Count); 
DragPolar_mod2 = zeros(Count,AoA_Count);
DragPolar_mod3 = zeros(Count,AoA_Count);


%% Loop through different configurations
for n = 1:Count
% /////////////////////////////////////////////////////////////////////////
% MODIFY THIS SECTION 
% /////////////////////////////////////////////////////////////////////////
    DragPolar_mod1(n,:)= WingLiftCurve{Count,:}.^2 .* InducedDrag_Data.k1_mod1(n) + Parasite_Drag_Data.CDo(n);
    DragPolar_mod2(n,:)= WingLiftCurve{Count,:}.^2 .* InducedDrag_Data.k1_mod2(n) + Parasite_Drag_Data.CDo(n);
    DragPolar_mod3(n,:)= WingLiftCurve{Count,:}.^2 .* InducedDrag_Data.k1_mod3(n) + Parasite_Drag_Data.CDo(n);
% /////////////////////////////////////////////////////////////////////////
% END OF SECTION TO MODIFY
% /////////////////////////////////////////////////////////////////////////
end
%% Convert to tables for output
AoA_Names = {'-5', '-4', '-3', '-2', '-1', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'};
DragPolar_mod1 = array2table(DragPolar_mod1); % Convert to table
DragPolar_mod1.Properties.VariableNames = AoA_Names; % Name column headers for clarity using vector defined above
DragPolar_mod2 = array2table(DragPolar_mod2); 
DragPolar_mod2.Properties.VariableNames = AoA_Names;
DragPolar_mod3 = array2table(DragPolar_mod3);
DragPolar_mod3.Properties.VariableNames = AoA_Names;

end