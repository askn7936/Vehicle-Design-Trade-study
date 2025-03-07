function [LD_mod1,LD_mod2,LD_mod3,LD_benchmark] =...
   LD(Benchmark,DragPolar_mod1,DragPolar_mod2,DragPolar_mod3,WingLiftCurve,AoA_Count,Count)
%% Lift over Drag Analysis function Summary
% This function creates four arrays to compare the lift over drag (LD)
% estimations for three different drag polar models (mod1, mod2, mod3) with
% the benchmark aircraft's truth data L/D.  Each column in the LD.mod1,2,3
% tables represents a differnt configuration from the Design Input
% spreadsheet. Each row of the LD.mod1,2,3 tables is for the discrete
% angles of attack evaluated (from the WingLiftDrag function).
% Note that once a drag polar model is selected for use in your desgin, you
% may comment out the other models if desired.
%% Outputs:
%
% LD_mod1/2/3:
%   Table containing L/D values for each induced drag model (1/2/3) where
%   each table has columns of AoA and rows of case inputs
%
% LD_benchmark:
%   Table containing L/D values for the benchmark data where each table has
%   columns of AoA and rows of case inputs
%% Preallocate variables of interest
% NOTE: These are being stored in a structure where the second level
% variables are the different models. The arrays within this second level
% are the arrays discussed above
LD_mod1 = zeros(Count,AoA_Count);
LD_mod2 = zeros(Count,AoA_Count);
LD_mod3 = zeros(Count,AoA_Count);
LD_benchmark = zeros(1,AoA_Count);
liftarray = zeros(Count,AoA_Count);
mod1array = zeros(Count,AoA_Count);
mod2array = zeros(Count,AoA_Count);
mod3array = zeros(Count,AoA_Count);
%% Loop through different configurations
for n = 1:Count
% /////////////////////////////////////////////////////////////////////////
% MODIFY THIS SECTION
% /////////////////////////////////////////////////////////////////////////
    liftarray(n,:) = table2array(WingLiftCurve(n,:));
    mod1array(n,:) = table2array(DragPolar_mod1(n,:));
    mod2array(n,:) = table2array(DragPolar_mod2(n,:));
    mod3array(n,:) = table2array(DragPolar_mod3(n,:));
   LD_mod1(n,:)= liftarray(n,:) ./ mod1array(n,:);
   LD_mod2(n,:)= liftarray(n,:) ./ mod2array(n,:);
   LD_mod3(n,:)= liftarray(n,:) ./ mod3array(n,:);
% /////////////////////////////////////////////////////////////////////////
% END OF SECTION TO MODIFY
% /////////////////////////////////////////////////////////////////////////
end
%% Convert to tables for output
AoA_Names = {'-5', '-4', '-3', '-2', '-1', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'};
LD_mod1 = array2table(LD_mod1); % Convert to table
LD_mod1.Properties.VariableNames = AoA_Names; % Name column headers for clarity using vector defined above
LD_mod2 = array2table(LD_mod2);
LD_mod2.Properties.VariableNames = AoA_Names;
LD_mod3 = array2table(LD_mod3);
LD_mod3.Properties.VariableNames = AoA_Names;
%% Calculate the L/D for each coefficient of lift for the benchmark aircraft
LD_benchmark(1,:)=Benchmark.CL(:)./Benchmark.CD(:);
LD_benchmark = array2table(LD_benchmark); % Convert to table for consistency
LD_benchmark.Properties.VariableNames = AoA_Names;
end

