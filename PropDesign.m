%% 6. PropDesign (Sizing and Optimizing Propeller for Mission Requirements (+ Hover))

% Call 5. Propulsion
Propulsion;
propIteration = "22";
propSweepsTemp = {'0rseqHVc.txt', '0rseqHS.txt', '0rseqC.txt';
                  '1rseqHVc.txt', '1rseqHS.txt', '1rseqC.txt';
                  '2rseqHVc.txt', '2rseqHS.txt', '2rseqC.txt';
                  '2rseqHVc.txt', '3rseqHS.txt', '3rseqC.txt';
                  '4rseqHVc.txt', '4rseqHS.txt', '4rseqC.txt';
                  '5rseqHVc.txt', '5rseqHS.txt', '5rseqC.txt';
                  '6rseqHVc.txt', '6rseqHS.txt', '6rseqC.txt';
                  '7rseqHVc.txt', '7rseqHS.txt', '7rseqC.txt';
                  '8rseqHVc.txt', '8rseqHS.txt', '8rseqC.txt';
                  '9rseqHVc.txt', '9rseqHS.txt', '9rseqC.txt';
                  '10rseqHVc.txt', '10rseqHS.txt', '10rseqC.txt';
                  '11rseqHVc.txt', '11rseqHS.txt', '11rseqC.txt';
                  '12rseqHVc.txt', '12rseqHS.txt', '12rseqC.txt';
                  '13rseqHVc.txt', '13rseqHS.txt', '13rseqC.txt';
                  '14rseqHVc.txt', '14rseqHS.txt', '14rseqC.txt';
                  '15rseqHVc.txt', '15rseqHS.txt', '15rseqC.txt';
                  '16rseqHVc.txt', '16rseqHS.txt', '16rseqC.txt';
                  '17rseqHVc.txt', '17rseqHS.txt', '17rseqC.txt';
                  '18rseqHVc.txt', '18rseqHS.txt', '18rseqC.txt';
                  '19rseqHVc.txt', '19rseqHS.txt', '19rseqC.txt';
                  '20rseqHVc.txt', '20rseqHS.txt', '20rseqC.txt';
                  '21rseqHVc.txt', '21rseqHS.txt', '21rseqC.txt';
                  '22rseqHVc.txt', '22rseqHS.txt', '22rseqC.txt';
               };
prefix = 'XROTORSweepData/';
propSweeps = cellfun(@(x) [prefix, x], propSweepsTemp, 'UniformOutput', false);

varTableProp = readtable("ThesisWingOptimizationVariables.xlsx", 'Sheet', 'PropellerOptimizationVariables', 'ReadVariableNames', true, 'VariableNamingRule', 'preserve');
iterationColProp = find(strcmp(varTableProp.Properties.VariableNames, propIteration));

% derived from T = 0.5*rho*A*V^2
V_airspeed_Cruise = V; % m/s
V_airspeed_climb = 1; % m/s
V_climb_assessed = 1;

% v_i_hover = 15; % m/s
propRadius_min_Cruise = sqrt(2*T_cruise / (pi*numEDF*rho_sel*V_airspeed_Cruise^2));

v_i_climb = V_climb_assessed / (sqrt(T_W_ratio_climb^2 - 1));
lambda = 1.28;  % Ducted fan thrust augmentation factor
propRadius_climb_ducted = sqrt(1/lambda * T_climb / ...
    (2*rho_sel*numEDF*pi*(v_i_climb + (V_climb_assessed/2))^2 - (V_climb_assessed/2)^2));

disp('Parameters for XRotor Design: ');
disp(['Max Unit Thrust from each EDF: ', num2str(T_mission_max_unit), ' m'])
disp(['Minimum Radius of Propeller: ', num2str(propRadius_climb_ducted), ' m'])

%% Move XRotor Prop Parameters to Excel Table
% Define input and output files
inputFile = 'DuctedPotentialFormulationSolutions.txt';
outputFile = 'ThesisWingOptimizationVariables.xlsx';
sheetName = 'PropellerOptimizationVariables';

% Desired labels
rowLabels = {'Vdisk/Vslip', 'Wake adv. ratio', 'no. blades', 'radius(m)', ...
    'adv. ratio' , 'thrust(N)' , 'power(W)' , 'Torque (N-m)' , 'Efficiency', ...
    'speed(m/s)', 'rpm', 'Eff induced', 'Eff ideal', 'Tcoef', 'Tnacel(N)', ...
    'hub rad.(m)', 'disp. rad.', 'Tvisc(N)', 'Pvisc(W)', 'rho(kg/m3)',...
    'Vsound(m/s)', 'mu(kg/m-s)', 'Ct', 'Cp', 'J', 'Tc', 'Pc', 'adv'};

colNames = {'Vdisk/Vslip', 'Wake adv. ratio', 'no. blades', 'radius(m)',... 
    'adv. ratio', 'thrust(N)', 'power(W)', '.  torque(N-m)', 'Efficiency',...
    'speed(m/s)', 'rpm', 'Eff induced', 'Eff ideal', 'Tcoef', 'Tnacel(N)',...
    'hub rad.(m)', 'disp. rad.', 'Tvisc(N)', 'Pvisc(W)', 'rho(kg/m3)',...
    'Vsound(m/s)', 'mu(kg/m-s)', 'Ct', 'Cp', 'J', 'Tc', 'Pc', 'adv'};

iterTableVarNames = {'Vdisk_slip', 'wake_J', 'numBlades', 'radius_m', 'J' ,...
    'T_N' , 'P_W' , 'torque_Nm' , 'eff', 'speed', 'rpm', 'eff_induced', ...
    'eff_ideal', 'Tcoef', 'Tnacel_N', 'radiusHub', 'disp_rad', 'Tvisc_N', ...
    'Pvisc_W', 'rho', 'Vsound', 'mu', 'Ct', 'Cp', 'J', 'Tc', 'Pc', 'adv'};


keyMap = containers.Map(rowLabels, colNames);
fileText = fileread(inputFile);
blocks = regexp(fileText, 'Ducted Potential Formulation Solution:.*?(?=Ducted Potential Formulation Solution:|$)', 'match');
currPropIteration = str2double(propIteration) + 1;
block = blocks{currPropIteration};
% Fix formatting for consistent key:value matching and match all "Label:
% Value" pairs
block = regexprep(block, '(\d)\s+([a-zA-Z\/\.\(\)%\- ]+?:)', '$1\n$2');
matches = regexp(block, '([a-zA-Z][a-zA-Z0-9 _\./\-\(\)%]*):\s*(-?\d*\.?\d+(?:[eE][-+]?\d+)?)', 'tokens');

% Build map of parsed data
dataMap = containers.Map();
for j = 1:length(matches)
    key = strtrim(matches{j}{1});
    val = str2double(matches{j}{2});
    if isKey(keyMap, key)
        dataMap(keyMap(key)) = val;
    else
        dataMap(key) = val;
    end
end

% Make values array and iterTable and fill excel table
values = nan(length(rowLabels), 1);
iterTable = struct();
for i = 1:length(rowLabels)
    label = rowLabels{i};
    if isKey(dataMap, label)
        values(i) = dataMap(label);
    end

    iterTable.(iterTableVarNames{i}) = values(i);

    rowMatch = strcmp(varTableProp.VariableName, label);
    if any(rowMatch)
        varTableProp{rowMatch, iterationColProp} = values(i);
    else
        newRow = cell(1, width(varTableProp));
        newRow{1} = label;
        newRow{iterationColProp} = values(i);
        varTableProp = [varTableProp; newRow];
    end
end

writetable(varTableProp, outputFile, 'Sheet', sheetName, 'WriteRowNames', true);


% Sanity Check to confirm that values from .txt file conform
propDia = iterTable.radius_m * 2; % m
n = iterTable.speed / (propDia* iterTable.J); % RPS
rpm_calc = n*60; % RPM
T_calc  = iterTable.Ct * iterTable.rho * n^2 * propDia^4;
P_calc = iterTable.Cp * iterTable.rho * n^3 * propDia^5;
fprintf('RPM Output: %.2f RPM, Calculated: %.2f RPM\n', iterTable.rpm, rpm_calc);
fprintf('Power Output: %.2f W, Calculated: %.2f W\n', iterTable.P_W, P_calc);
fprintf('Thrust Output: %.2f N, Calculated: %.2f N\n', iterTable.T_N, T_calc);
disp('--')

%% XRotor Analaysis
radius_sim = varTableProp{strcmp(varTableProp.VariableName, "radius(m)"), propIteration}; 

[~, numSweepsPerIteratios] = size(propSweeps);
currPropIteration = str2double(propIteration)+1;
filenameSetProp = propSweeps(currPropIteration, :); 

colMap = cell(1, numel(filenameSetProp));
dataAll = cell(1, numel(filenameSetProp)); 

for i = 1:numel(filenameSetProp)
    fid = fopen(filenameSetProp{i}, 'r');
    header = '';
    while ~feof(fid)
        line = fgetl(fid);
        if contains(line, 'n') && contains(line, 'eff')
            header = strtrim(line);
            break;
        end
    end
    while ~feof(fid)
        line = fgetl(fid);
        if contains(line, '')
            break;
        end
    end

    colNames = strsplit(strtrim(header));
    colMap{i} = containers.Map(colNames, 1:length(colNames));
    data = [];
    expectedCols = length(colNames);
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if isempty(line)
            continue;
        end
        rowVals = sscanf(line, '%f');
        if length(rowVals) ~= expectedCols
            warning('Skipping malformed row in %s: %s', filenameSetProp{i}, line);
            continue;
        end
        data = [data; rowVals'];  
    end

    fclose(fid);
    dataAll{i} = data;
end

% Extract from first file (rseqH)
colMap_rseqHVc = colMap{1};
colMap_rseqHS   = colMap{2};
colMap_rseqC = colMap{3};
%colMap_vseq = colMap{4};
%colMap_aseq = colMap{5};
data_rseqHVc   = dataAll{1};
data_rseqHS   = dataAll{2};
data_rseqC   = dataAll{3};

%data_vseq   = dataAll{4};
%data_aseq   = dataAll{5};

% Named column arrays
n_rseqHVc       = data_rseqHVc(:, colMap_rseqHVc('n'));
V_wR_rseqHVc    = data_rseqHVc(:, colMap_rseqHVc('V/wR'));
Btip_rseqHVc    = data_rseqHVc(:, colMap_rseqHVc('Btip'));
V_rseqHVc       = data_rseqHVc(:, colMap_rseqHVc('V'));
rpm_rseqHVc     = data_rseqHVc(:, colMap_rseqHVc('rpm'));
rho_rseqHVc     = data_rseqHVc(:, colMap_rseqHVc('rho'));
mu_rseqHVc      = data_rseqHVc(:, colMap_rseqHVc('mu*1e5')) * 1e-5;  % Convert back to actual units
Vsound_rseqHVc  = data_rseqHVc(:, colMap_rseqHVc('Vsound'));
h_rseqHVc       = data_rseqHVc(:, colMap_rseqHVc('h'));
P_rseqHVc       = data_rseqHVc(:, colMap_rseqHVc('P(kW)')) * 1e3;    % Convert kW to W
T_rseqHVc       = data_rseqHVc(:, colMap_rseqHVc('T(N)'));
Q_rseqHVc       = data_rseqHVc(:, colMap_rseqHVc('Q(N-m)'));
eff_rseqHVc     = data_rseqHVc(:, colMap_rseqHVc('eff'));

n_rseqHS       = data_rseqHS(:, colMap_rseqHS('n'));
V_wR_rseqHS    = data_rseqHS(:, colMap_rseqHS('V/wR'));
Btip_rseqHS    = data_rseqHS(:, colMap_rseqHS('Btip'));
V_rseqHS       = data_rseqHS(:, colMap_rseqHS('V'));
rpm_rseqHS     = data_rseqHS(:, colMap_rseqHS('rpm'));
rho_rseqHS     = data_rseqHS(:, colMap_rseqHS('rho'));
mu_rseqHS      = data_rseqHS(:, colMap_rseqHS('mu*1e5')) * 1e-5;  % Convert back to actual units
Vsound_rseqHS  = data_rseqHS(:, colMap_rseqHS('Vsound'));
h_rseqHS       = data_rseqHS(:, colMap_rseqHS('h'));
P_rseqHS       = data_rseqHS(:, colMap_rseqHS('P(kW)')) * 1e3;    % Convert kW to W
T_rseqHS       = data_rseqHS(:, colMap_rseqHS('T(N)'));
Q_rseqHS       = data_rseqHS(:, colMap_rseqHS('Q(N-m)'));
eff_rseqHS     = data_rseqHS(:, colMap_rseqHS('eff'));

n_rseqC       = data_rseqC(:, colMap_rseqC('n'));
V_wR_rseqC    = data_rseqC(:, colMap_rseqC('V/wR'));
Btip_rseqC    = data_rseqC(:, colMap_rseqC('Btip'));
V_rseqC       = data_rseqC(:, colMap_rseqC('V'));
rpm_rseqC     = data_rseqC(:, colMap_rseqC('rpm'));
rho_rseqC     = data_rseqC(:, colMap_rseqC('rho'));
mu_rseqC      = data_rseqC(:, colMap_rseqC('mu*1e5')) * 1e-5;  % Convert back to actual units
Vsound_rseqC  = data_rseqC(:, colMap_rseqC('Vsound'));
h_rseqC       = data_rseqC(:, colMap_rseqC('h'));
P_rseqC       = data_rseqC(:, colMap_rseqC('P(kW)')) * 1e3;    % Convert kW to W
T_rseqC       = data_rseqC(:, colMap_rseqC('T(N)'));
Q_rseqC       = data_rseqC(:, colMap_rseqC('Q(N-m)'));
eff_rseqC     = data_rseqC(:, colMap_rseqC('eff'));

%{
n_vseq       = data_vseq(:, colMap_vseq('n'));
V_wR_vseq    = data_vseq(:, colMap_vseq('V/wR'));
Btip_vseq    = data_vseq(:, colMap_vseq('Btip'));
V_vseq       = data_vseq(:, colMap_vseq('V'));
rpm_vseq     = data_vseq(:, colMap_vseq('rpm'));
rho_vseq     = data_vseq(:, colMap_vseq('rho'));
mu_vseq      = data_vseq(:, colMap_vseq('mu*1e5')) * 1e-5;  % Convert back to actual units
Vsound_vseq  = data_vseq(:, colMap_vseq('Vsound'));
h_vseq       = data_vseq(:, colMap_vseq('h'));
P_vseq       = data_vseq(:, colMap_vseq('P(kW)')) * 1e3;    % Convert kW to W
T_vseq       = data_vseq(:, colMap_vseq('T(N)'));
Q_vseq       = data_vseq(:, colMap_vseq('Q(N-m)'));
eff_vseq     = data_vseq(:, colMap_vseq('eff'));

n_aseq       = data_aseq(:, colMap_aseq('n'));
V_wR_aseq    = data_aseq(:, colMap_aseq('V/wR'));
Btip_aseq    = data_aseq(:, colMap_aseq('Btip'));
V_aseq       = data_aseq(:, colMap_aseq('V'));
rpm_aseq     = data_aseq(:, colMap_aseq('rpm'));
rho_aseq     = data_aseq(:, colMap_aseq('rho'));
mu_aseq      = data_aseq(:, colMap_aseq('mu*1e5')) * 1e-5;  % Convert back to actual units
Vsound_aseq  = data_aseq(:, colMap_aseq('Vsound'));
h_aseq       = data_aseq(:, colMap_aseq('h'));
P_aseq       = data_aseq(:, colMap_aseq('P(kW)')) * 1e3;    % Convert kW to W
T_aseq       = data_aseq(:, colMap_aseq('T(N)'));
Q_aseq       = data_aseq(:, colMap_aseq('Q(N-m)'));
eff_aseq     = data_aseq(:, colMap_aseq('eff'));



rhos = [rho_rseqC(1), rho_rseqHVc(1), rho_rseqHS(1), rho_aseq(1), rho_vseq(1)];

if ~all(rhos == rhos(1)) || ~all(rho_rseqC == rho_rseqC(1)) || ~all(rho_rseqHVc == rho_rseqHVc(1)) || ~all(rho_rseqHS == rho_rseqHS(1)) || ~all(rho_vseq == rho_vseq(1)) || ~all(rho_aseq == rho_aseq(1))
    error('Rho must be consistent!');
end

%}

% Preparing for Interpolation to find power at target thrusts
[T_sorted, sortIdx] = sort(T_rseqC);
P_sorted = P_rseqC(sortIdx);
eff_sorted = eff_rseqC(sortIdx);

% Find RPM for flight conditions
Q_cruise = interpolate(T_cruise_unit, T_rseqC, Q_rseqC);
Q_Vc = interpolate(T_climb_unit, T_rseqHVc, Q_rseqHVc);
Q_HS = interpolate(T_hover_min_unit, T_rseqHS, Q_rseqHS); % T = W

% Find Q for flight conditions
cruiseRPM = interpolate(T_cruise_unit, T_rseqC, rpm_rseqC); % N-m
hoverRPM_Vc = interpolate(T_climb_unit, T_rseqHVc, rpm_rseqHVc); % N-m
hoverRPM_Still = interpolate(T_hover_min_unit, T_rseqHS, rpm_rseqHS); % N-m

MotorSelection;

% Cruise Data
powerCruise_unit = interpolate(cruiseRPM, rpm_rseqC, P_rseqC);
J_cruise = findJ(V, cruiseRPM, radius_sim); % will be 0 for purely vertical hover
prop_E_cruise = interpolate(cruiseRPM, rpm_rseqC, eff_rseqC);
power_true_unit_cruise = powerCruise_unit / motor_E_curr_cruise;

% P_rseq already includes the effect of low efficiency, because XROTOR computes 
% the actual work required to achieve thrust T at velocity V, given the blade aerodynamics.
% So no prop_E to determine true cruise power
power_true_total_Cruise = power_true_unit_cruise * numEDF;

% Hover Data

v_i_hoverStill = sqrt((T_hover_min_unit)/(2*rho_rseqHS(1)*(radius_sim^2*pi)));
v_i_hoverVc = sqrt((V_climb_assessed / 2)^2 + T_climb_unit / (2 * rho_rseqHVc(1) * (pi * radius_sim^2))) - V_climb_assessed / 2;

powerActual_unit_hoverVc = (interpolate(hoverRPM_Vc, rpm_rseqHVc, P_rseqHVc)) / motor_E_curr_hoverVc;
powerActual_unit_hoverStill = (interpolate(hoverRPM_Still, rpm_rseqHS, P_rseqHS)) / motor_E_curr_hoverStill;

powerIdealUnit_hoverVc = idealHoverPower(T_climb_unit, V_climb_assessed, rho_rseqHVc(1), radius_sim);
powerIdealUnit_hoverStill = idealHoverPower(T_hover_min_unit, 0, rho_rseqHS(1), radius_sim);

FM_hoverVc = powerIdealUnit_hoverVc / powerActual_unit_hoverVc;
FM_hoverStill = powerIdealUnit_hoverStill / powerActual_unit_hoverStill;

power_true_total_hoverVc = powerActual_unit_hoverVc * numEDF;
power_true_total_hoverStill = powerActual_unit_hoverStill * numEDF;

% Sweep number of operating EDFs while maintaining the same weight to
% determine optimal cruise operation (with determined Antigravity T-motor)
for i = 1:numel(numEDFSweep)
    T_climb_ES = W * T_W_ratio_climb;
    T_hoverStill_ES = W;
    hoverRPM_Vcs = interpolate(T_climb_ES/ numEDFSweep(i), T_rseqHVc, rpm_rseqHVc);
    hoverRPM_Stills = interpolate(T_hoverStill_ES / numEDFSweep(i), T_rseqHS, rpm_rseqHS); % T = W
    
    cruiseRPMs = interpolate(T_cruises(i) / numEDFSweep(i), T_rseqC, rpm_rseqC);
    
    [bestMotorIndex_ES, motor_E_cruiseS(i), motor_E_HSS(i), motor_E_HVcS(i)] = ...
         getMotorEffsSingle(cruiseRPMs, hoverRPM_Stills, hoverRPM_Vcs, ...
         motorData, 1, throttle_max, voltage);

    powerCruise_units = interpolate(cruiseRPMs, rpm_rseqC, P_rseqC);
    prop_E_cruises = interpolate(cruiseRPMs, rpm_rseqC, eff_rseqC);
    power_total_Cruises = powerCruise_units * numEDFSweep(i);
    power_true_total_Cruises(i) = power_total_Cruises / (motor_E_cruiseS(i));

    powerActual_unit_hoverVcs = interpolate(hoverRPM_Vcs, rpm_rseqHVc, P_rseqHVc);
    powerActual_unit_hoverStills = interpolate(hoverRPM_Stills, rpm_rseqHS, P_rseqHS);
    
    power_total_hoverVcs = powerActual_unit_hoverVcs * numEDFSweep(i);
    power_total_hoverStills = powerActual_unit_hoverStills * numEDFSweep(i);

    power_true_total_hoverVcs(i) = power_total_hoverVcs / motor_E_HVcS(i);
    power_true_total_hoverStills(i) = power_total_hoverStills / motor_E_HSS(i);
end

figure;
% Subplot 1: Power Required 
subplot(1,2,1);
bar(numEDFSweep, power_true_total_Cruises, 'b');
title('Power Required vs. Number of Operating EDFs');
xticks(numEDFSweep);
xlabel('Number of Operating EDFs');
ylabel('Power Required (W)');
grid on;

% Subplot 2: Motor Efficiency 
subplot(1,2,2);
plot(numEDFSweep, motor_E_cruiseS, 'ro-');
title('Motor Efficiency vs. Number of Operating EDFs');
xlabel('Number of Operating EDFs');
ylabel('Motor Efficiency');
ylim([0 1]); % since efficiency is between 0 and 1
grid on;
sgtitle('Sweep of Number of Operating EDFs in a 6 EDF System - Cruise Condition');

for i = 1:numel(numEDFSweep)
    W_edfSweep = (W - (totalEDFWeightEst*g))  + (edf_mass_unit * g * numEDFSweep(i));
    T_climb_ES = W_edfSweep * T_W_ratio_climb;
    T_hoverStill_ES = W_edfSweep;
    hoverRPM_Vcs = interpolate(T_climb_ES/ numEDFSweep(i), T_rseqHVc, rpm_rseqHVc);
    hoverRPM_Stills = interpolate(T_hoverStill_ES / numEDFSweep(i), T_rseqHS, rpm_rseqHS); % T = W
    
    cruiseRPMs = interpolate(T_cruises(i) / numEDFSweep(i), T_rseqC, rpm_rseqC);
    
    [bestMotorIndex_ES, motor_E_cruiseS(i), motor_E_HSS(i), motor_E_HVcS(i)] = ...
        getMotorEffs(cruiseRPMs, hoverRPM_Stills, hoverRPM_Vcs, motorData, throttle_max, voltage);

    powerCruise_units = interpolate(cruiseRPMs, rpm_rseqC, P_rseqC);
    prop_E_cruises = interpolate(cruiseRPMs, rpm_rseqC, eff_rseqC);
    power_total_Cruises = powerCruise_units * numEDFSweep(i);
    power_true_total_Cruises(i) = power_total_Cruises / (motor_E_cruiseS(i));

    powerActual_unit_hoverVcs = interpolate(hoverRPM_Vcs, rpm_rseqHVc, P_rseqHVc);
    powerActual_unit_hoverStills = interpolate(hoverRPM_Stills, rpm_rseqHS, P_rseqHS);
    
    power_total_hoverVcs = powerActual_unit_hoverVcs * numEDFSweep(i);
    power_total_hoverStills = powerActual_unit_hoverStills * numEDFSweep(i);

    power_true_total_hoverVcs(i) = power_total_hoverVcs / motor_E_HVcS(i);
    power_true_total_hoverStills(i) = power_total_hoverStills / motor_E_HSS(i);
end
figure;
plot(numEDFSweep, power_true_total_hoverVcs, 'LineWidth', 2);
hold on;
plot(numEDFSweep, power_true_total_hoverStills, 'LineWidth', 2);
hold on;
plot(numEDFSweep, power_true_total_Cruises, 'LineWidth', 2);
hold on;
legend('Vertical Climb at Speed V_c', 'Pure Hover', 'Cruise');
title('Total Power Required vs. Number of EDFs on System')
xlabel('Number of EDFs');
ylabel('Total Power Required (W)');
grid on;


filename = 'ThesisWingOptimizationVariables.xlsx';
sheetname = 'PropellerOptimizationVariables';

% Read Table 
tbl = readtable(filename, 'Sheet', sheetname, 'ReadVariableNames', false);
cruise_rpm = tbl{strcmp(tbl{:,1}, 'RPM_Cruise'), 2:end};
cruise_power_unit = tbl{strcmp(tbl{:,1}, 'P_Cruise (Unit) (W)'), 2:end};
cruise_efficiency = tbl{strcmp(tbl{:,1}, 'Prop Efficiency for Cruise'), 2:end};
rpm_hover = tbl{strcmp(tbl{:,1}, 'RPM_Hover Still'), 2:end};
power_hover_unit = tbl{strcmp(tbl{:,1}, 'P_Hover Still (Unit) (W)'), 2:end};
fm_hover = tbl{strcmp(tbl{:,1}, 'FM for Hover Still'), 2:end};
rpm_climb = tbl{strcmp(tbl{:,1}, 'RPM_Hover at Vc'), 2:end};
power_climb_unit = tbl{strcmp(tbl{:,1}, 'P_Hover at Vc (Unit) (W)'), 2:end};
fm_climb = tbl{strcmp(tbl{:,1}, 'FM for Hover at Vc'), 2:end};


[sorted_rpm_cruise, idx_cruise] = sort(cruise_rpm);
sorted_power_cruise_unit = cruise_power_unit(idx_cruise);
sorted_eff_cruise = cruise_efficiency(idx_cruise);
[sorted_rpm_hover, idx_hover] = sort(rpm_hover);
sorted_power_hover_unit = power_hover_unit(idx_hover);
sorted_fm_hover = fm_hover(idx_hover);
[sorted_rpm_climb, idx_climb] = sort(rpm_climb);
sorted_power_climb_unit = power_climb_unit(idx_climb);
sorted_fm_climb = fm_climb(idx_climb);

figure;
% Climb Plot 
subplot(1,3,1);
yyaxis left;
plot(sorted_rpm_climb, sorted_power_climb_unit, 'b-o', 'LineWidth', 1.5);
ylabel('Power Required for Unit EDF (W)');
xlabel('RPM');
title('Climb');
grid on;
yyaxis right;
plot(sorted_rpm_climb, sorted_fm_climb, 'r--s', 'LineWidth', 1.5);
ylabel('Figure of Merit (FM)');
legend('Power', 'FM', 'Location', 'best');

% Hover Plot 
subplot(1,3,2);
yyaxis left;
plot(sorted_rpm_hover, sorted_power_hover_unit, 'b-o', 'LineWidth', 1.5);
ylabel('Power Required for Unit EDF (W)');
xlabel('RPM');
title('Pure Hover');
grid on;
yyaxis right;
plot(sorted_rpm_hover, sorted_fm_hover, 'r--s', 'LineWidth', 1.5);
ylabel('Figure of Merit (FM)');
legend('Power', 'FM', 'Location', 'best');

% Cruise Plot 
subplot(1,3,3);
yyaxis left;
plot(sorted_rpm_cruise, sorted_power_cruise_unit, 'b-o', 'LineWidth', 1.5);
ylabel('Power Required for Unit EDF (W)');
xlabel('RPM');
title('Cruise');
grid on;
yyaxis right;
plot(sorted_rpm_cruise, sorted_eff_cruise, 'm--s', 'LineWidth', 1.5);
ylabel('Propeller Efficiency');
legend('Power', 'Efficiency', 'Location', 'best');

sgtitle('RPM vs Power and Efficiency/FM for Flight Conditions');

% Total power = sum of all three
total_power_ES = power_true_total_hoverVcs + power_true_total_hoverStills + power_true_total_Cruises;

% Find index of minimum total power
[~, minIndex] = min(total_power_ES);
optimal_num_edfs = numEDFSweep(minIndex);

opt_power_true_total_Vc_ES = power_true_total_hoverVcs(minIndex);
opt_power_true_total_Still_ES = power_true_total_hoverStills(minIndex);
opt_power_true_total_Cruise_ES = power_true_total_Cruises(minIndex);
opt_motor_E_cruise = motor_E_cruiseS(minIndex);
opt_motor_E_Vc = motor_E_HVcS(minIndex);
opt_motor_E_Still = motor_E_HSS(minIndex);

varTableProp{strcmp(varTableProp.VariableName, "Selected Wing"), iterationColProp} = selectedWing;
varTableProp{strcmp(varTableProp.VariableName, "T for Hover Still (N)"), iterationColProp} = T_hover_min_unit;
varTableProp{strcmp(varTableProp.VariableName, "T_Hover (Unit) at V_c"), iterationColProp} = T_climb_unit;
varTableProp{strcmp(varTableProp.VariableName, "T_Cruise (Unit)"), iterationColProp} = T_cruise_unit;
varTableProp{strcmp(varTableProp.VariableName, "T_Glide (Unit)"), iterationColProp} = T_glide_unit;
varTableProp{strcmp(varTableProp.VariableName, "Number of EDFs"), iterationColProp} = numEDF;
varTableProp{strcmp(varTableProp.VariableName, "V_c"), iterationColProp} = V_climb_assessed;
varTableProp{strcmp(varTableProp.VariableName, "v_i at V_c"), iterationColProp} = v_i_hoverVc;
varTableProp{strcmp(varTableProp.VariableName, "RPM_Hover at Vc"), iterationColProp} = hoverRPM_Vc;
varTableProp{strcmp(varTableProp.VariableName, "P_Hover at Vc (Unit) (W)"), iterationColProp} = powerActual_unit_hoverVc;
varTableProp{strcmp(varTableProp.VariableName, "FM for Hover at Vc"), iterationColProp} = FM_hoverVc;
varTableProp{strcmp(varTableProp.VariableName, "FM for Hover at Vc"), iterationColProp} = FM_hoverVc;
varTableProp{strcmp(varTableProp.VariableName, "P_Hover (Total) at Vc (W)"), iterationColProp} = power_true_total_hoverVc;
varTableProp{strcmp(varTableProp.VariableName, "v_i for Hover Still"), iterationColProp} = v_i_hoverStill;
varTableProp{strcmp(varTableProp.VariableName, "RPM_Hover Still"), iterationColProp} = hoverRPM_Still;
varTableProp{strcmp(varTableProp.VariableName, "P_Hover Still (Unit) (W)"), iterationColProp} = powerActual_unit_hoverStill;
varTableProp{strcmp(varTableProp.VariableName, "FM for Hover Still"), iterationColProp} = FM_hoverStill;
varTableProp{strcmp(varTableProp.VariableName, "P_Hover (Total) Still (W)"), iterationColProp} = power_true_total_hoverStill;
varTableProp{strcmp(varTableProp.VariableName, "RPM_Cruise"), iterationColProp} = cruiseRPM;
varTableProp{strcmp(varTableProp.VariableName, "P_Cruise (Unit) (W)"), iterationColProp} = power_true_unit_cruise;
varTableProp{strcmp(varTableProp.VariableName, "Prop Efficiency for Cruise"), iterationColProp} = prop_E_cruise;
varTableProp{strcmp(varTableProp.VariableName, "P_Cruise (Total) (W)"), iterationColProp} = power_true_total_Cruise;

varTableProp{strcmp(varTableProp.VariableName, "Selected Motor Name"), iterationColProp} = bestMotorIndex;
varTableProp{strcmp(varTableProp.VariableName, "Selected Motor_E Cruise"), iterationColProp} = motor_E_curr_cruise;
varTableProp{strcmp(varTableProp.VariableName, "Selected Motor_E Hover Vc"), iterationColProp} = motor_E_curr_hoverVc;
varTableProp{strcmp(varTableProp.VariableName, "Selected Motor_E Hover Still"), iterationColProp} = motor_E_curr_hoverStill;
varTableProp{strcmp(varTableProp.VariableName, "Optimal Number of EDFs"), iterationColProp} = optimal_num_edfs;
varTableProp{strcmp(varTableProp.VariableName, "Optimal Motor Name"), iterationColProp} = bestMotorIndex_ES;
varTableProp{strcmp(varTableProp.VariableName, "Optimal Motor_E Cruise"), iterationColProp} = opt_motor_E_cruise;
varTableProp{strcmp(varTableProp.VariableName, "Optimal Motor_E Hover Vc"), iterationColProp} = opt_motor_E_Vc;
varTableProp{strcmp(varTableProp.VariableName, "Optimal Motor_E Hover Still"), iterationColProp} = opt_motor_E_Still;

varTableProp{strcmp(varTableProp.VariableName, "Optimal True P_Cruise (Total) (W)"), iterationColProp} = opt_power_true_total_Cruise_ES;
varTableProp{strcmp(varTableProp.VariableName, "Optimal True P_Hover (Total) at Vc (W)"), iterationColProp} = opt_power_true_total_Vc_ES;
varTableProp{strcmp(varTableProp.VariableName, "Optimal True P_Hover (Total) Still (W)"), iterationColProp} = opt_power_true_total_Still_ES;

writetable(varTableProp, "ThesisWingOptimizationVariables.xlsx", 'Sheet', sheetName, 'WriteRowNames', true);
disp(varTableProp);

propDesign = false;

%% Functions
function col2_output = interpolate(col1_target, col1Array, col2Array)
    if col1_target < min(col1Array) || col1_target > max(col1Array)
        warning(['Target value ', inputname(1), ' is out of range for interpolation!']);
        col2_output = min(col2Array);
    else
        % Linear interpolation to find corresponding col2_output
        col2_output = interp1(col1Array, col2Array, col1_target, "linear");
    end
end

function J = findJ(V_freestream, RPM, radius)
    omega = RPM * (2*pi/60); % rad/s
    J = V_freestream / (omega*radius);
end

function col2Target = findValuefromTargetWithFit(col1Target, col1Array, col2Array)
    p = polyfit(col2Array, col1Array, 2); 
    fine = linspace(min(col2Array), max(col2Array), 10000);
    fit = polyval(p, fine);
    [~, idx] = min(abs(fit - col1Target));
    col2Target = fine(idx);
end

function powerIdealUnit = idealHoverPower(thrust_unit_target, vc, rho, radius_m)
    A = pi * radius_m^2; % m^2 (area of ONE rotor)
    powerIdealUnit = thrust_unit_target * ...
        ((vc / 2) + sqrt((vc/2)^2 + (thrust_unit_target / (2 * rho * A))));
end

function eta = eta_func_rpm(rpm, volt, I0, Res, Kv)
    denom = volt - rpm / Kv;
    if denom <= 1e-3
        eta = 0;
        return;
    end
    
    eta = (1 - (I0 * Res) / denom) * (rpm / Kv) / (volt);  % Efficiency as ratio of mechanical power to electrical power
    
    % Clamp efficiency between 0 and 1
    eta = max(min(eta, 1), 0);
end

function [bestMotorIndex, motor_E_curr_cruise, motor_E_curr_hoverStill, motor_E_curr_hoverVc] = ...
    getMotorEffs(cruiseRPM, hoverRPM_Still, hoverRPM_Vc, motorData, throttle_max, voltage)

% Extract motor parameters
motor_names = motorData.MotorName;
Kv = motorData.K_v;
I0 = motorData.I_0_A_;
Res = motorData.R_ohms_;

motors = {};
for i = 1:height(motorData)
    motors = [motors; {motor_names{i}, Kv(i), I0(i), Res(i)}];
end

voltage_supplied = throttle_max * voltage;

% Use RPM 
target_rpms = [cruiseRPM, hoverRPM_Vc, hoverRPM_Still];
best_total_eta = -Inf;

for i = 1:height(motorData)
    Kv_val = motors{i, 2};
    I0_val = motors{i, 3};
    R_val = motors{i, 4};
    motor_name = motor_names{i};

    rpm_max = voltage_supplied * Kv_val;
    rpm_range = linspace(0, rpm_max * 1.2, 200);
    effs = zeros(size(rpm_range));

    for j = 1:length(rpm_range)
        rpm = rpm_range(j);  % Current RPM value
        effs(j) = eta_func_rpm(rpm, voltage_supplied, I0_val, R_val, Kv_val);  % Call the function with correct arguments
    end

    eta_at_targets = interp1(rpm_range, effs, target_rpms, 'linear', 0);
    total_eta = sum(eta_at_targets);

    if total_eta > best_total_eta
        bestMotorIndex = i;
        best_total_eta = total_eta;
        best_motor_name = motor_name;
        best_eta_values = eta_at_targets;
        motor_E_curr_cruise = best_eta_values(1);
        motor_E_curr_hoverVc = best_eta_values(2);
        motor_E_curr_hoverStill = best_eta_values(3);
    end
end

% Print results
fprintf("Best Motor: %s\n", best_motor_name);
fprintf("Efficiencies at [%d, %d, %d, %d] RPM: [%.4f, %.4f, %.4f, %.4f]\n", ...
    target_rpms(1), target_rpms(2), target_rpms(3), ...
    best_eta_values(1), best_eta_values(2), best_eta_values(3));
end

function [motorIndex, motor_E_curr_cruise, motor_E_curr_hoverStill, motor_E_curr_hoverVc] = ...
    getMotorEffsSingle(cruiseRPM, hoverRPM_Still, hoverRPM_Vc, motorData, motorIdx, throttle_max, voltage)

% Extract motor parameters
motor_names = motorData.MotorName;
Kv = motorData.K_v;
I0 = motorData.I_0_A_;
Res = motorData.R_ohms_;

% Get selected motor
Kv_val = Kv(motorIdx);
I0_val = I0(motorIdx);
R_val = Res(motorIdx);
motor_name = motor_names{motorIdx};

voltage_supplied = throttle_max * voltage;

% Define target RPMs
target_rpms = [cruiseRPM, hoverRPM_Vc, hoverRPM_Still];

% Simulate over a range to build efficiency curve
rpm_max = voltage_supplied * Kv_val;
rpm_range = linspace(0, rpm_max * 1.2, 200);
effs = zeros(size(rpm_range));

for j = 1:length(rpm_range)
    rpm = rpm_range(j);
    effs(j) = eta_func_rpm(rpm, voltage_supplied, I0_val, R_val, Kv_val);
end

% Interpolate efficiencies at target RPMs
eta_at_targets = interp1(rpm_range, effs, target_rpms, 'linear', 0);

% Output results
motorIndex = motorIdx;
motor_E_curr_cruise = eta_at_targets(1);
motor_E_curr_hoverVc = eta_at_targets(2);
motor_E_curr_hoverStill = eta_at_targets(3);

% Print results
fprintf("Selected Motor: %s\n", motor_name);
fprintf("Efficiencies at [%d, %d, %d] RPM: [%.4f, %.4f, %.4f]\n", ...
    target_rpms(1), target_rpms(2), target_rpms(3), ...
    motor_E_curr_cruise, motor_E_curr_hoverVc, motor_E_curr_hoverStill);

end
