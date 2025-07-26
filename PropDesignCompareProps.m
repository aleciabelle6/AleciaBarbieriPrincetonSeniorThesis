%% 6. PropDesign (Sizing and Optimizing Propeller for Mission Requirements (+ Hover))

% Call 5. Propulsion
Propulsion;
propIterationRange = 15:22;
for m = 1:numel(propIterationRange)
    propIterationStr = string(propIterationRange(m));
    propIterationVal = propIterationRange(m);  % numeric
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
prefix = 'XROTORSweepData/';  % Folder name for polar data
propSweeps = cellfun(@(x) [prefix, x], propSweepsTemp, 'UniformOutput', false);
    % Create table and write to Excel
    varTableProp = readtable("ThesisWingOptimizationVariables.xlsx", 'Sheet', 'PropellerOptimizationVariables', 'ReadVariableNames', true, 'VariableNamingRule', 'preserve');
iterationColProp = find(strcmp(varTableProp.Properties.VariableNames, propIterationStr));

    
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
    
    %% XRotor Analaysis
radius_sim = varTableProp{strcmp(varTableProp.VariableName, "radius(m)"), propIterationVal}; 
currPropIteration = propIterationVal + 1;
filenameSetProp = propSweeps(currPropIteration, :); 

    colMap = cell(1, numel(filenameSetProp));
    dataAll = cell(1, numel(filenameSetProp));  % store data separately for each file
    
    % Loop through all files
    for i = 1:numel(filenameSetProp)
        fid = fopen(filenameSetProp{i}, 'r');
    
        % Skip until header line (e.g., starts with 'n' and includes 'eff')
        header = '';
        while ~feof(fid)
            line = fgetl(fid);
            if contains(line, 'n') && contains(line, 'eff')
                header = strtrim(line);
                break;
            end
        end
    
        % Skip the dashed separator line
        while ~feof(fid)
            line = fgetl(fid);
            if contains(line, '---')
                break;
            end
        end
    
        % Extract column names and create mapping
        colNames = strsplit(strtrim(header));
        colMap{i} = containers.Map(colNames, 1:length(colNames));
    
        % Read and store data
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
            data = [data; rowVals'];  % Store as row
        end
    
        fclose(fid);
        dataAll{i} = data;
    end
    
    % Extract from first file (rseqH)
    colMap_rseqHVc = colMap{1};
    colMap_rseqHS   = colMap{2};
    colMap_rseqC = colMap{3};

    data_rseqHVc   = dataAll{1};
    data_rseqHS   = dataAll{2};
    data_rseqC   = dataAll{3};
    
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
    
    rhos = [rho_rseqC(1), rho_rseqHVc(1), rho_rseqHS(1)];
    
    if ~all(rhos == rhos(1)) || ~all(rho_rseqC == rho_rseqC(1)) || ~all(rho_rseqHVc == rho_rseqHVc(1)) || ~all(rho_rseqHS == rho_rseqHS(1)) 
        error('Rho must be consistent!');
    end


    powerCruise_m = zeros(1, numel(numEDFSweep));
powerHoverVcs_m = zeros(1, numel(numEDFSweep));
powerHoverStills_m = zeros(1, numel(numEDFSweep));

    
    for i = 1:numel(numEDFSweep)

% System parameters
voltage = 43.2;

% peak power is 13.29kW
peak_power = 13400;
peak_power_unit = peak_power / numEDF;
I_max_motor = 52;
P_max_motor = I_max_motor*voltage;
throttle_max = peak_power_unit / (P_max_motor);
voltage_supplied = voltage * throttle_max;

% Read motor data
motorData = readtable('ThesisWingOptimizationVariables.xlsx', 'Sheet', 'MotorSelection');
motor_names = motorData.MotorName;
Kv = motorData.K_v;
I0 = motorData.I_0_A_;
Res = motorData.R_ohms_;

        W_edfSweep = (W - (totalEDFWeightEst*g))  + (edf_mass_unit * g * numEDFSweep(i));
        T_climb_ES = W_edfSweep * T_W_ratio_climb;
        T_hoverStill_ES = W_edfSweep;
        hoverRPM_Vcs = interpolate(T_climb_ES/ numEDFSweep(i), T_rseqHVc, rpm_rseqHVc);
        hoverRPM_Stills = interpolate(T_hoverStill_ES / numEDFSweep(i), T_rseqHS, rpm_rseqHS); % T = W
        
        cruiseRPMs = interpolate(T_cruises(i) / numEDFSweep(i), T_rseqC, rpm_rseqC);
        
        [bestMotorIndex_ES(i,m), motor_E_cruiseS(i), motor_E_HSS(i), motor_E_HVcS(i)] = ...
            getMotorEffs(cruiseRPMs, hoverRPM_Stills, hoverRPM_Vcs, motorData, throttle_max, voltage);
    
        powerCruise_units = interpolate(cruiseRPMs, rpm_rseqC, P_rseqC);
        prop_E_cruises = interpolate(cruiseRPMs, rpm_rseqC, eff_rseqC);
        power_total_Cruises = powerCruise_units * numEDFSweep(i);



        powerActual_unit_hoverVcs = interpolate(hoverRPM_Vcs, rpm_rseqHVc, P_rseqHVc);
        powerActual_unit_hoverStills = interpolate(hoverRPM_Stills, rpm_rseqHS, P_rseqHS);
        
        power_total_hoverVcs = powerActual_unit_hoverVcs * numEDFSweep(i);
        power_total_hoverStills = powerActual_unit_hoverStills * numEDFSweep(i);
    
powerCruise_m(i) = power_total_Cruises / motor_E_cruiseS(i);
powerHoverVcs_m(i) = power_total_hoverVcs / motor_E_HVcS(i);
powerHoverStills_m(i) = power_total_hoverStills / motor_E_HSS(i);

power_true_total_Cruises{m} = powerCruise_m;
power_true_total_hoverVcs{m} = powerHoverVcs_m;
power_true_total_hoverStills{m} = powerHoverStills_m;

    end
end
  
figure;

% Colormap and legend labels
cmap = lines(numel(propIterationRange));  % Generate distinct colors
legendEntries = arrayfun(@(x) ['Propeller ', num2str(x)], propIterationRange, 'UniformOutput', false);

% Store plot handles for legend
cruiseHandles = gobjects(1, numel(propIterationRange));


% === Cruise Subplot ===
subplot(3,1,1); hold on;
for m = 1:numel(propIterationRange)
    cruiseHandles(m) = plot(numEDFSweep, power_true_total_Cruises{m}, ...
        'Color', cmap(m,:), 'LineWidth', 1.5);
end
title('Cruise');
xticks(floor(min(numEDFSweep)):ceil(max(numEDFSweep)));
ylabel('Total Power (W)');
grid on;

% === Climb Subplot (Hover with V_c) ===
subplot(3,1,2); hold on;
for m = 1:numel(propIterationRange)
    plot(numEDFSweep, power_true_total_hoverVcs{m}, ...
        'Color', cmap(m,:), 'LineWidth', 1.5);
end
title('Climb');
xticks(floor(min(numEDFSweep)):ceil(max(numEDFSweep)));
ylabel('Total Power (W)');
grid on;

% === Hover Subplot (Still Hover) ===
subplot(3,1,3); hold on;
for m = 1:numel(propIterationRange)
    plot(numEDFSweep, power_true_total_hoverStills{m}, ...
        'Color', cmap(m,:), 'LineWidth', 1.5);
end
title('Pure Hover');
xticks(floor(min(numEDFSweep)):ceil(max(numEDFSweep)));
xlabel('Number of EDFs'); ylabel('Total Power (W)');
grid on;

% === Unified Legend (based on cruiseHandles) ===
legend(cruiseHandles, legendEntries, 'Location', 'eastoutside');

% Title for full figure
sgtitle('Power Required vs EDF Count Across Best Propeller Designs');

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

% === Use RPM ===
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
