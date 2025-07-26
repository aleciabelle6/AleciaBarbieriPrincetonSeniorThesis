%% Motor Efficiency and Torque Calculations (RPM units) 

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
motorData = readtable('ThesisWingOptimizationVariables.xlsx', 'Sheet', 'MotorSelection Final');
motor_names = motorData.MotorName;
Kv = motorData.K_v;
I0 = motorData.I_0_A_;
Res = motorData.R_ohms_;

motors = {};
for i = 1:height(motorData)
    motors = [motors; {motor_names{i}, Kv(i), I0(i), Res(i)}];
end

% Target RPMs
target_rpms = [cruiseRPM, hoverRPM_Vc, hoverRPM_Still];
num_motors = height(motorData);
best_total_eta = -Inf;

% Plot Efficiency vs RPM 
colors = lines(num_motors);
figure('Name', 'Motor Efficiency vs RPM'); hold on;

for i = 1:num_motors
    Kv_val = motors{i,2};
    I0_val = motors{i,3};
    R_val = motors{i,4};
    motor_name = motor_names{i};

    rpm_max = voltage_supplied * Kv_val;
    rpm_range = linspace(0, rpm_max * 1.2, 200);
    effs = zeros(size(rpm_range));

    for j = 1:length(rpm_range)
        rpm = rpm_range(j);
        effs(j) = eta_func_rpm(rpm, voltage_supplied, I0_val, R_val, Kv_val);
    end

    eta_at_targets = interp1(rpm_range, effs, target_rpms, 'linear', 0);
    total_eta = sum(eta_at_targets);
    
    if total_eta > best_total_eta
        bestMotorIndex = i;
        best_total_eta = total_eta;
        best_motor_name_curr = motor_name;
        best_eta_values = eta_at_targets;
        motor_E_curr_cruise = best_eta_values(1);
        motor_E_curr_hoverVc = best_eta_values(2);
        motor_E_curr_hoverStill = best_eta_values(3);
    end

    plot(rpm_range, effs, 'Color', colors(i, :), 'DisplayName', motor_name, 'LineWidth', 1);
end

xlabel('RPM'); ylabel('Efficiency');
title('Motor Efficiency vs RPM'); legend show; grid on;

fprintf("Best Motor: %s\n", best_motor_name_curr);
fprintf("Efficiencies at target RPMs: [%.4f, %.4f, %.4f, %.4f, %.4f]\n", best_eta_values);

%% Motor Torque vs RPM 
figure('Name', 'Motor Torque vs RPM'); hold on;


% Define motor parameters as a cell array of structs for better readability
motors = {
    struct('name', 'T-Motor MN8014 Antigravity Type Kv100', 'Kv', 100, 'I0', 1.4, 'R', 0.065),
    struct('name', '100KV Brushless Heavy Lift Drone Motor', 'Kv', 100, 'I0', 0.8, 'R', 0.07)
};

voltage_supplied = 43.2; % [V]
num_motors = length(motors);

figure;
hold on;

for i = 1:num_motors
    Kv_val = motors{i}.Kv;
    I0_val = motors{i}.I0;
    R_val = motors{i}.R;
    motor_name = motors{i}.name;

    rpm_range = linspace(0, voltage_supplied * Kv_val * 1.2, 200);
    torques = zeros(size(rpm_range));

    for j = 1:length(rpm_range)
        rpm = rpm_range(j);
        volt_backemf = rpm / Kv_val;
        denom = voltage_supplied - volt_backemf;

        if denom <= 1e-3
            torques(j) = 0;
            continue;
        end

        I = denom / R_val;
        Kt = 60 / (2 * pi * Kv_val);  % [Nm/A]
        tq = Kt * (I - I0_val);
        torques(j) = max(tq, 0);
    end

    plot(rpm_range, torques, 'DisplayName', motor_name);
end
hold on;
plot(cruiseRPM, Q_cruise, 'ro', 'MarkerSize', 4, 'LineWidth', 1, ...
    'MarkerFaceColor', 'r', 'DisplayName', 'Cruise Operation');
hold on;

plot(hoverRPM_Still, Q_HS, 'bo', 'MarkerSize', 4, 'LineWidth', 1, ...
    'MarkerFaceColor', 'b', 'DisplayName', 'Pure Hover Operation');
hold on;

plot(hoverRPM_Vc, Q_Vc, 'mo', 'MarkerSize', 4, 'LineWidth', 1, ...
    'MarkerFaceColor', 'm', 'DisplayName', 'Vertical Climb Operation');
xlim([0, 5000]);
xlabel('RPM');
ylabel('Torque (N-m)');
title('Motor Torque vs RPM');
legend show;
grid on;

%% Efficiency Function (RPM-based) 
function eta = eta_func_rpm(rpm, volt, I0, Res, Kv)
    V_backemf = rpm / Kv;
    denom = volt - V_backemf;
    if denom <= 1e-3
        eta = 0;
        return;
    end
    eta = (1 - I0 * Res / denom) * rpm / (volt * Kv);
    eta = max(min(eta, 1), 0);
end
