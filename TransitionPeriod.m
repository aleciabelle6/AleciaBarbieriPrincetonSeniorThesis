fprintf(['The fitted polynomial function is: C_D = ' ...
    '%.4f + %.4f*C_L + %.4f*C_L^2 + %.4f*C_L^3\n'], ...
    pCD_CL(4), pCD_CL(3), pCD_CL(2), pCD_CL(1));

% Inputs
rho = rho_curr; % Air density at cruise
n = 3; % Load factor during pull-up
CD0 = C_d0_level_WH; % Parasitic drag coefficient at cruise
eta_prop = prop_E_cruise;  % Propeller efficiency at cruise

% Velocities
V_glide = velo_glide_optimal;
V_cruise = V;
V_avg = (V_glide + V_cruise) / 2;  % Average velocity during pull-up

% Step 1: Compute CL at increased lift using average velocity
CL = (2 * n * W) / (rho * V_avg^2 * S);

% Step 2: Total drag coefficient
CD = AeroCoeffs.YfromAlpha(C_L{1}, C_D{1}, CL);

% Step 3: Drag force at V_avg
D = 0.5 * rho * V_avg^2 * S * CD;

% Step 4: Power required to overcome drag at V_avg
P_peak = (D * V_avg) / eta_prop;   % [W]

% Display result
fprintf('Peak power required during pull-up: %.2f kW\n', P_peak / 1000);

% Compute minimum pull-up radius to stay under n_max g
R_min = V_avg^2 / (g * (n - 1)); % Uses V_avg now
theta_deg = glide_angle_degrees; % Entry angle into the turn in degrees
theta_rad = deg2rad(theta_deg);

% Compute transition time using arc length and average velocity
t_transition = (R_min * theta_rad) / V_avg;

% Time step and time vector
dt = 0.01;
t = 0:dt:t_transition;
n_points = length(t);
half_idx = round(n_points / 2);

% === Power Profile: 0 W → Peak → Cruise ===
P_profile = zeros(size(t));

% Phase 1: Ramp from 0 to peak power
P_profile(1:half_idx) = linspace(0, P_peak, half_idx);

% Phase 2: Ramp from peak to cruise power
P_profile(half_idx+1:end) = linspace(P_peak, power_true_total_Cruise, n_points - half_idx);

% === Cumulative Energy Consumption (in Wh) ===
E_wh_G2C_Array = cumtrapz(t, P_profile) / 3600;  % Convert from Ws to Wh
E_wh_G2C = E_wh_G2C_Array(end);

% === Plotting ===

% Plot 1: Power Profile
figure;
subplot(2,1,1)
plot(t, P_profile, 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Power (W)');
title('Power Profile: Pull-Up Transition from Unpowered Glide to Cruise');
grid on;

% Plot 2: Cumulative Energy
subplot(2,1,2)
plot(t, E_wh_G2C_Array, 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Cumulative Energy (Wh)');
title('Energy Consumption vs. Time');
grid on;

% === Display Total Energy Consumed ===
fprintf('Total energy consumed during pull-up: %.4f Wh\n', E_wh_G2C_Array(end));
fprintf('Total energy consumed during pull-up: %.4f J\n', E_wh_G2C_Array(end)*3600);


% === Inputs ===
P_cruise = power_true_total_Cruise;  % Known cruise power [W]
V_glide = velo_glide_optimal;
V_cruise = V;
V_avg = (V_glide + V_cruise) / 2;     % Average velocity during transition
theta_deg = glide_angle_degrees;     % Same flight path angle in degrees
theta_rad = deg2rad(theta_deg);

% === Use same pull-up radius, assuming symmetric arc
R_min = V_avg^2 / (g * (n - 1));     % Based on V_avg and same n as before

% === Compute transition time (arc length / avg velocity)
t_transition = (R_min * theta_rad) / V_avg;

% === Time step and vector
dt = 0.01;
t = 0:dt:t_transition;
n_points = length(t);

% === Power profile: Cruise → 0 W
P_profile = linspace(P_cruise, 0, n_points);  % Linear ramp down

% === Cumulative energy in Wh
E_wh_C2G_Array = cumtrapz(t, P_profile) / 3600;
E_wh_C2G = E_wh_C2G_Array(end);

% === Plotting ===
figure;
subplot(2,1,1)
plot(t, P_profile, 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Power (W)');
title('Power Profile: Cruise to Unpowered Glide');
grid on;

subplot(2,1,2)
plot(t, E_wh_C2G_Array, 'r-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Cumulative Energy (Wh)');
title('Energy Consumption vs. Time');
grid on;

% === Output energy
fprintf('Total energy consumed during cruise-to-glide: %.4f Wh\n', E_wh_C2G_Array(end));
fprintf('Total energy consumed during cruise-to-glide: %.4f J\n', E_wh_C2G_Array(end)*3600);
