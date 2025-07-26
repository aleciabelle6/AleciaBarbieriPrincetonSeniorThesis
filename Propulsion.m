%% 5. Propulsion (Determines Thrust Required for each Segment: Cruise, Optimal Glide, Hover)

% Calls 4. Mission
Mission;

g = 9.81; % m/s^2

% Cruise
T_cruise = D_cruise_sel;

% Glide
T_glide = abs(D_glide_sel - W*sind(opt_glide_angle_deg_sel));

% Hover
T_W_ratio_climb = 1.015;
T_climb = W * T_W_ratio_climb; % for incrementing when T/W > 1
T_Hover_min = W; % for still hover after being in the air

thrustsSegmented = [T_climb, T_Hover_min, T_cruise, T_glide];
disp(thrustsSegmented);

maxThrust = max(thrustsSegmented);
T_mission_max_tot = max(thrustsSegmented); % N

T_climb_unit = T_climb / numEDF;
T_cruise_unit = T_cruise / numEDF;
T_glide_unit = T_glide / numEDF;

T_hover_min_unit = T_Hover_min / numEDF;
thrustsSegmented_unit = thrustsSegmented ./ numEDF;
T_mission_max_unit  = max(thrustsSegmented_unit); % N