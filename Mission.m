%% 4. Mission (Defines the mission - time intervals and distances)

% Calls Flight Envelope
FlightEnvelope;

selectedWing = iteration;

V_sel = varTable{strcmp(varTable.VariableName, "Exhaust Velocity (m/s)"), selectedWing}; 
S_sel = varTable{strcmp(varTable.VariableName, "Wing area (m^2)"), selectedWing}; 
MAC_wing_sel = varTable{strcmp(varTable.VariableName, "Wing MAC (m)"), selectedWing};
rho_sel = varTable{strcmp(varTable.VariableName, "rho (kg/m^3)"), selectedWing};
C_Dtot_level_sel = varTable{strcmp(varTable.VariableName, "Level C_D"), selectedWing};
T_Hover_min_sel = varTable{strcmp(varTable.VariableName, "'Min. Thrust for Hover (N)'"), selectedWing};
D_cruise_sel = varTable{strcmp(varTable.VariableName, "Level D (N)"), selectedWing};
D_glide_sel = varTable{strcmp(varTable.VariableName, "Glide D (N)"), selectedWing};
opt_glide_angle_deg_sel = varTable{strcmp(varTable.VariableName, "Optimal Glide Angle (deg)"), selectedWing};
opt_glide_velocity_sel = varTable{strcmp(varTable.VariableName, "Optimal Glide Velocity (m/s)"), selectedWing};
C_L_glide_sel = varTable{strcmp(varTable.VariableName, "Glide C_L"), selectedWing};
C_Dtot_glide_sel = varTable{strcmp(varTable.VariableName, "Glide C_D"), selectedWing};


optimal_glide_angle = opt_glide_angle_deg_sel;
optimal_glide_velocity = opt_glide_velocity_sel; % m/s

% Givens
mountainH = 500;
h_clear_tree = 50;
h_cruise = 500;
tot_hor_dist = 25000;

dist_start2clearhill = 17000;


glide1_hor_dist = h_clear_tree / tand(optimal_glide_angle);
glide2_hor_dist = h_cruise / tand(optimal_glide_angle);

cruise_hor_dist = tot_hor_dist - (glide1_hor_dist + glide2_hor_dist);

glide1_true_dist = sqrt(h_clear_tree^2 + glide1_hor_dist^2);
glide2_true_dist = sqrt(h_cruise^2 + glide2_hor_dist^2);

t_cruise = cruise_hor_dist / V_sel;
t_glide1 = glide1_true_dist / optimal_glide_velocity;
t_glide2 = glide2_true_dist / optimal_glide_velocity;

disp(['Total minutes of cruise: ', num2str(t_cruise / 60)]);
disp(['Distance in level flight: ', num2str(cruise_hor_dist)]);