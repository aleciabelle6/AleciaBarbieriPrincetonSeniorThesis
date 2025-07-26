%% Aero (Initial Mass Estimates and Wing Geometry Optimization)
clear; close all;

skip = [0,1,18,19,21,22,28,33];
%for iterationTemp = setdiff(2:42, skip)
    %iteration = string(iterationTemp);
    iteration = "42";
    cfShellThickness = 0.0007; % m
    CL_max = 1.7;
    AoA_stall = 18; % 15-18 deg, but 18 gives lowest W/S (safest benchmark for max W/S)
    
    polarValuesTemp = {'polar_i_WH_P.txt', 'polar_i_WP_V.txt';
                   'polar_1_WH_P.txt', 'polar_1_WP_V.txt';
                   'polar_2_WH_P.txt', 'polar_2_WP_V.txt';                  
                   'polar_3_WH_P.txt', 'polar_3_WP_V.txt'; 
                   'polar_4_WH_P.txt', 'polar_4_WP_V.txt';
                   'polar_5_WH_P.txt', 'polar_5_WP_V.txt';
                   'polar_6_WH_P.txt', 'polar_6_WP_V.txt';
                   'polar_7_WH_P.txt', 'polar_7_WP_V.txt';
                   'polar_8_WH_P.txt', 'polar_8_WP_V.txt';
                   'polar_9_WH_P.txt', 'polar_9_WP_V.txt';
                   'polar_10_WH_P.txt', 'polar_10_WP_V.txt';
                   'polar_11_WH_P.txt', 'polar_11_WP_V.txt';
                   'polar_12_WH_P.txt', 'polar_12_WP_V.txt';
                   'polar_13_WH_P.txt', 'polar_13_WP_V.txt';
                   'polar_14_WH_P.txt', 'polar_14_WP_V.txt';
                   'polar_15_WH_P.txt', 'polar_15_WP_V.txt';
                   'polar_16_WH_P.txt', 'polar_16_WP_V.txt';
                   'polar_17_WH_P.txt', 'polar_17_WP_V.txt';
                   'polar_18_WH_P.txt', 'polar_18_WP_V.txt';
                   'polar_19_WH_P.txt', 'polar_19_WP_V.txt';
                   'polar_20_WH_P.txt', 'polar_20_WP_V.txt';
                   'polar_21_WH_P.txt', 'polar_21_WP_V.txt';
                   'polar_22_WH_P.txt', 'polar_22_WP_V.txt';
                   'polar_23_WH_P.txt', 'polar_23_WP_V.txt';
                   'polar_24_WH_P.txt', 'polar_24_WP_V.txt';
                   'polar_25_WH_P.txt', 'polar_25_WP_V.txt';
                   'polar_26_WH_P.txt', 'polar_26_WP_V.txt';
                   'polar_27_WH_P.txt', 'polar_27_WP_V.txt';
                   'polar_28_WH_P.txt', 'polar_28_WP_V.txt';
                   'polar_29_WH_P.txt', 'polar_29_WP_V.txt';
                   'polar_30_WH_P.txt', 'polar_30_WP_V.txt';
                   'polar_31_WH_P.txt', 'polar_31_WP_V.txt';
                   'polar_32_WH_P.txt', 'polar_32_WP_V.txt';
                   'polar_33_WH_P.txt', 'polar_33_WP_V.txt';
                   'polar_34_WH_P.txt', 'polar_34_WP_V.txt';
                   'polar_35_WH_P.txt', 'polar_35_WP_V.txt';
                   'polar_36_WH_P.txt', 'polar_36_WP_V.txt';
                   'polar_37_WH_P.txt', 'polar_37_WP_V.txt';
                   'polar_38_WH_P.txt', 'polar_38_WP_V.txt';
                   'polar_39_WH_P.txt', 'polar_39_WP_V.txt';
                   'polar_40_WH_P.txt', 'polar_40_WP_V.txt';
                   'polar_41_WH_P.txt', 'polar_41_WP_V.txt';
                   'polar_42_WH_P.txt', 'polar_42_WP_V.txt';
                   'polar_43_WH_P.txt', 'polar_43_WP_V.txt';
                   };
    prefix = 'PolarData/'; 
    polarValues = cellfun(@(x) [prefix, x], polarValuesTemp, 'UniformOutput', false);
    
    varTable = readtable("ThesisWingOptimizationVariables.xlsx", 'ReadVariableNames', true, 'VariableNamingRule', 'preserve');
    AR_proj = varTable{strcmp(varTable.VariableName, "AR"), iteration}; 
    wingspan_true = varTable{strcmp(varTable.VariableName, "Wingspan (m)"), iteration}; 
    wingspan_proj = varTable{strcmp(varTable.VariableName, "Wingspan_proj"), iteration}; 
    S = varTable{strcmp(varTable.VariableName, "Wing area (m^2)"), iteration}; 
    MAC_wing = varTable{strcmp(varTable.VariableName, "Wing MAC (m)"), iteration}; 
    root_c = varTable{strcmp(varTable.VariableName, "Wing Section Root C (m)"), iteration}; 
    wing_sweep_deg = varTable{strcmp(varTable.VariableName, "Wing Sweep (deg)"), iteration}; 
    wingWeight = varTable{strcmp(varTable.VariableName, "Wing Mass (kg)"), iteration};
    
    h_curr = varTable{strcmp(varTable.VariableName, "h (m)"), iteration}; % m
    V = varTable{strcmp(varTable.VariableName, "Exhaust Velocity (m/s)"), iteration}; % m/s
    
    AR_true = wingspan_true^2 / S;
    
    g = 9.81; % m/s^2
    
    % Weights in pounds
    myWeightLB = 120; %lb
    %batteryWeightLB = 5; %lb
    
    % Weights in kg
    myWeight = myWeightLB * 0.453592; %kg
    %batteryWeightEst = batteryWeightLB * 0.453592; %kg
    batteryWeightEst = 8.34; %kg
    
    carbonFiberDensity = 1400; % kg/m^3
    
    % EDF Calculations
    numEDF = 6;
    num_blades = 12;
    radius = 0.3; % m
    hub_radius = 0.05;       % m
    duct_thickness_tot = 0.003;  % m
    blade_thickness = 0.002; % m
    blade_chord = 0.0175;      % m

    hub_length = 0.12;  % m
    rho_kevlar = 1440;  % kg/m^3
    motor_mass_unit = 0.315; %kg
    r_motor = 0.0439; % m
    h_motor = 0.0315; % m
    rho_rohacell_foam = 32; % kg/m^3
    
    blade_area = blade_chord * (radius - hub_radius);  % m^2
    blade_volume = blade_area * blade_thickness;       % m^3
    blade_mass = blade_volume * rho_kevlar;           
    total_blades_mass = num_blades * blade_mass;
    
    % hub calculations
    % after 2/3 of hub height, the circle tapers conically to a point
    L_cyl = (2/3) * hub_length;           
    L_cone = (1/3) * hub_length;    
   
    t_inner_kev = 0.0015;
    t_foam     = 0.028;
    t_outer_cf = 0.0035;
    
    r_inner = r_motor;
    r_kev   = r_inner + t_inner_kev;
    r_foam  = r_kev   + t_foam;
    r_cf2   = r_foam  + t_outer_cf;
    
    vol_kev_inner_cyl = pi * (r_kev^2 - r_inner^2) * L_cyl;
    vol_foam_cyl     = pi * (r_foam^2 - r_kev^2) * L_cyl;
    vol_cf_outer_cyl = pi * (r_cf2^2 - r_foam^2) * L_cyl;
    
    % Inner CF layer volume
    vol_kev_inner_cone = (1/3)*pi * L_cone * ...
        ((r_kev)^2 + r_kev*r_inner + (r_inner)^2);
    vol_inner_core_cone = (1/3)*pi * L_cone * ...
        ((r_inner)^2 + r_inner^2 + r_inner^2);  % same radius, essentially 0 volume
    vol_kev_inner_cone = vol_kev_inner_cone - vol_inner_core_cone;
    
    % Foam layer volume
    vol_foam_outer_cone = (1/3)*pi * L_cone * ...
        ((r_foam)^2 + r_foam*r_kev + (r_kev)^2);
    vol_foam_inner_cone = (1/3)*pi * L_cone * ...
        ((r_kev)^2 + r_kev^2 + r_kev^2);
    vol_foam_cone = vol_foam_outer_cone - vol_foam_inner_cone;
    
    % Outer CF layer volume
    vol_cf_outer_cone = (1/3)*pi * L_cone * ...
        ((r_cf2)^2 + r_cf2*r_foam + (r_foam)^2);
    vol_outer_foam_cone = (1/3)*pi * L_cone * ...
        ((r_foam)^2 + r_foam^2 + r_foam^2);
    
    vol_cf_outer_cone = vol_cf_outer_cone - vol_outer_foam_cone;
    
    hub_mass = rho_kevlar  * (vol_kev_inner_cyl + vol_kev_inner_cone) + ...
               rho_rohacell_foam * (vol_foam_cyl + vol_foam_cone) + ...
               carbonFiberDensity * (vol_cf_outer_cyl + vol_cf_outer_cone);
    
    duct_length = 0.4;  %m
    outer_radius = radius;
    kevThicknessDuct  = 0.0005;
    cfThicknessDuct  = 0.0005;
    foamThicknessDuct = duct_thickness_tot - kevThicknessDuct - cfThicknessDuct;
    
    inner_radius = radius - duct_thickness_tot;
    duct_volume = pi * (outer_radius^2 - inner_radius^2) * duct_length;
    duct_mass = duct_volume * (rho_kevlar * (kevThicknessDuct/duct_thickness_tot) +...
                               carbonFiberDensity * (cfThicknessDuct/duct_thickness_tot) + ...  
                               rho_rohacell_foam * (foamThicknessDuct/duct_thickness_tot));
    
    edf_mass_unit = total_blades_mass + hub_mass + duct_mass + motor_mass_unit;
    totalEDFWeightEst = edf_mass_unit * numEDF;
    
    % Misc. Weight Estimates
    miscElectronicsWeight = 1; %kg
    harnessGearWeight = 2; %kg
    
    nonWingMass = totalEDFWeightEst + batteryWeightEst + myWeight + miscElectronicsWeight;
    
    % wing spar calculations
    wall_t_spar = 1.5e-3;  % Spar wall thickness
    [r_outer_spar, r_inner_spar, diameter_spar] = findWingSparDimsMin((nonWingMass + wingWeight) * g, S, wall_t_spar);
    L_spar = wingspan_true - 0.2; % m, slightly less than full span to account fot tip losses
    A_spar = pi * (r_outer_spar^2 - r_inner_spar^2);
    V_spar = A_spar * L_spar;
    wingSparWeight = V_spar * carbonFiberDensity;
    
    wingMass = wingWeight + wingSparWeight;
    totalMass = nonWingMass + wingMass; %kg
    
    disp(['Component Masses:']);
    disp(['Total EDF Mass: ', num2str(totalEDFWeightEst), ' kg']);
    disp(['Battery Mass: ', num2str(batteryWeightEst), ' kg']);
    disp(['Human Mass: ', num2str(myWeight), ' kg']);
    disp(['Wing Mass: ', num2str(wingMass), ' kg']);
    disp(['Misc. Electronics Mass: ', num2str(miscElectronicsWeight), ' kg']);
    disp(['Harness Gear Mass: ', num2str(harnessGearWeight), ' kg']);
    disp(['-----------------------']);
    disp(['Total Estimated Mass: ', num2str(totalMass), ' kg']);
    disp(['-----------------------']);
    
    foamDensity = 20; %kg/m^3
    
    cfShellThickness_TaB = cfShellThickness*2; % top and bottom estimation
    
    % honeycomb wing
    % Calculate weighted average wing density using area values for half of the wing's
    % sections (this works because wing is symmetric)
    t_over_c = {
        varTable{strcmp(varTable.VariableName, "Wing t/c 0"), iteration}, ...
        varTable{strcmp(varTable.VariableName, "Wing t/c 1"), iteration}, ...
        varTable{strcmp(varTable.VariableName, "Wing t/c 2"), iteration}, ...
        varTable{strcmp(varTable.VariableName, "Wing t/c 3"), iteration}, ...
        varTable{strcmp(varTable.VariableName, "Wing t/c 4"), iteration}
    };
    
    S_sections = {
        varTable{strcmp(varTable.VariableName, "Section Area 1"), iteration}, ...
        varTable{strcmp(varTable.VariableName, "Section Area 2"), iteration}, ...
        varTable{strcmp(varTable.VariableName, "Section Area 3"), iteration}, ...
        varTable{strcmp(varTable.VariableName, "Section Area 4"), iteration}
    };
    
    % Convert to numeric arrays and filter out NaNs
    t_over_c_values = cell2mat(t_over_c);
    t_over_c_values = t_over_c_values(~isnan(t_over_c_values));
    
    area_values = cell2mat(S_sections);
    area_values = area_values(~isnan(area_values));
    
    densityTimesArea = zeros(numel(area_values), 1);
    wingThickness = zeros(numel(area_values), 1);
    for i = 1:numel(area_values)
        t_over_c_avg = (t_over_c_values(i) + t_over_c_values(i+1)) / 2;
        wingThickness(i) = t_over_c_avg * MAC_wing;
        foamThickness = wingThickness(i) - (cfShellThickness_TaB);
        wingDensity = carbonFiberDensity * (cfShellThickness_TaB/wingThickness(i)) + foamDensity * (foamThickness/wingThickness(i));
        densityTimesArea(i) = wingDensity * area_values(i);
    end
    totalWingDensity = sum(densityTimesArea) / sum(area_values);
    
    wingVolume = wingWeight / totalWingDensity;
    
    W = totalMass * g; % N
    T = W;
    T_unit = T / numEDF;
    
    wingloading_curr = W/S;
    
    disp(['For hover, thrust >= ', num2str(T), ' N']);
    disp(['For hover, thrust per EDF >= ', num2str(T_unit), ' N']);
    disp(['-----------------------']);
    
    % Calculate variables at altitude
    p_sl = 101325; %Pa
    rho_sl = 1.225; %kg/m3
    temp_sl = 288.15; %K
    R = 287; %J/(kg-K)
    lapse_rate = 0.0065; %K/m
    gamma = 1.4; % for dry air 
    
    temp_h = temp_sl - h_curr * lapse_rate; % K
    press_h = p_sl * (temp_h / temp_sl)^(g / (lapse_rate*R)); % Pa
    rho_h = press_h / (R * temp_h); % kg/m^3
    
    visc_sl = 1.716e-5; % kg/(ms) (reference dynamic viscosity at sea level)
    S_visc = 110.4; % K (sutherland constant for air)
    % Sutherland equation
    visc_h = visc_sl * ((temp_sl + S_visc) / (temp_h + S_visc)) * ((temp_h/temp_sl)^(3/2)); % kg/(ms)
    
    %v = sqrt(T_unit / (rho_sl * A_fan)); % m/s
    
    a = sqrt(gamma * R * temp_h);
    
    Mach = V / a;
    
    rho_curr = rho_h;
    
    Re = (rho_curr * V * MAC_wing) / visc_sl;
    COF_air = 1/sqrt(Re);
    q_TAS = 0.5 * rho_curr * V^2; % kg / ms^2 = Pa
    
    wingWeightGoal = (0.36 * S * q_TAS * AR_true) / totalMass;
    disp(['Wing Mass Goal = ', num2str(wingWeightGoal), ' kg']);
    disp(['Actual Wing Weight = ', num2str(wingWeight), ' kg']);
    
    disp(['Variables to continue iteration:']);
    disp(['Wing Density = ', num2str(totalWingDensity), ' kg/m^3']);
    disp(['Mach = ', num2str(Mach)]);
    disp(['Re = ', num2str(Re)]);
    disp('-------------------------');
    
    V_stall_actual = sqrt( (2 * wingloading_curr) / (rho_curr * CL_max) );
    
    iterationCol = find(strcmp(varTable.Properties.VariableNames, iteration));
    % Mass Values
    varTable{strcmp(varTable.VariableName, 'AoA Stall (deg)'), iterationCol} = AoA_stall;
    varTable{strcmp(varTable.VariableName, 'Total EDF Mass (kg)'), iterationCol} = totalEDFWeightEst;
    varTable{strcmp(varTable.VariableName, 'Battery mass (kg)'), iterationCol} = batteryWeightEst;
    varTable{strcmp(varTable.VariableName, 'Human Mass (kg)'), iterationCol} = myWeight;
    varTable{strcmp(varTable.VariableName, 'Wing Mass (including wing spar) (kg)'), iterationCol} = wingMass;
    varTable{strcmp(varTable.VariableName, 'Wing Density (kg/m^3)'), iterationCol} = totalWingDensity;
    varTable{strcmp(varTable.VariableName, 'Misc. Electronics Mass (kg)'), iterationCol} = miscElectronicsWeight;
    varTable{strcmp(varTable.VariableName, 'Harness Gear Mass (kg)'), iterationCol} = harnessGearWeight;
    varTable{strcmp(varTable.VariableName, 'Total Mass in System (kg)'), iterationCol} = totalMass;
    varTable{strcmp(varTable.VariableName, 'Min. Thrust for Hover (N)'), iterationCol} = T;
    varTable{strcmp(varTable.VariableName, 'Min. Thrust for hover per EDF (N)'), iterationCol} = T_unit;
    varTable{strcmp(varTable.VariableName, 'Wing Thickness (m)'), iterationCol} = strjoin(string(wingThickness), ',');
    varTable{strcmp(varTable.VariableName, 'Mach'), iterationCol} = Mach;
    varTable{strcmp(varTable.VariableName, 'C_L_Max'), iterationCol} = CL_max;
    varTable{strcmp(varTable.VariableName, 'Current W/S (N/m^2)'), iterationCol} = wingloading_curr;
    varTable{strcmp(varTable.VariableName, 'Re'), iterationCol} = Re;
    varTable{strcmp(varTable.VariableName, 'h (m)'), iterationCol} = h_curr;
    varTable{strcmp(varTable.VariableName, 'rho (kg/m^3)'), iterationCol} = rho_curr;
    varTable{strcmp(varTable.VariableName, 'Non-wing Weight (kg)'), iterationCol} = nonWingMass;
    varTable{strcmp(varTable.VariableName, 'V_stall Current (m/s)'), iterationCol} = V_stall_actual;
    
    
    writetable(varTable, 'ThesisWingOptimizationVariables.xlsx');
    
    %% 2. Flight Envelope (Power Required for Cruise and Optimal Glide)
    
    [~, numIterationPolars] = size(polarValues);
    currIteration = str2double(iteration)+1;
    filenameSet = polarValues(currIteration, :);
    
    % From Pod Man from OpenVSP
    C_d0_hum = 0.02133;

    for k = 1:numIterationPolars
        [alpha, C_L{k}, C_D_temp{k}, C_di{k}, C_D0_temp{k}, L_over_D{k}, E{k}] = ...
            getPolarArrays(filenameSet{k});
        format long
        C_D0{k} = C_D0_temp{k} + C_d0_hum;
        C_D{k} = C_D_temp{k} + C_d0_hum;
    end
    
    % Calculate AoA for level flight
    C_L_level_calc = W / (q_TAS*S);  
    pCL = polyfit(alpha, C_L{1}, 3);
    poly_func = @(alphaTarget) polyval(pCL, alphaTarget) - C_L_level_calc;
    level_alpha_degrees = fzero(poly_func, 10);
    disp('-----------------------------------------');
    
    % Calculate the glide angle (using AoA where L/D is maximized)
    Glide;

    V_stall_target = find_Vstall_target(V, level_alpha_degrees, AoA_stall);
    wingloading_max = 0.5 * rho_curr * CL_max * V_stall_target^2;
    if (V < V_stall_actual * 1.23)
        error('V should be less than V_stall actual!');
    end
    if V_stall_actual > V_stall_target
        error('Actual stall speed %.2f m/s exceeds target %.2f m/s!', V_stall_actual, V_stall_target);
    end
    disp(['W/S should be less than ', num2str(wingloading_max), ' N/m^2']);
    disp(['W/S is ', num2str(wingloading_curr), ' N/m^2']);
    
    % Fits for plotting
    [pCL_WH, pCD_WH, pCDi_WH, pCD0_WH, pLD_WH, pEff_WH, ... 
    C_L_fit_WH, C_D_fit_WH, C_di_fit_WH, C_d0_fit_WH, L_D_fit_WH, eff_fit_WH] = ...
        AeroFits.AeroFitsFromAlpha(alpha, C_L{1}, C_D{1}, C_di{1}, C_D0{1}, L_over_D{1}, E{1});
    [pCL_WP, pCD_WP, pCDi_WP, pCD0_WP, pLD_WP, pEff_WP, ... 
    C_L_fit_WP, C_D_fit_WP, C_di_fit_WP, C_d0_fit_WP, L_D_fit_WP, eff_fit_WP] = ...
        AeroFits.AeroFitsFromAlpha(alpha, C_L{2}, C_D{2}, C_di{2}, C_D0{2}, L_over_D{2}, E{2});
    
    % For Level Flight for WH
    [C_L_level_WH, C_D_level_WH, C_di_level_WH, C_d0_level_WH, L_D_level_WH, eff_level_WH] = ...
        AeroCoeffs.AeroCoeffsFromAlpha(level_alpha_degrees, alpha, C_L{1}, C_D{1}, C_di{1}, C_D0{1}, L_over_D{1}, E{1});
    % For Level Flight for WP
    [C_L_level_WP, C_D_level_WP, C_di_level_WP, C_d0_level_WP, L_D_level_WP, eff_level_WP] = ...
        AeroCoeffs.AeroCoeffsFromAlpha(level_alpha_degrees, alpha, C_L{2}, C_D{2}, C_di{2}, C_D0{2}, L_over_D{2}, E{2});
    % For Glide Flight for WH
    [C_L_glide_WH, C_D_glide_WH, C_di_glide_WH, C_d0_glide_WH, L_D_glide_WH, eff_glide_WH] = ...
        AeroCoeffs.AeroCoeffsFromAlpha(alpha_glide, alpha, C_L{1}, C_D{1}, C_di{1}, C_D0{1}, L_over_D{1}, E{1});
    % For Glide Flight for WP
    [C_L_glide_WP, C_D_glide_WP, C_di_glide_WP, C_d0_glide_WP, L_D_glide_WP, eff_glide_WP] = ...
        AeroCoeffs.AeroCoeffsFromAlpha(alpha_glide, alpha, C_L{2}, C_D{2}, C_di{2}, C_D0{2}, L_over_D{2}, E{2});
    
    % Prints for Confirmation of Aeordynamic Coefficients 
    disp('For Glide Flight for WH:');
    disp(['C_L = ', num2str(C_L_glide_WH)]);
    disp(['C_D = ', num2str(C_D_glide_WH)]);
    disp(['C_D0 = ', num2str(C_d0_glide_WH)]);
    disp(['C_Di = ', num2str(C_di_glide_WH)]);
    disp('-----------------------------------------');
    disp('For Level Flight for WH:');
    disp(['Calculated C_L = ', num2str(C_L_level_calc)]);
    disp(['Alpha = ', num2str(level_alpha_degrees)]);
    disp(['C_D = ', num2str(C_D_level_WH)]);
    disp(['C_D0 = ', num2str(C_d0_level_WH)]);
    disp(['C_di = ', num2str(C_di_level_WH)]);
    disp('-----------------------------------------');

    % WH: Lift Coefficient (C_L) vs angle of attack (α)
    figure;
    plot(alpha, C_L{1}, 'bo', 'MarkerFaceColor', 'b'); 
    hold on;
    plot(alpha, C_L_fit_WH, 'r-', 'LineWidth', 2); 
    hold on;
    plot(level_alpha_degrees*ones(numel(alpha), 1), C_L_fit_WH, 'k--');
    hold on;
    plot(alpha_glide*ones(numel(alpha), 1), C_L_fit_WH, 'r--');
    hold on;
    title('Lift Coefficient (C_L) vs angle of attack (α)');
    xlabel('α');
    ylabel('C_L');
    legend('C_L Data', 'C_L Fitted Line', ['Level flight α = ', num2str(level_alpha_degrees), '°'], ['Steady glide α = ', num2str(glide_angle_degrees), '°']);
    grid on;
    
    % WH: Drag and Parasitic Drag Coefficients (C_D and CD0) vs angle of attack (α)'
    figure;
    plot(alpha, C_D{1}, 'bo', 'MarkerFaceColor', 'b'); 
    hold on;
    plot(alpha, C_D_fit_WH, 'r-', 'LineWidth', 2); 
    hold on;
    plot(alpha, C_D0{1}, 'mo', 'MarkerFaceColor', 'm'); 
    plot(alpha, C_d0_fit_WH, 'g-', 'LineWidth', 2); 
    plot(level_alpha_degrees*ones(numel(alpha), 1), C_D_fit_WH, 'k--');
    hold on;
    plot(alpha_glide*ones(numel(alpha), 1), C_D_fit_WH, 'r--');
    hold on;
    title('Drag and Parasitic Drag Coefficients (C_D and CD0) vs angle of attack (α)');
    xlabel('α');
    ylabel('Coefficient Magnitude');
    legend('C_D Data', 'C_D Fitted Line', 'CD0 Data', 'CD0 Fitted Line', 'Level flight α', 'Steady glide α');
    grid on;
    
    % WH: Plot W/P vs W/S
    wingloadingArray = 0:10:(wingloading_max+50);
    for k = 1:numel(wingloadingArray)
        % Calculate AoA for level flight
        C_L_level_calc_plot = (2 * (wingloadingArray(k))) / (rho_curr * V^2);  
        C_L_level_calc_plots(k) = C_L_level_calc_plot; 
        poly_func = @(alpha) polyval(pCL_WH, alpha) - C_L_level_calc_plot;
        level_alpha_degrees_plot = fzero(poly_func, 10);  
        level_alpha_degrees_plots(k) = level_alpha_degrees_plot;
        [C_L_level_WH_plot, C_D_level_WH_plot, C_di_level_WH_plot(k), C_d0_level_WH_plot(k), L_D_level_WH_plot, eff_level_WH_plot] = ...
        AeroCoeffs.AeroCoeffsFromAlpha(level_alpha_degrees_plot, alpha, C_L{1}, C_D{1}, C_di{1}, C_D0{1}, L_over_D{1}, E{1});
    
        % Power Stuff
        % For level flight (T = D and L = W)
        C_D_level_plot_WH(k) = C_d0_level_WH_plot(k) + C_di_level_WH_plot(k);
        D_level_plot(k) = q_TAS * S * C_D_level_plot_WH(k);
        P_req_level(k) = D_level_plot(k) * V;
        W_P(k) = W / P_req_level(k);
        %W_P_new(k) = wingloading_curr / (0.5*rho_curr*V^2*C_d0_level_WH + );
    end
    
    % For steady glide
    C_D_glide_calc = C_d0_glide_WH + C_di_glide_WH;
    D_glide = q_TAS * S * C_D_glide_calc;
    % derived from C_L @ max L/D = sqrt(Cdo/k2)
    k_2_glide =  C_d0_glide_WH / C_L_glide_WH^2;
    velo_glide_optimal = sqrt(2*wingloading_curr / rho_curr) * ((k_2_glide/(3*C_d0_glide_WH))^(1/4));
    
    % For level flight (T = D and L = W)
    C_D_level_calc = C_d0_level_WH + C_di_level_WH;
    D_level = q_TAS * S * C_D_level_calc;
    P_req_level_curr = D_level * V;
    P_avail_level_curr = P_req_level_curr;
    
    % Sweeping number of EDFs to assess drag at cruise:
    numEDFSweep = 2:6;
    for i = 1:numel(numEDFSweep)
        W_edfSweep = (W - (totalEDFWeightEst*g))  + (edf_mass_unit * g * numEDFSweep(i));
        C_L_level_calcs = W_edfSweep / (q_TAS*S);  
        pCLs = polyfit(alpha, C_L{1}, 3);
        poly_func = @(alphaTarget) polyval(pCLs, alphaTarget) - C_L_level_calcs;
        level_alpha_degrees_S = fzero(poly_func, 10);
        [C_L_level_WHs, C_D_level_WHs, C_di_level_WHs, C_d0_level_WHs, L_D_level_WHs, eff_level_WHs] = ...
        AeroCoeffs.AeroCoeffsFromAlpha(level_alpha_degrees_S, alpha, C_L{1}, C_D{1}, C_di{1}, C_D0{1}, L_over_D{1}, E{1});
        C_D_level_calcs = C_d0_level_WHs + C_di_level_WHs;
        D_levels = q_TAS * S * C_D_level_calcs;
        T_cruises(i) = D_levels;
    end
    
    wingloadingYRange1 = 0:max(W_P) / (numel(wingloadingArray) - 1):max(W_P);
figure;
plot(wingloadingArray, W_P, 'b-', 'LineWidth', 1.5);
hold on;

% Vertical lines for wingloading_max and wingloading_curr
plot(ones(numel(wingloadingArray), 1) * wingloading_max, wingloadingYRange1, 'r-', 'LineWidth', 2);
hold on;
plot(ones(numel(wingloadingArray), 1) * wingloading_curr, wingloadingYRange1, 'm--', 'LineWidth', 2);
hold on;
% Shade area under curve and to the left of wingloading_max
idx = wingloadingArray <= wingloading_max;
x_fill = [wingloadingArray(idx), fliplr(wingloadingArray(idx))];
y_fill = [W_P(idx), zeros(1, sum(idx))]; 
fill(x_fill, y_fill, [0.8 0.9 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

legendEntries = {
    'Cruise Condition', ...
    sprintf('Stall Wingloading = %.1f N/m^2', wingloading_max), ...
    sprintf('Current Wingloading = %.1f N/m^2', wingloading_curr), ...
    'Flight Envelope'
};

legend(legendEntries, 'Location', 'best');

title('Flight Envelope');
xlabel('W/S (N/m^2)');
ylim([0, 1]);
ylabel('W/P (N/W)');
grid on;
    wingloadingYRange2 = 0:max(P_req_level) / (numel(wingloadingArray) - 1):max(P_req_level);
    figure;
    plot(wingloadingArray, P_req_level);
    hold on;
    plot(ones(numel(wingloadingArray), 1) * wingloading_max, wingloadingYRange2, 'LineWidth', 2);
    hold on;
    plot(ones(numel(wingloadingArray), 1) * wingloading_curr, wingloadingYRange2, 'LineWidth', 2);
    title('Power vs W/S for Level Flight');
    xlabel('W/S (N/m^2)');
    ylabel('Power (W)');
    grid on;
    
eff_WH_cruise = C_L_level_WH / (pi * AR_proj * C_di_level_WH);
eff_WH = C_L{1}.^2 ./ (pi .* AR_proj .* C_di{1});
eff_WP = C_L{2}.^2 ./ (pi .* AR_proj .* C_di{2});
    

    % Verification of VLM method vs Panel Method
    if iteration == "42"
        verificationFilesTemp = {'polar_42_W_VLM.txt', 'polar_42_W_P.txt'};
        prefix = 'PolarData/';  
        verificationFiles = cellfun(@(x) [prefix, x], verificationFilesTemp, 'UniformOutput', false);
        figure;
        for i = 1:2
            [alphaver, C_Lver{i}, C_Dver{i}, C_diver{i}, C_D0ver{i}, L_over_Dver{i}, Ever{i}] = getPolarArrays(verificationFiles{i});
            format long
            % Find C_L vs alpha curve

    [pCLver, pCDver, pCDiver, pCD0ver, pLDver, pEffver, ...
        C_L_fit_ver{i}, C_D_fit_ver{i}, C_di_fit_ver{i}, C_d0_fit_ver{i}, L_D_fit_ver{i}, eff_fit_ver{i}] = ...
        AeroFits.AeroFitsFromAlpha(alphaver,  C_Lver{i}, C_Dver{i},  C_diver{i}, C_D0ver{i}, L_over_Dver{i}, Ever{i});
    
    pCL = polyfit(alphaver, C_Lver{i}, 3);
    poly_func = @(alphaTarget) polyval(pCL, alphaTarget) - C_L_level_calc;
    level_alpha_degrees = fzero(poly_func, 10);
    
    C_d0_level_ver(i) = AeroCoeffs.YfromAlpha(alphaver, C_D0ver{i}, level_alpha_degrees);
    C_d0_glide_ver(i) = AeroCoeffs.YfromAlpha(alphaver, C_D0ver{i}, glide_angle_degrees);
    
    C_d_level_ver(i) = AeroCoeffs.YfromAlpha(alphaver, C_Dver{i}, level_alpha_degrees);
    C_d_glide_ver(i) = AeroCoeffs.YfromAlpha(alphaver, C_Dver{i}, glide_angle_degrees);
    
    C_L_level_ver(i) = AeroCoeffs.YfromAlpha(alphaver, C_Lver{i}, level_alpha_degrees);
    C_L_glide_ver(i) = AeroCoeffs.YfromAlpha(alphaver, C_Lver{i}, glide_angle_degrees);
    
    C_di_level_ver(i) = AeroCoeffs.YfromAlpha(alphaver, C_diver{i}, level_alpha_degrees);
    C_di_glide_ver(i) = AeroCoeffs.YfromAlpha(alphaver, C_diver{i}, glide_angle_degrees);
       
    e_level_ver(i) = AeroCoeffs.YfromAlpha(alphaver, Ever{i}, level_alpha_degrees);
    e_glide_ver(i) = AeroCoeffs.YfromAlpha(alphaver, Ever{i}, glide_angle_degrees);
       
    
            % Find CD vs CL polar
            [pCD, C_Dver_fit] = AeroFits.FitfromAlpha(C_Lver{i}, C_Dver{i});
            plot(C_Lver{i}, C_Dver{i}, 'bo', 'MarkerFaceColor', 'b'); 
            hold on;
            plot(C_Lver{i}, C_Dver_fit, 'r-', 'LineWidth', 2); 
            hold on;
        end
        title('C_D vs C_L Polar Curve');
        xlabel('C_L');
        ylabel('C_D');
        legend('VLM Data', 'VLM Fitted Curve', 'Panel Data', 'Panel Fitted Curve');
        grid on;
figure;

% Row 1: C_di 
subplot(4,3,1); 
for i = 1:2
    plot(alphaver, C_diver{i}, 'o-', 'LineWidth', 1.5); hold on;
    plot(alphaver, C_di_fit_ver{i}, '-', 'LineWidth', 2);
end
title('C_{di} vs \alpha');
ylabel('C_{di}'); grid on;
legend('VLM Data', 'VLM Fit', 'Panel Data', 'Panel Fit', 'FontSize', 7);

subplot(4,3,2); 
Cdi_diff = abs(C_diver{1} - C_diver{2});
plot(alphaver, Cdi_diff, 'k*-', 'LineWidth', 1.5);
title('|\Delta C_{di}|');
ylabel('|\Delta|'); grid on;

subplot(4,3,3); 
Cdi_max_total = max([abs(C_diver{1}), abs(C_diver{2})], [], 'all');
Cdi_pct_diff = 100 * Cdi_diff ./ Cdi_max_total;
plot(alphaver, Cdi_pct_diff, 'm*-', 'LineWidth', 1.5);
title('% \Delta C_{di} (global norm)');
ylabel('% Diff'); grid on;

% Row 2: C_L 
subplot(4,3,4); 
for i = 1:2
    plot(alphaver, C_Lver{i}, 'o-', 'LineWidth', 1.5); hold on;
    plot(alphaver, C_L_fit_ver{i}, '-', 'LineWidth', 2);
end
title('C_L vs \alpha');
ylabel('C_L'); grid on;
legend('VLM Data', 'VLM Fit', 'Panel Data', 'Panel Fit', 'FontSize', 7);

subplot(4,3,5); 
CL_diff = abs(C_Lver{1} - C_Lver{2});
plot(alphaver, CL_diff, 'k*-', 'LineWidth', 1.5);
title('|\Delta C_L|');
ylabel('|\Delta|'); grid on;

subplot(4,3,6); 
CL_max_total = max([abs(C_Lver{1}), abs(C_Lver{2})], [], 'all');
CL_pct_diff = 100 * CL_diff ./ CL_max_total;
plot(alphaver, CL_pct_diff, 'm*-', 'LineWidth', 1.5);
title('% \Delta C_L (global norm)');
ylabel('% Diff'); grid on;

% Row 3: C_D0 
subplot(4,3,7); 
for i = 1:2
    plot(alphaver, C_D0ver{i}, 'o-', 'LineWidth', 1.5); hold on;
    plot(alphaver, C_d0_fit_ver{i}, '-', 'LineWidth', 2);
end
title('C_{D0} vs \alpha');
ylabel('C_{D0}'); grid on;
legend('VLM Data', 'VLM Fit', 'Panel Data', 'Panel Fit', 'FontSize', 7);

subplot(4,3,8); 
CD0_diff = abs(C_D0ver{1} - C_D0ver{2});
plot(alphaver, CD0_diff, 'k*-', 'LineWidth', 1.5);
title('|\Delta C_{D0}|');
ylabel('|\Delta|'); grid on;

subplot(4,3,9); 
CD0_max_total = max([abs(C_D0ver{1}), abs(C_D0ver{2})], [], 'all');
CD0_pct_diff = 100 * CD0_diff ./ CD0_max_total;
plot(alphaver, CD0_pct_diff, 'm*-', 'LineWidth', 1.5);
title('% \Delta C_{D0} (global norm)');
ylabel('% Diff'); grid on;

% Row 4: L/D 
subplot(4,3,10); 
hold on;
h1 = plot(alphaver, L_over_Dver{1}, 'o-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5); % Blue - VLM Data
h2 = plot(alphaver, L_over_Dver{1}, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2); % Red - VLM Fit (placeholder)
h3 = plot(alphaver, L_over_Dver{2}, 'o-', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5); % Yellow - Panel Data
h4 = plot(alphaver, L_over_Dver{2}, '-', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2); % Purple - Panel Fit (placeholder)
title('L/D vs \alpha');
xlabel('\alpha (deg)'); ylabel('L/D'); grid on;
legend([h1 h2 h3 h4], 'VLM Data', 'VLM Fit', 'Panel Data', 'Panel Fit', 'FontSize', 7);

subplot(4,3,11); 
LD_diff = abs(L_over_Dver{1} - L_over_Dver{2});
plot(alphaver, LD_diff, 'k*-', 'LineWidth', 1.5);
title('|\Delta L/D|');
xlabel('\alpha (deg)'); ylabel('|\Delta|'); grid on;

subplot(4,3,12); 
LD_max_total = max([abs(L_over_Dver{1}), abs(L_over_Dver{2})], [], 'all');
LD_pct_diff = 100 * LD_diff ./ LD_max_total;
plot(alphaver, LD_pct_diff, 'm*-', 'LineWidth', 1.5);
title('% \Delta L/D (global norm)');
xlabel('\alpha (deg)'); ylabel('% Diff'); grid on;
sgtitle('Difference in Aerodynamic Parameters of Wing Only: Panel Method vs VLM');
    
    eff_W_P = C_Lver{2}.^2 ./ (pi .* AR_proj .* C_diver{2});
    eff_W_VLM = C_Lver{1}.^2 ./ (pi .* AR_proj .* C_diver{1});

% Create subplot figure
figure;

% Left subplot - Panel Method
subplot(1, 2, 1);
plot(alphaver, eff_W_P, '-o', 'LineWidth', 1.5);
title('Panel Method');
xlabel('\alpha (deg)');
ylabel('Efficiency');
grid on;

% Right subplot - VLM Method
subplot(1, 2, 2);
plot(alphaver, eff_W_VLM, '-o', 'LineWidth', 1.5);
title('VLM');
xlabel('\alpha (deg)');
ylabel('Efficiency');
grid on;
sgtitle('Comparison of Oswald Efficiency of Wing Only: Panel Method vs VLM');

% Create the plot
figure;
plot(alphaver, C_D0_temp{1}, '-o');     
hold on;
plot(alphaver, C_D0ver{2}, '-s'); 
hold off;
xlabel('\alpha (deg)');
ylabel('$C_{D_0}$', 'Interpreter', 'latex');
title('Parasitic Drag vs AoA simulated with Panel Method');
legend('Wing and Realistic Human', 'Wing Only'); grid on;
    
figure;

% Col 1: Raw panel method values 
subplot(1,3,1); % 1 row, 3 columns, position 1
plot(alphaver, C_D0_temp{1}, '-o', 'LineWidth', 1.5); hold on;
plot(alphaver, C_D0ver{2}, '-s', 'LineWidth', 1.5);
xlabel('\alpha (deg)');
ylabel('$C_{D_0}$', 'Interpreter', 'latex');
title('Parasitic Drag vs AoA');
legend('Wing and Human', 'Wing Only', 'FontSize', 7);
grid on;

% Col 2: Absolute difference 
subplot(1,3,2); % position 2
CD0_diff = abs(C_D0_temp{1} - C_D0ver{2});
plot(alphaver, CD0_diff, 'k*-', 'LineWidth', 1.5);
xlabel('\alpha (deg)');
ylabel('$|\Delta C_{D_0}|$', 'Interpreter', 'latex');
title('Absolute Difference');
grid on;

% Col 3: Percent difference (global normalized) 
subplot(1,3,3); 
CD0_max_total = max([abs(C_D0_temp{1}), abs(C_D0ver{2})], [], 'all');
CD0_pct_diff = 100 * CD0_diff ./ CD0_max_total;
plot(alphaver, CD0_pct_diff, 'm*-', 'LineWidth', 1.5);
xlabel('\alpha (deg)');
ylabel('$\%\ \Delta C_{D_0}$', 'Interpreter', 'latex');
title('Percent Difference (Global Norm)');
grid on;
sgtitle('Parasitic Drag vs AoA simulated with Panel Method');

    end
    
    % Get Steady Glide Variables into varTable
    varTable{strcmp(varTable.VariableName, 'Max W/S (N/m^2)'), iterationCol} = wingloading_max;
    varTable{strcmp(varTable.VariableName, 'V_stall Target (m/s)'), iterationCol} = V_stall_target;
    
    varTable{strcmp(varTable.VariableName, 'Max L/D'), iterationCol} = L_over_D_max;
    varTable{strcmp(varTable.VariableName, 'Optimal Glide Angle (deg)'), iterationCol} = glide_angle_degrees;
    varTable{strcmp(varTable.VariableName, 'Steady Glide AoA (deg)'), iterationCol} = alpha_glide;
    varTable{strcmp(varTable.VariableName, 'Glide C_L'), iterationCol} = C_L_glide_WH;
    varTable{strcmp(varTable.VariableName, 'Glide C_D0'), iterationCol} = C_d0_glide_WH;
    varTable{strcmp(varTable.VariableName, 'Glide C_Di Simulated'), iterationCol} = C_di_glide_WH;
    %varTable{strcmp(varTable.VariableName, 'Glide C_Di Calculated'), iterationCol} = C_di_glide_calc;
    varTable{strcmp(varTable.VariableName, 'Glide C_D'), iterationCol} = C_D_glide_WH;
    varTable{strcmp(varTable.VariableName, 'Glide D (N)'), iterationCol} = D_glide;
    varTable{strcmp(varTable.VariableName, 'Optimal Glide Velocity (m/s)'), iterationCol} = velo_glide_optimal;
    % Get Level Flight Variables into varTable
    varTable{strcmp(varTable.VariableName, 'Level Flight AoA (deg)'), iterationCol} = level_alpha_degrees;
    varTable{strcmp(varTable.VariableName, 'Calculated level C_L'), iterationCol} = C_L_level_calc;
    varTable{strcmp(varTable.VariableName, 'Level C_D0'), iterationCol} = C_d0_level_WH;
    varTable{strcmp(varTable.VariableName, 'Level C_Di Simulated'), iterationCol} = C_di_level_WH;
    %varTable{strcmp(varTable.VariableName, 'Level C_Di Calculated'), iterationCol} = C_di_level_calc;
    varTable{strcmp(varTable.VariableName, 'Level C_D'), iterationCol} = C_D_level_WH;
    varTable{strcmp(varTable.VariableName, 'Level D (N)'), iterationCol} = D_level;
    varTable{strcmp(varTable.VariableName, 'Power Required for Level Flight (W)'), iterationCol} = P_req_level_curr;
    %varTable{strcmp(varTable.VariableName, 'Efficiency for Wing'), iterationCol} = eff_W;
    %varTable{strcmp(varTable.VariableName, 'Efficiency for Wing / Human'), iterationCol} = eff_WH;
    
    writetable(varTable, 'ThesisWingOptimizationVariables.xlsx');
    disp(varTable);
%clear; close all;
%end

%% Functions
function [r_outer, r_inner, best_diameter] = findWingSparDimsMin(W, S, t)
    M_max = (W * S) / 4;  % Max bending moment
    
    % Carbon fiber properties
    sigma_max = 400e6;  % Max allowable stress (SF = 2, actual ~800 MPa)
    
    % Finding the minimum required diameter
    for d_outer = 0.025:0.001:0.1  % diameters from 35mm to 80mm
        r_outer = d_outer / 2;
        r_inner = r_outer - t;
        I = (pi / 4) * (r_outer^4 - r_inner^4);  % Second moment of area
        sigma = (M_max * r_outer) / I;  % Bending stress
        
        if sigma < sigma_max
            best_diameter = d_outer; % m
            break;
        end
    end
end

function V_stall = find_Vstall_target(CurrentV, CurrentAoA, AoAstall)
    a = CurrentV^2 * CurrentAoA;
    V_stall = sqrt(a/AoAstall);
end

function range = createRange(start_val, end_val, numEls)
    if numEls < 2
        range = start_val; 
    else
        step = (end_val - start_val) / (numEls - 1);
        range = start_val + step * (0:numEls - 1);
    end
end

function [alpha, C_L, C_D, C_di, C_D0, L_over_D, E] = getPolarArrays(filename)

    % Gets aerodynamic coefficients for WHOLE SYSTEM CONFIGURATION
    fid = fopen(filename, 'r');
    headerLine = fgetl(fid); 
    headerCellsTemp = strsplit(headerLine);
    headerCells = headerCellsTemp(~cellfun('isempty', headerCellsTemp));
    % Identify the column indices
    alphaIdx = find(strcmpi(headerCells, 'AoA'));
    ClIdx = find(strcmpi(headerCells, 'CL'));
    CdIdx = find(strcmpi(headerCells, 'CDtot'));
    Cd0Idx = find(strcmpi(headerCells, 'CDo'));
    CdiIdx = find(strcmpi(headerCells, 'CDi'));
    LoverDIdx = find(strcmpi(headerCells, 'L/D'));
    EIdx = find(strcmpi(headerCells, 'E'));
    % Read the rest of the data (skip header) 
    data = textscan(fid, repmat('%f', 1, numel(headerCells)), 'Delimiter', '\t'); 
    fclose(fid);
    
    % Extract the data for the coefficients 
    alpha = transpose(data{alphaIdx}); 
    C_L = transpose(data{ClIdx}); 
    C_D = transpose(data{CdIdx}); 
    C_di = transpose(data{CdiIdx});
    C_D0 = transpose(data{Cd0Idx});
    L_over_D = transpose(data{LoverDIdx});
    E = transpose(data{EIdx});
end