%% 7. PowerRequired (Calculates Power Requirements for our sizing our Power Source)

% Call 6. PropDesign
PropDesign;

TransitionPeriod;

% Calculate Power Required for Transitional Periods between Glide & Cruise

V = 43.2; % V

segmentNames = ["Climb", "Still Hover", "Cruise"];
powerValues = [power_true_total_hoverVc, power_true_total_hoverStill, 2437.7];
timeSeconds = [50, 20, 60*10];

E_J_G2C = E_wh_G2C * 3600;
E_J_C2G = E_wh_C2G * 3600;

energiesNonTransition = powerValues .* timeSeconds;

Joules = sum(energiesNonTransition) + E_J_G2C + E_J_C2G;

energy_WH = Joules / 3600; % s -> hr

C = Joules / V;

mAh = C * 0.277778;
disp(['Battery rating: ', num2str(mAh), ' mAh']) 

current_mA = mAh / (sum(timeSeconds) * 3600);
current_A = current_mA / 1000;

% Throttle Calculations
throttle_Cruise =  power_true_unit_cruise / P_max_motor;
throttle_hoverStill = powerActual_unit_hoverStill / P_max_motor;
throttle_hoverVc = powerActual_unit_hoverVc / P_max_motor;
disp(['Cruise Throttle: ', num2str(throttle_Cruise*100), '%']) 
disp(['Hover Throttle: ', num2str(throttle_hoverStill*100), '%']) 
disp(['Climb Throttle: ', num2str(throttle_hoverVc*100), '%']) 
