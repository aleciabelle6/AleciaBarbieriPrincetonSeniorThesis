%% Glide Angle Calculations
% Calculate glide angle

% Fiting L/D vs alpha curve
[xData, yData] = prepareCurveData(alpha, L_over_D{2});
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.999896350047154;
[fitresult, gof] = fit( xData, yData, ft, opts );
d_fitted_curve = differentiate(fitresult, xData);

% Find the maximum DC_ value and the corresponding C_L
alpha_range = linspace(min(xData), max(xData), 100);
L_over_D_vals = feval(fitresult, alpha_range);
[L_over_D_max, idx_max] = max(L_over_D_vals);
alpha_glide = alpha_range(idx_max);

disp(['Maximum L/D occurs at α = ', num2str(alpha_glide)]);
disp(['(L/D)_max = ', num2str(L_over_D_max)]);



% M=0.08, slip angle, beta=0
% Plot fit with data.
figure;
plot(alpha, L_over_D{2}, 'bo', 'MarkerFaceColor', 'b');
hold on;
hPlot = plot( fitresult, xData, yData);
hold on; 
hGlide = plot(alpha_glide*ones(numel(alpha),1), createRange(0, 25, numel(alpha)), 'r--', 'LineWidth', 1.5);
hold on;
hLevel = plot(level_alpha_degrees*ones(numel(alpha), 1), createRange(0, 25, numel(alpha)), 'k--', 'LineWidth', 1.5);
title("L/D vs. Angle of Attack (α)");
legend([hPlot; hGlide; hLevel], ...
       {'Data', ...
        'Fitted Curve', ...
        ['Steady glide α = ', num2str(alpha_glide), '°'], ...
        ['Level flight α = ', num2str(level_alpha_degrees), '°']}, ...
       'Location', 'NorthEast', 'Interpreter', 'none');
set(hPlot, 'LineWidth', 2);
xlabel( 'α (degrees)', 'Interpreter', 'none' );
ylabel( 'L/D', 'Interpreter', 'none' );
grid on

% Fit a 3rd degree polynomial to C_L vs alpha
[C_L_target] = AeroCoeffs.YfromAlpha(alpha, C_L{1}, alpha_glide);
[pCL, C_L_fit2] = AeroFits.FitfromAlpha(alpha, C_L{1});
disp(['C_L for (L/D)_max = ', num2str(C_L_target)]);

% Fit a 3rd degree polynomial to C_L vs C_D
[C_D_target] = AeroCoeffs.YfromAlpha(C_L{1}, C_D{1}, C_L_target);
[pCD, C_D_fit] = AeroFits.FitfromAlpha(C_L{1}, C_D{1});
disp(['C_D for (L/D)_max = ', num2str(C_D_target)]);

[C_D0_target] = AeroCoeffs.YfromAlpha(alpha, C_D0{1}, C_L_target);
[pCD0, C_D0fit3] = AeroFits.FitfromAlpha(alpha, C_D0{1});
disp(['C_D0 for (L/D)_max = ', num2str(C_D0_target)]);

% Find glide angle with tangent
pCD_CL = polyfit(C_L{1}, C_D{1}, 3);
dp = polyder(pCD_CL);
slope = polyval(dp, C_L_target);  % Evaluate the derivative at C_L_target
y0 = polyval(pCD_CL, C_L_target);  % Evaluate the original polynomial at C_L_target
x_tangent = linspace(C_L_target-0.5, C_L_target+0.7, 50);  
y_tangent_line = slope * (x_tangent - C_L_target) + y0;
glide_angle_radians = atan(slope);  
glide_angle_degrees = rad2deg(glide_angle_radians);
disp(['glide angle (in degrees) = ', num2str(glide_angle_degrees)]);


C_L_opp_range_polar = createRange(-0.1, C_L_target, 10);
C_D_opp_range_polar = createRange(0, C_D_target, 10);

figure;
plot(C_L{1}, C_D{1}, 'bo', 'MarkerSize', 7); 
hold on;
plot(C_L{1}, C_D_fit, 'r-', 'LineWidth', 2);
hold on;
plot(C_L_target, C_D_target, 'bo', 'MarkerSize', 7, 'MarkerFaceColor', 'b');
hold on;
plot(x_tangent, y_tangent_line, 'LineWidth', 2, 'Color', [0.92, 0.72, 0.18]);  % Plot the tangent
hold on;
plot(C_L_target*ones(numel(alpha), 1), C_D_opp_range_polar, 'r--', 'LineWidth', 1);
hold on;
plot(C_L_opp_range_polar, C_D_target*ones(numel(alpha), 1), 'r--', 'LineWidth', 1);
hold on;
title('C_D vs C_L Polar Curve');
xlabel('C_L');
ylabel('C_D');
xlim([-0.1 1.3]);
legend('Data', ...
       'Fitted Curve', ...
       ['(C_L, C_D) at max (L/D) = (', num2str(C_L_target, '%.3f'), ', ',...
       num2str(C_D_target, '%.4f'), ')'], ...
       'Tangent Line at max (L/D)');
grid on;

fprintf('The fitted polynomial function is: C_D = %.4f + %.4f*C_L + %.4f*C_L^2 + %.4f*C_L^3\n', pCD_CL(4), pCD_CL(3), pCD_CL(2), pCD_CL(1));


function range = createRange(start_val, end_val, numEls)
    if numEls < 2
        range = start_val;
    else
        step = (end_val - start_val) / (numEls - 1);
        range = start_val + step * (0:numEls - 1);
    end
end