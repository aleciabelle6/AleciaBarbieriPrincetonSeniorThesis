syms T_sl W_to P_s q_TAS S b a k2 n CD0 VTAS

% Define the equation using == for equality
lhs = T_sl / W_to;
rhs = (b / a) * ( ( ((q_TAS * S)/(b * W_to) )) * (k2 * ((n * b * W_to) / (q_TAS * S))^2 + CD0) + (P_s / VTAS));

% Set up the equation lhs == rhs
eq = lhs == rhs;

% Solve for W_to
P_to_expr = solve(eq, P_s);
P_to_expr_simplify = simplify(P_to_expr);
% Solve for W_to / P_s
W_over_P_s = W_to / P_to_expr;

% Simplify the expression
W_over_P_s_simplified = simplify(W_over_P_s);

%W_to_over_P_s_numeric = subs(W_to_over_P_s_simplified, ...
%    {W_to, T_sl, q_TAS, S, b, a, k2, n, CD0, VTAS}, ...
%    {W_to_val, T_sl_val, q_TAS_val, S_val, beta_val, alpha_val, k2_val, n_val, CD0_val, VTAS_val});

%result = double(W_to_over_P_s_numeric);