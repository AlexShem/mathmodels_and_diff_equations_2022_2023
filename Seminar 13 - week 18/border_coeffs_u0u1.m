clear
clc

%% Declare symbolic variables
syms x t u(t, x) f(t, x)
syms D C h tau positive
syms alpha beta [3 4]
syms nu mu

alpha_mat = alpha;
alpha_mat(2, 1:2) = 0;
alpha_mat(3, :) = 0;
alpha = nonzeros(alpha_mat);

beta_mat = beta;
beta_mat(1, [3:end]) = 0;
beta_mat(2, 1:end) = 0;
beta_mat(3, :) = 0;
beta = nonzeros(beta_mat);

%% Define test functions
test_funs = [x^2, x^3, x^4, ...
    t*x^2, t*x^3,...
    t^2*x^2]; % Rank =

% h-tau offset grid
[h_mat, tau_mat] = meshgrid((0:3)*h, (2:-1:0)*tau);

% Rod differential operator
P = @(U) diff(U, t, 2) - D*diff(diff(U, x, 2), t, 2) + C*diff(U, x, 4);

%% Build 2 systems
coeffs = sym('Eqs', [length(test_funs) + 2, 1]);

for k = 1 : numel(test_funs)
    u(t, x) = test_funs(k);
    f(t, x) = P(u);
    
    u_compact = sum(alpha_mat .* u(tau_mat, h_mat), 'all');
    f_compact = sum(beta_mat .* f(tau_mat, h_mat), 'all');
    
    coeffs(k, 1) = u_compact == f_compact;
end

%% Compute the coeffs of compact border
coeff_1 = coeffs;
coeff_2 = coeffs;

% Introduce the normalisation
coeff_1(end-1) = alpha1_1 == 1;
coeff_1(end) = alpha1_2 == 0;
coeff_2(end-1) = alpha1_1 == 0;
coeff_2(end) = alpha1_2 == 1;

border_1 = solve(coeff_1, [alpha; beta]);
border_2 = solve(coeff_2, [alpha; beta]);
