%% Initialization
clear
clc
syms a b c % Coeffs of the left side (approximating u)
syms 'p' 'q' [1 2] % Coeffs of the right side (approximating f)
syms D h tau nu positive % Positive parameters
syms x t u(t, x) f(t, x)

coefs_u_app = [a b c];
coefs_f_app = [p q];

%% Test functions
x_test = x.^[0 2 4];
t_test = t.^[0 1 2];
test_funs = transpose(x_test.' * t_test);
test_funs = test_funs(:);
test_funs = test_funs(1:end-2); % Keep 7 test functions

%% Differential equation
du_fun = @(U) diff(U, t) - D * diff(U, x, 2);

%% System of equations
coef_eqs = sym('Eqs', [length(test_funs), 1]);
for k = 1 : numel(test_funs)
    u(t, x) = test_funs(k);
    f(t, x) = du_fun(u(t, x));
    
    u_comact = 1*u(t+tau, x) + ...
        a * (u(t+tau, x+h) + u(t+tau, x-h)) + ...
        b * u(t, x) + ...
        c * (u(t, x+h) + u(t, x-h));
    u_comact = subs(u_comact, [t, x], [0, 0]);
    f_comact = q1 * f(t, x) + ...
        q2 * f(t+tau, x) + ...
        p1 * (f(t, x+h) + f(t, x-h)) + ...
        p2 * (f(t+tau, x+h) + f(t+tau, x-h));
    f_comact = subs(f_comact, [t, x], [0, 0]);
    
    coef_eqs(k, 1) = u_comact == f_comact;
end

comp_eqs = solve(coef_eqs, [coefs_u_app, coefs_f_app]);

%% Solution analysis
cfs = struct2array(comp_eqs).';
cfs = simplify(subs(cfs, D*tau, nu*h^2));
% nu * h^2 = D * tau
alpha = 4*(6*nu + 5);
cfs_norm = simplify(cfs * alpha);
