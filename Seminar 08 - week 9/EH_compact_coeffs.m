clear

%% Initialization
syms a b c p q r [1 2]
syms h tau positive % Positive parameters
syms x t u(t, x) f(t, x)

coefs_u_app = [a b c];
coefs_f_app = [p q r];

%% Test functions
u_test = [sym(1), sym(0), sym(0), ...
    t, t^2, x, ...
    x^2, t*x, t*x^2, ...
    x*t^2, x^2*t^2].';
f_test = [sym(0), sym(1), t, ...
    x, 2*t*x, 0, ...
    0, x^2/2, x^3/3, ...
    t*x^2, 2*t*x^3/3].';

%% Differential equation
% du_fun = @(U) diff(U, t) - diff(f, x);

%% Template
u_compact_scheme = [a; b; c].';
f_compact_scheme = [p; q; r].';

%% System of equations
coef_eqs = sym('Eqs', [length(u_test) + 1, 1]);
[x_mesh, t_mesh] = meshgrid([-h, 0, h], [0, tau]);

for k = 1 : numel(u_test)
    u(t, x) = u_test(k);
    f(t, x) = f_test(k);

    u_comact = u(t_mesh, x_mesh) .* u_compact_scheme;
    f_comact = f(t_mesh, x_mesh) .* f_compact_scheme;

    coef_eqs(k, 1) = sum(u_comact, 'all') == sum(f_comact, 'all');
end
coef_eqs(end) = p2 == -3*tau/(2*h);
coef_eqs(end+1) = p2 == -r2;

%% Solution
comp_eqs = solve(coef_eqs, [coefs_u_app, coefs_f_app])
