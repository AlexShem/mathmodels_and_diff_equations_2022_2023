clear

%% Initialization
syms a b c p q r [1 2]
syms D h tau nu positive % Positive parameters
syms x t u(t, x) f(t, x)

coefs_u_app = [a b c];
coefs_f_app = [p q r];

%%Test functions
u_test = [sym(1), sym(0), sym(0), ...
    t, t^2, x, ...
    x^2, x^3, t*x, ...
    t*x^2, x*t^2].';
f_test = [sym(0), sym(1), t, ...
    x, 2*t*x, 0, ...
    -2*D*x, -3*D*x^2, x^2/2, ...
    x^3/3-2*D*t*x, t*x^2].';

%%Differential equation
du_fun = @(U) diff(U, t) - D*diff(U, x, 2) - diff(f, x);

%%Template
u_compact_scheme = [a; b; c].';
f_compact_scheme = [p; q; r].';

%%System of equations
coef_eqs = sym('Eqs', [length(u_test) + 1, 1]);
[x_mesh, t_mesh] = meshgrid([-h, 0, h], [0, tau]);

for k = 1 : numel(u_test)
    u(t, x) = u_test(k);
    f(t, x) = f_test(k);

    u_comact = u(t_mesh, x_mesh) .* u_compact_scheme;
    f_comact = f(t_mesh, x_mesh) .* f_compact_scheme;

    coef_eqs(k, 1) = sum(u_comact, 'all') == sum(f_comact, 'all');
end
coef_eqs(end) = b2 == -8/3 - 4*D*tau/h^2;
% coef_eqs(end) = b2 == 4 + 6*D*tau/h^2;


[A, b] = equationsToMatrix(coef_eqs, [coefs_u_app, coefs_f_app]);
% [(1:1:length(coef_eqs)).', [u_test; 'norm'], [f_test; 'cond'], coef_eqs]
disp(['Rank of the system: ', num2str(rank(A))]);
A

comp_eqs = solve(coef_eqs, [coefs_u_app, coefs_f_app])
comp_eqs2 = structfun(@(f) simplify(subs(f, D*tau, nu*h^2)), comp_eqs)
