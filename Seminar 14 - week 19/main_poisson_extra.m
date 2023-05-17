%% Parameters

% scheme = 'compact'; params.r = 0.0;
scheme = 'standard';

L = 2*pi;

Nx = 101;
Ny = Nx;

h = L / Nx;

x = linspace(0, L, Nx);
y = linspace(0, L, Ny);
[x, y] = meshgrid(x, y);

T = 1;
Nt = 1000;
tau = T/(Nt - 1);
u0 = zeros(size(x));

params.x = x;
params.y = y;
params.h = h;
params.Nx = Nx;
params.Ny = Ny;
params.T = T;
params.tau = tau;
params.Nt = Nt;

%% Right-hand side
syms xs ys
Lu = @(u) diff(u, xs, 2) + diff(u, ys, 2);

u_sym = sin(xs)*sin(ys);
% u_sym = exp(-(xs-L/2)^2 - (ys-L/2)^2);

f_sym = Lu(u_sym);
u_true = matlabFunction(u_sym, 'Vars', [xs, ys]);
f = matlabFunction(f_sym, 'Vars', [xs, ys]);
clear xs ys u_sym f_sym

%% Integration
tol = 1e-3;
% rng(1);

u = system_poisson_dirichlet(scheme, params, f);
% u = integrate_system_poisson_dirichlet(scheme, params, u0, f);
% [u, Ni] = seidel_system_poisson_dirichlet(scheme, params, u0, f, tol, 1);

%% Visualisation
figure(1)
% contour(x, y, u, 'ShowText', 'on');
surf(x, y, u);
xl = xlabel('$x$'); xl.Interpreter = 'latex'; xl.FontSize = 16;
yl = ylabel('$y$'); yl.Interpreter = 'latex'; yl.FontSize = 16;

%% Visualisation 2
figure(2)
contour(x, y, log10(abs(u - u_true(x, y))), 18, 'ShowText', 'on');
xl = xlabel('$x$'); xl.Interpreter = 'latex'; xl.FontSize = 16;
yl = ylabel('$y$'); yl.Interpreter = 'latex'; yl.FontSize = 16;
ttl = title(['$\log_{10} |u - u^*|$']); ttl.Interpreter = 'latex'; ttl.FontSize = 16;
