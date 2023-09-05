%% Parameters
scheme = 'compact'; params.r = 0.1;
% scheme = 'standard';

L = 2*pi;

Nx = 11;
Ny = Nx;

hx = L / Nx;
hy = L / Ny;

x = linspace(0, L, Nx);
y = linspace(0, L, Ny);
[x, y] = meshgrid(x, y);

u0 = zeros(size(x));

params.x = x;
params.y = y;
params.hx = hx;
params.hy = hy; 
params.Nx = Nx;
params.Ny = Ny;
params.L = L;

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
mg_lvls = 4;
N_sim = 10;
% N_iter = NaN(N_sim, 1);

tol = 1e-4;
tic
for k = 1 : N_sim
    [u_sol, Ni] = seidel_system_poisson_dirichlet(scheme, params, u0, f, tol, 1, mg_lvls);
    N_iter(k) = Ni;
end
toc

N_avg = mean(N_iter);
[~, N_const] = seidel_system_poisson_dirichlet(scheme, params, u0, f, tol, 0);

disp(['N simulations: ', num2str(N_sim)]);
disp(['Average N: ', num2str(N_avg)]);
disp(['Constant N: ', num2str(N_const)]);
disp('-------------------------');

%% Visualisation
figure(1);
plot(1:N_sim, N_iter, 'o-');
line([1 N_sim], [N_const N_const], 'Color', 'red');
line([1 N_sim], [N_avg N_avg], 'Color', 'black', 'LineStyle', '--');
xlabel('Number of simulations');
ylabel('Number of Seidel iterations');
