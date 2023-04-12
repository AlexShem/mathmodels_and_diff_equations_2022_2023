%% Main parameters
D = .5;
L = 2*pi;
T = .1;
nu = .1; % nu = D*tau/h^2;
Nx = [25 50 100 200];

correction = 'explicit';
% correction = 'compact';

u_ref = @(t, x) -2*D*cos(x) ./ (2*exp(D*t) + sin(x));
u_0 = @(x) u_ref(0, x);
f = @(u) -u.^2/2;

%%
C_norm = zeros(size(Nx));
figure(1);

for j = 1 : length(Nx)
    x = linspace(0, L, Nx(j) + 1);
    h = x(2) - x(1);
    tau = nu * h.^2 / D;

    [U, ~, ~, x, ~, t] = B_integration(L, T, D, Nx(j), tau, nu, u_0, correction);
    [x_mesh, t_mesh] = meshgrid(x, t);
    Cn = max(abs(u_ref(t_mesh, x_mesh) - U), [], 2);
    semilogy(t, Cn, LineWidth = 1); hold on;

    C_norm(j) = Cn(end);
end
hold off;
xlabel('$t$, time', Interpreter = 'latex', FontSize = 14)
legend(num2str(Nx.'), FontSize = 14, Location = 'southeast')

%% Visualisation: Order
figure(2);
loglog(Nx, C_norm, '-o', LineWidth = 1);
xlabel('$N_x$', Interpreter = 'latex', FontSize = 14);
title('$C$ norm', Interpreter = 'latex', FontSize = 16);

%% Order
order_comp = (log10(C_norm(end)) - log10(C_norm(1))) / ...
    (log10(Nx(end)) - log10(Nx(1)));
disp(['[ORDER] ', correction, ' scheme: ', num2str(-order_comp), '. nu = ', num2str(nu)])
