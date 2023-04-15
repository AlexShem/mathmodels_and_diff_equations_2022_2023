%% Main parameters
c = 1;
L = 2*pi;
T = .5;
nu = .5; % nu = c*tau/h^2;
Nx_set = [25 50 100 200];
Nx = Nx_set(3);

x = linspace(0, L, Nx + 1);
h = x(2) - x(1);
tau = nu*h^2/c;
Nt = ceil(T/tau) + 1;
T = tau*(Nt - 1);
t = 0 : tau : T;

scheme = 'compact';
% scheme = 'explicit';
% comp_correction = false;
comp_correction = true;

anim = true;
% anim = false;

%% Reference solution
A = 1;
B = 2;
lam = 1;

% u_ref = @(t, x, A, B, lam) lam/c * (lam*t + x) + A;
u_ref = @(t, x, A, B, lam) (x - A).^2./(6*c*(B - t));

u_0 = @(x) u_ref(0, x, A, B, lam);
f = @(u) u.^2/2;

%% Integration
U = zeros(Nt, Nx);
u_0_val = u_0(x(1:end-1));
U(1, :) = u_0_val;

for k = 2 : Nt
    u = U(k-1, :);
    u_ex = explicit_euler(u, f, nu);
    %     u_ex = lax_wendroff(u, f, tau, h);
    %     u_ex = maccormack(u, f, tau, h);

    if comp_correction
        comp_eps = compact_correction(u, u_ex, c, h, tau);
    else
        comp_eps = zeros(length(u), 1);
    end
    u_k = u_ex + comp_eps.';
    U(k, :) = u_k;
end
U = [U, U(:, 1)];

%% Visualisation
ax = [0 L min(U, [], 'all') max(U, [], 'all')];

if anim
    figure(1)
    for k = [1 : floor(Nt/200) : Nt, Nt]
        plot(x, U(k, :));
        axis(ax);
        title(['t = ', num2str((k-1)*tau)]);
        drawnow;
    end
end
