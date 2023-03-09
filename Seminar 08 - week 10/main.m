clear

%% Main parameters
L = 2*pi;
T = 1.5;
tau = 0.0061359;
Nx = 101;

comp_correction = true;

%% Secondary parameters
x = linspace(0, L, Nx + 1);
h = x(2) - x(1);

Nt = ceil(T/tau) + 1;
T = tau*(Nt - 1);
t = 0 : tau : T;

% u_0 = x.*sin(x).^3; ax = [0 L -6 2];
u_0 = sin(x); ax = [0 L -2 2];
f = @(u) u.^2/2;
% figure(1)
% plot(x, u_0)

%% Integration
U = zeros(Nt, Nx);
U(1, :) = u_0(1 : end-1);

for k = 2 : Nt
    u = U(k-1, :);
    u_ex = explicit_euler(u, f, tau, h);
%     u_ex = lax_wendroff(u, f, tau, h);
%     u_ex = maccormack(u, f, tau, h);

    if comp_correction
        comp_eps = compact_correction(u, u_ex, h, tau);
    else
        comp_eps = zeros(length(u), 1);
    end
    U(k, :) = u_ex + comp_eps.';
end
U = [U, U(:, 1)];

%% Visualisation
figure(2)
for k = [1 : floor(Nt/200) : Nt, Nt]
    plot(x, U(k, :));
    axis(ax);
    title(['t = ', num2str((k-1)*tau)]);
    drawnow;
end
