%% Parameters
C = 1;
nu = .1;
Nx = 90;
T = 1;
L = 1;

%% Secondary praramters
x = linspace(0, L, Nx);
% L / (Nx - 1)
h = x(2) - x(1);

% nu = C * tau / h
tau = nu * h / C;
Nt = ceil(T / tau) + 1;
T_fin = (Nt - 1) * tau;
t = linspace(0, T_fin, Nt);

u0 = cos(2*pi*x/L);
utau = u0;

%% Transition matrices
U_next = (1 + nu^2)*eye(Nx) + ...
    diag(-nu^2/2*ones(Nx-1, 1), 1) + ...
    diag(-nu^2/2*ones(Nx-1, 1), -1);
U_now = -2*eye(Nx);
U_prev = U_next;

%% Border conditions
U_now([1, end], :) = 0;
U_next([1, end], :) = 0; U_next(1, 1) = 1; U_next(end, end) = 1;
U_next(1, 2) = -1; U_next(end, end-1) = -1;
U_prev([1, end], :) = 0;

%% Integration
u = zeros(Nt, Nx);
u(1, :) = u0;
u(2, :) = utau;

for k = 3 : Nt
    rhs = -U_now*(u(k-1, :).') - U_prev*(u(k-2, :).');
    rhs([1, end]) = 0;
    u(k, :) = U_next \ rhs;
end

%% Visualisation
[xm, tm] = meshgrid(x, t);

figure(2)
surf(xm, tm, u);
xlabel('x');
ylabel('t');

%% Animation
figure(3)
for k = 1 : 1 : Nt
    plot(x, u(k, :));
    axis([0 1 -1 1]);
    t = tau*(k-1);
    title(['t = ' num2str(t)]);
    drawnow;
end
