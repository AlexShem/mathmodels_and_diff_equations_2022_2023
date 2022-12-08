%% Main parameters
C = 1;
L = 1;
T = 1;

Nx = 121;
nu = .01; % nu = C * tau / h

%% Secondary parameters
x = linspace(0, L, Nx + 1);
h = x(2) - x(1);

tau = nu * h / C;
Nt = ceil(T/tau) + 1;

u_0 = sin(x*2*pi/L);
u_tau = u_0;
% figure(1)
% plot(x, u_0)

%% Transition matrices
U_next = eye(Nx);
U_now = (2-2*nu^2)*eye(Nx) + ...
    diag(nu^2 * ones(Nx-1, 1), 1) + ...
    diag(nu^2 * ones(Nx-1, 1),-1);
U_prev = -U_next;

%% Periodic border condition
U_now(1, end) = nu^2;
U_now(end, 1) = nu^2;

%% Integration
U = zeros(Nt, Nx);
U(1, :) = u_0(1 : end-1);
U(2, :) = u_tau(1 : end-1);

for k = 3 : Nt
    U(k, :) = U_next \ (U_now * U(k - 1, :).' + U_prev * U(k - 2, :).');
end
U = [U, U(:, 1)];

%% Visualisation
figure(2)

for k = [1 : floor(Nt/100) : Nt, Nt]
    plot(x, U(k, :));
    axis([0 L -1 1]);
    title(['t = ', num2str((k-1)*tau)]);
    drawnow;
end
