%% Main parameters
D = .1;
L = 1;
T = 1;

Nx = 9;
nu = 1.5; % nu = D * tau / h^2

%% Secondary parameters
x = linspace(0, L, Nx + 1);
h = x(2) - x(1);

tau = nu * h^2 / D;
Nt = ceil(T/tau) + 1;

u_0 = exp(-(x - .5).^2 ./ .05);
% figure(1)
% plot(x, u_0)

%% Transition matrices
U_next = (-2/nu-2)*eye(Nx) + ...
    diag(ones(Nx-1, 1), 1) + ...
    diag(ones(Nx-1, 1),-1);
U_now = (2/nu-2)*eye(Nx) + ...
    diag(ones(Nx-1, 1), 1) + ...
    diag(ones(Nx-1, 1),-1);

%% Periodic border condition
U_next(1, end) = 1;
U_next(end, 1) = 1;
U_now(1, end) = 1;
U_now(end, 1) = 1;

%% Integration
U = zeros(Nt, Nx);
U(1, :) = u_0(1 : end-1);

for k = 2 : Nt
    U(k, :) = -U_next \ U_now * U(k - 1, :).';
end
U = [U, U(:, 1)];

%% Visualisation
figure(2)

for k = [1 : floor(Nt/100) : Nt, Nt]
    plot(x, U(k, :));
    axis([0 L 0 1]);
    title(['t = ', num2str((k-1)*tau)]);
    drawnow;
end
