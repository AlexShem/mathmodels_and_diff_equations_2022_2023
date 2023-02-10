%% Main parameters
D = .1;
L = 1;
T = 1;

Nx = 100;
nu = 1; % nu = D * tau / h^2

%% Secondary parameters
x = linspace(0, L, Nx + 1);
h = x(2) - x(1);

tau = nu * h^2 / D;
Nt = ceil(T/tau) + 1;

u_0 = exp(-(x - .5).^2 ./ .05);
u_0 = sin(2*pi*x/L);

% figure(1)
% plot(x, u_0)

%% Transition matrices
alpha = 4*(6*nu + 5);
a = 2*(1 - 6*nu)/alpha;
b = -4*(5 - 6*nu)/alpha;
c = -2*(1 + 6*nu)/alpha;

U_now = b*eye(Nx) + ...
    diag(c * ones(Nx-1, 1), 1) + ...
    diag(c * ones(Nx-1, 1),-1);
U_next = eye(Nx) + ...
    diag(a * ones(Nx-1, 1), 1) + ...
    diag(a * ones(Nx-1, 1),-1);

%% Periodic border condition
U_next(1, end) = a;
U_next(end, 1) = a;
U_now(1, end) = c;
U_now(end, 1) = c;

%% Integration
U = zeros(Nt, Nx);
U(1, :) = u_0(1 : end-1);

for k = 2 : Nt
    U(k, :) = -U_next \ U_now * U(k - 1, :).';
end
U = [U, U(:, 1)];

%% Visualisation
figure(2)

for k = [1 : floor(Nt/200) : Nt, Nt]
    plot(x, U(k, :));
    axis([0 L -1 1]);
    title(['t = ', num2str((k-1)*tau)]);
    drawnow;
end
