%% Main parameters
C = 2;
L = 1;
T = 1;

Nx = 10;
nu = .1; % nu = C * tau / h

%% Secondary parameters
x = linspace(0, L, Nx + 1);
h = x(2) - x(1);

tau = nu * h / C;
Nt = ceil(T/tau) + 1;
T = tau*(Nt - 1);

u_0 = sin(x*2*pi/L);
u_tau = u_0;
% u_tau = u_0 * cos(C*2*pi/L*tau);

% figure(1)
% plot(x, u_0)

%% Transition matrices
U_next = (1+nu^2)*eye(Nx) + ...
    diag(-nu^2/2 * ones(Nx-1, 1), 1) + ...
    diag(-nu^2/2 * ones(Nx-1, 1),-1);
U_now = -2*eye(Nx);

%% Periodic border condition
U_now(1, end) = -nu^2;
U_now(end, 1) = -nu^2;

U_prev = U_next;

%% Integration
U = zeros(Nt, Nx);
U(1, :) = u_0(1 : end-1);
U(2, :) = u_tau(1 : end-1);

for k = 3 : Nt
    U(k, :) = -U_next \ (U_now * U(k - 1, :).' + U_prev * U(k - 2, :).');
end
U = [U, U(:, 1)];

%% Visualisation
figure(3)

for k = [1 : floor(Nt/100) : Nt, Nt]
    plot(x, U(k, :));
    t = (k-1)*tau;
    hold on;
    plot(x, sin(x/L*2*pi) * cos(C*2*pi/L*t), '--r');
    hold off;
    axis([0 L -1 1]);
    title(['t = ', num2str((k-1)*tau)]);
    drawnow;
end

%% Error
t = (0 : tau : T).';
U_true = cos(C*2*pi/L*t) * sin(x/L*2*pi);

C_norm = max(abs(U_true - U), [], 2);
L2_norm = sqrt(h*sum((U_true - U).^2, 2));

figure(4);
plot(t, C_norm);
hold on;
plot(t, L2_norm);
hold off;
legend('C', 'L^2', 'Location', 'best')
xlabel('t');
