%% Main parameters
D = 1;
L = 2*pi;
T = 6;
nu = .45; % nu = D * tau / h^2

% Nx = [12, 25, 50, 100, 200];
Nx = ceil(logspace(1, 2.1, 8));
% Nx = [12, 25, 50];

scheme_set = ["euler_ex", "cn", "compact"];
% scheme_set = ["cn", "compact"];
scheme_n = length(scheme_set);
% Orders of 3 different schemes
order_C = cell(scheme_n, 1);
order_L2 = cell(scheme_n, 1);

%% Order computation across all schemes
for i = 1 : length(Nx)
    nx = Nx(i);
    x = linspace(0, L, nx + 1);
    x = x(1:end-1);
    h = x(2) - x(1);

    tau = nu * h^2 / D;
    Nt = ceil(T/tau) + 1;
    t = linspace(0, T, Nt).';
% 
%     u_true = @(t, x) exp(-D*t)*cos(x);
%     f = @(t, x) zeros(length(t), length(x));
    u_true = @(t, x) sin(t) * cos(x);
    f = @(t, x) (D*sin(t) + cos(t))*cos(x);

    u_0 = u_true(0, x);
    u_tau = u_true(tau, x);

    for scheme = 1 : scheme_n
        switch scheme_set(scheme)
            case "euler_ex"
                % Left part
                a = 0;
                b = -1 + 2*nu;
                c = -nu;
                % Right part
                p1 = 0;
                p2 = 0;
                q1 = tau;
                q2 = 0;
            case "cn"
                alpha = 1 + nu;
                % Left part
                a = -nu/2/alpha;
                b = (-1 + nu)/alpha;
                c = -nu/2/alpha;
                % Right part
                p1 = 0;
                p2 = 0;
                q1 = tau/2/alpha;
                q2 = tau/2/alpha;
            case "compact"
                alpha = 4*(6*nu + 5);
                % Left part
                a = 2*(1 - 6*nu)/alpha;
                b = -4*(5 - 6*nu)/alpha;
                c = -2*(1 + 6*nu)/alpha;
                % Right part
                p1 = tau/alpha;
                p2 = tau/alpha;
                q1 = 10*tau/alpha;
                q2 = 10*tau/alpha;
            otherwise
                error(['Scheme ' scheme ' is not supported yet']);
        end

        U_next = eye(nx) + ...
            diag(a * ones(nx-1, 1), 1) + ...
            diag(a * ones(nx-1, 1),-1);
        U_now = b*eye(nx) + ...
            diag(c * ones(nx-1, 1), 1) + ...
            diag(c * ones(nx-1, 1),-1);

        f_next = diag(q2*ones(nx, 1)) + ...
            diag(p2*ones(nx-1, 1), 1) + diag(p2*ones(nx-1, 1), -1);
        f_now = diag(q1*ones(nx, 1)) + ...
            diag(p1*ones(nx-1, 1), 1) + diag(p1*ones(nx-1, 1), -1);

        % Periodic border condition
        U_next(1, end) = a;
        U_next(end, 1) = a;
        U_now(1, end) = c;
        U_now(end, 1) = c;
        f_next(1, end) = p2;
        f_next(end, 1) = p2;
        f_now(1, end) = p1;
        f_now(end, 1) = p1;

        % Initial conditions
        U = zeros(Nt, nx);
        U(1, :) = u_0;
        U(2, :) = u_tau;
        f_val = f(t, x);

        % Integration
        for k = 3 : Nt
            right_part = -U_now*U(k-1, :).' + ...
                f_next*f_val(k, :).' + f_now*f_val(k-1, :).';
            U(k, :) = U_next \ right_part;
        end

        c_norm = max(abs(U(end, :) - u_true(t(end), x)));
        order_C{scheme}(i) = c_norm;
    end
end

%% Visualisation
figure(1)
loglog(Nx, order_C{1}, '-o');
hold on;
for i = 2 : scheme_n
    loglog(Nx, order_C{i}, '-o');
end
hold off;

legend('Euler (explicit)', 'CN', 'Compact');
xlabel('$N_x$', Interpreter = 'latex');
title('C norm')
subtitle(['$T$ = ' num2str(t(end)) ', $\nu$ = ' num2str(nu)], Interpreter = 'latex')
