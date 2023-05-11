%% Primary parameters
rho = 7900;
R = 10^-2;
E = 210e9;
L = 2*pi;
T = 1;

Nx = 101;
x = linspace(0, L, Nx);
h = x(2) - x(1);
nu = 0.1;

D = R^2;
C = E*R^2/rho;
mu = D/h^2;

% scheme = "CN";
% scheme = "353";
scheme = "555";

% border = "naive";
border = "compact";

%% Secondary parameters

%nu = C tau^2 / h^4
tau = sqrt(nu*h^4 / C);
Nt = ceil(T/tau) + 1;
T_fin = tau * (Nt-1);
t = linspace(0, T_fin, Nt);

%% Reference soltuion
u_ref = @(t, x) cos(2*pi*t).*sin(x).^2;
f_ref = @(t, x) 4*cos(2*pi*t).*(pi^2*cos(2*x)/2 - pi^2/2 - 2*C*cos(2*x) + 2*D*pi^2*cos(2*x));

%% Computational scheme
switch scheme
    case "CN"
        alpha = 1 + 3*nu + 2*mu;
        a = -(2*nu + mu)/alpha;
        b = -(2 + 4*mu)/alpha;
        c = 2*mu/alpha;
        d = 0;
        e = nu/2/alpha;
        p0 = 0;
        p1 = 0;
        q0 = 0;
        q1 = 0;
        r0 = 0;
        r1 = tau^2/alpha;
    case "555"
        alpha = 72*(41 + 30*nu + 90*mu);
        a = 96*(7 - 15*nu - 30*mu) / alpha;
        b = -144*(41 - 150*nu + 90*mu) / alpha;
        c = -192*(7 + 75*nu - 30*mu) / alpha;
        d = -24*(1 - 150*nu - 30*mu) / alpha;
        e = 12*(1 + 30*nu - 30*mu) / alpha;
        p0 = 10*tau^2 / alpha;
        q0 = 560*tau^2 / alpha;
        r0 = 2460*tau^2 / alpha;
        p1 = tau^2 / alpha;
        q1 = 56*tau^2 / alpha;
        r1 = 246*tau^2 / alpha;
    case "353"
        alpha = 96*(1 + 3*mu);
        a = 24*(1 - 6*mu) / alpha;
        b = -96*(2 - 9*nu + 6*mu) / alpha;
        c = -48*(1 + 12*nu - 6*mu) / alpha;
        d = 144*nu / alpha;
        e = 0;
        p0 = -10*tau^2*(nu - mu) / alpha;
        q0 = 20*tau^2*(1 + 2*nu - 2*mu) / alpha;
        r0 = 20*tau^2*(4 - 3*nu + 3*mu) / alpha;
        p1 = -tau^2*(nu - mu) / alpha;
        q1 = 2*tau^2*(1 + 2*nu - 2*mu) / alpha;
        r1 = 2*tau^2*(4 - 3*nu + 3*mu) / alpha;
    otherwise
        error("This scheme is not supported yet.")
end

U_next = eye(Nx) + diag(a*ones(Nx-1, 1), 1) + diag(a*ones(Nx-1, 1), -1) + ...
    diag(e*ones(Nx-2, 1), 2) + diag(e*ones(Nx-2, 1), -2);
U_now = b*eye(Nx) + diag(c*ones(Nx-1, 1), 1) + diag(c*ones(Nx-1, 1), -1) + ...
    diag(d*ones(Nx-2, 1), 2) + diag(d*ones(Nx-2, 1), -2);

F_next = r1*eye(Nx) + diag(q1*ones(Nx-1, 1), 1) + diag(q1*ones(Nx-1, 1), -1) + ...
    diag(p1*ones(Nx-2, 1), 2) + diag(p1*ones(Nx-2, 1), -2);
F_now = r0*eye(Nx) + diag(q0*ones(Nx-1, 1), 1) + diag(q0*ones(Nx-1, 1), -1) + ...
    diag(p0*ones(Nx-2, 1), 2) + diag(p0*ones(Nx-2, 1), -2);

%% Dirichlet boundary conditions
switch border
    case "naive"
        U_next(1, 1:4) = [1 0 0 0];
        U_next(2, 1:4) = [0 1 0 0];
        U_next(end-1, end-3:end) = [0 0 1 0];
        U_next(end, end-3:end) = [0 0 0 1];

        U_now(1, 1:4) = 0;
        U_now(2, 1:4) = 0;
        U_now(end-1, end-3:end) = 0;
        U_now(end, end-3:end) = 0;

        F_next(1, 1:4) = 0;
        F_next(2, 1:4) = 0;
        F_next(end-1, end-3:end) = 0;
        F_next(end, end-3:end) = 0;


        F_now(1, 1:4) = 0;
        F_now(2, 1:4) = 0;
        F_now(end-1, end-3:end) = 0;
        F_now(end, end-3:end) = 0;

        U_prev = U_next;
        F_prev = F_next;

        U_prev([1 2 end-1 end], :) = 0;
    case "compact"
        U_next(1, 1:4) = [1 0 0 0];
        U_next(2, 1:4) = [0 1 -.5 1/9];
        U_next(end-1, end-3:end) = [1/9 -.5 1 0];
        U_next(end, end-3:end) = [0 0 0 1];

        U_now(1, 1:4) = 0;
        U_now(2, 1:4) = 0;
        U_now(end-1, end-3:end) = 0;
        U_now(end, end-3:end) = 0;

        F_next(1, 1:4) = 0;
        F_next(2, 1:4) = [-(h^2*(- h^2 + 2*D))/(12*C), D*h^2/(6*C), 0, 0];
        F_next(end-1, end-3:end) = [0, 0, D*h^2/(6*C), -(h^2*(- h^2 + 2*D))/(12*C)];
        F_next(end, end-3:end) = 0;


        F_now(1, 1:4) = 0;
        F_now(2, 1:4) = 0;
        F_now(end-1, end-3:end) = 0;
        F_now(end, end-3:end) = 0;

        U_prev = U_next;
        F_prev = F_next;

        U_prev([1 2 end-1 end], :) = 0;
        F_prev([1 2 end-1 end], :) = 0;
    otherwise
        error("This border approximation is not supported yet.")
end

%% Integration
u = zeros(Nt, Nx);
u(1, :) = u_ref(0, x);
u(2, :) = u_ref(tau, x);
[x_gr, t_gr] = meshgrid(x, t);

f = f_ref(t_gr, x_gr);

for k = 3 : Nt
    rhs = -U_now*u(k-1, :).' - U_prev*u(k-2, :).' + ...
        F_next*f(k, :).' + F_now*f(k-1, :).' + F_prev*f(k-2, :).';
    u(k, :) = U_next \ rhs;
end

%% Visualisation
fr_sh = .005;
fn = ceil(Nt * fr_sh);
frames = [floor(linspace(1, Nt, fn)), Nt];

figure(2)
for k = frames
    cur_t = (k-1)*tau;
    if k == 1
        p1 = plot(x, u(k, :));
        hold on;
        p2 = plot(x, u_ref(cur_t, x), '--r');
        hold off;
        xlabel("x");
        ylabel("u");
        ttl = title("t = " + num2str(round(cur_t, 2)));
        axis([0 L -1 1]);
        legend("$u\,(t, x)$", "$u^*(t, x)$", interpreter = "latex", location = "northeast");
        drawnow;
    else
        p1.YDataSource = 'y1';
        p2.YDataSource = 'y2';
        y1 = u(k, :);
        y2 = u_ref(cur_t, x);
        refreshdata;
        ttl.String = "t = " + num2str(round(cur_t, 2));
        drawnow;
    end
end

%% Error
u_star = u_ref(t_gr, x_gr);
C_norm = max(abs(u - u_star), [], 2);

figure(3)
semilogy(t, C_norm);
hold on;
xlabel("$t$", Interpreter="latex");
title("$C$-norm", Interpreter="latex")
