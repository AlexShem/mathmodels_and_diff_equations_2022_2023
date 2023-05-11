%% Primary parameters
rho = 7900;
R = 10^-2;
E = 210e9;
L = 2*pi;
T = 1;

Nx = 100;
x = linspace(0, L, Nx + 1);
h = x(2) - x(1);
nu = 0.1;

D = R^2;
C = E*R^2/rho;
mu = D/h^2;

scheme = "CN";
% scheme = "555";

%% Secondary parameters

%nu = C tau^2 / h^4
tau = sqrt(nu*h^4 / C);
Nt = ceil(T/tau) + 1;
T_fin = tau * (Nt-1);
t = linspace(0, T_fin, Nt);

%% Reference soltuion
u_ref = @(t, x) sin(x + t*sqrt(C/(D+1)));
f_ref = @(t, x) zeros(size(x));

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
%         alpha = 1 + 3*nu + 2*mu;
%         a = -(2*nu + mu)/alpha;
%         b = -(2 + 4*mu)/alpha;
%         c = 2*mu/alpha;
%         d = 0;
%         e = nu/2/alpha;
%         p0 = 0;
%         p1 = 0;
%         q0 = 0;
%         q1 = 0;
%         r0 = 0;
%         r1 = tau^2/alpha;
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

% Periodic boundary conditions
U_next(1, end-1:end) = [e a];
U_next(2, end) = e;
U_next(end-1, 1) = e;
U_next(end, 1:2) = [a e];

U_now(1, end-1:end) = [d c];
U_now(2, end) = d;
U_now(end-1, 1) = d;
U_now(end, 1:2) = [c d];

F_next(1, end-1:end) = [p1 q1];
F_next(2, end) = p1;
F_next(end-1, 1) = p1;
F_next(end, 1:2) = [q1 p1];

F_now(1, end-1:end) = [p0 q0];
F_now(2, end) = p0;
F_now(end-1, 1) = p0;
F_now(end, 1:2) = [q0 p0];

U_prev = U_next;
F_prev = F_next;

%% Integration
u = zeros(Nt, Nx);
u(1, :) = u_ref(0, x(1:end-1));
u(2, :) = u_ref(tau, x(1:end-1));
[x_gr, t_gr] = meshgrid(x(1:end-1), t);

f = f_ref(t_gr, x_gr);

for k = 3 : Nt
    rhs = -U_now*u(k-1, :).' - U_prev*u(k-2, :).' + ...
        F_next*f(k, :).' + F_now*f(k-1, :).' + F_prev*f(k-2, :).';
    u(k, :) = U_next \ rhs;
end
u = [u, u(:, 1)];

%% Visualisation
figure(2)
for k = 1 : 100 : Nt
    plot(x, u(k, :));
    axis([0 L -1 1]);
    pause(0.1);
end
