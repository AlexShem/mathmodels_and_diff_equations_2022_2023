function u = integrate_system_poisson_dirichlet(scheme, params, u0, f)
h = params.h;
N = numel(params.x);
Nx = params.Nx;
Ny = params.Ny;

if strcmp(scheme, 'standard')
    a = 1;
    b = 0;
    c = -4;
    p = h^2;
    q = 0;
    r = 0;
elseif strcmp(scheme, 'compact')
    a = -1/5;
    b = -1/20;
    c = 1;
    r = params.r;
    p = -h^2/5 + 4*r;
    q = -h^2/40 - 2*r;
else
    error(['Scheme "', scheme, '" is not supported yet.']);
end

A = c*eye(N);
F = p*eye(N);

%% X Left border
A(1:Ny, 1:Ny) = eye(Ny);
F(1:Ny, :) = 0;

%% Intermediate X steps
for k = Ny+1 : N-Ny
    if mod(k, Ny) == 1 % upper border
        A(k, k) = 1; % right
        F(k, :) = 0; % bl
    elseif mod(k, Ny) == 0 % bottom border
        A(k, k) = 1; % right
        F(k, :) = 0; % bl
    else % interior point
        A(k, k + Ny) = a; % right
        A(k, k + 1) = a; % down
        A(k, k - Ny) = a; %left
        A(k, k - 1) = a; %up
        A(k, k + Ny - 1) = b; % tr
        A(k, k + Ny + 1) = b; % br
        A(k, k - Ny + 1) = b; % bl
        A(k, k - Ny - 1) = b; % tl
        
        F(k, k + Ny) = q; % right
        F(k, k + 1) = q; % down
        F(k, k - Ny) = q; %left
        F(k, k - 1) = q; %up
        F(k, k + Ny - 1) = r; % tr
        F(k, k + Ny + 1) = r; % br
        F(k, k - Ny + 1) = r; % bl
        F(k, k - Ny - 1) = r; % tl
    end
end

%% X Right border
A(N-Ny+1 : N, N-Ny+1 : N) = eye(Ny);
F(N-Ny+1 : N, :) = 0;

%% Integration
tau = params.tau;
Nt = params.Nt;

f_val = f(params.x, params.y);
rhs = sparse(F)*f_val(:);
A_sp = sparse(A);
u = u0(:) + tau*(-A_sp*u0(:) + rhs);

for k = 3:Nt
    u0 = u;
    u = u0 + tau*(-A_sp*u0 + rhs);
end
u = reshape(u, Ny, Nx);
end
