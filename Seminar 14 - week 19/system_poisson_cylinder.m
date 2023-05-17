function u = system_poisson_dirichlet(scheme, params, f)
h = params.h;
N = numel(params.x);
Nx = params.Nx;
Ny = params.Ny;

if strcmp(scheme, 'standard')
    a = -4;
    b = 1;
    c = 0;
    p = h^2;
    q = 0;
    r = 0;
elseif strcmp(scheme, 'compact')
    a = 1;
    b = -1/5;
    c = -1/20;
    r = params.r;
    p = -h^2/5 + 4*r;
    q = -h^2/40 - 2*r;
else
    error(['Scheme "', scheme, '" is not supported yet.']);
end

% Right border - periodic conditions
A = a*eye(N - Ny);
F = p*eye(N - Ny);

%% Top Bottom Dirichlet border
A(1:Ny:end, 1:Ny:end) = eye(Ny-1);
F(1:Ny:end, :) = 0;

A(Ny:Ny:end, Ny:Ny:end) = eye(Ny-1);
F(Ny:Ny:end, :) = 0;

%% Left periodic border
for k = 2 : Ny-1
    A(k, k + Ny) = b; % right
    A(k, k + 1) = b; % down
    A(k, N - 2*Ny + k) = b; %left
    A(k, k - 1) = b; %up
    A(k, k + Ny - 1) = c; % tr
    A(k, k + Ny + 1) = c; % br
    A(k, N - 2*Ny + k + 1) = c; % bl
    A(k, N - 2*Ny + k - 1) = c; % tl

    F(k, k + Ny) = q; % right
    F(k, k + 1) = q; % down
    F(k, N - 2*Ny + k) = q; %left
    F(k, k - 1) = q; %up
    F(k, k + Ny - 1) = r; % tr
    F(k, k + Ny + 1) = r; % br
    F(k, N - 2*Ny + k + 1) = r; % bl
    F(k, N - 2*Ny + k - 1) = r; % tl
    F(k, N - 2*Ny + k + 1) = r; % bl
    F(k, N - 2*Ny + k - 1) = r; % tl
end

%% Intermediate X steps
for k = Ny+1 : N-2*Ny
    if mod(k, Ny) == 1 % upper border
        continue;
    elseif mod(k, Ny) == 0 % bottom border
        continue;
    else % interior point
        A(k, k + Ny) = b; % right
        A(k, k + 1) = b; % down
        A(k, k - Ny) = b; %left
        A(k, k - 1) = b; %up
        A(k, k + Ny - 1) = c; % tr
        A(k, k + Ny + 1) = c; % br
        A(k, k - Ny + 1) = c; % bl
        A(k, k - Ny - 1) = c; % tl

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

%% Right periodic border
for k = N-2*Ny+2 : N-Ny-1
    kl = mod(k, Ny);
    A(k, kl) = b; % right
    A(k, k + 1) = b; % down
    A(k, k - Ny) = b; %left
    A(k, k - 1) = b; %up
    A(k, kl - 1) = c; % tr
    A(k, kl + 1) = c; % br
    A(k, k - Ny + 1) = c; % bl
    A(k, k - Ny - 1) = c; % tl

    F(k, kl) = q; % right
    F(k, k + 1) = q; % down
    F(k, k - Ny) = q; %left
    F(k, k - 1) = q; %up
    F(k, kl - 1) = r; % tr
    F(k, kl + 1) = r; % br
    F(k, k - Ny + 1) = r; % bl
    F(k, k - Ny - 1) = r; % tl
end

%% Solution
f_val = f(params.x(:, 1:end-1), params.y(:, 1:end-1));
rhs = F*f_val(:);
u = sparse(A) \ rhs;

u = reshape(u, Ny, Nx-1);
u = [u, u(:, 1)];
end
