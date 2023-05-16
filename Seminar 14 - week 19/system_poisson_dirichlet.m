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
    b = -1/5;
    c = -1/20;
    a = 1;
    r = params.r;
    p = -h^2/5 + 4*r;
    q = -h^2/40 - 2*r;
else
    error(['Scheme "', scheme, '" is not supported yet.']);
end

A = a*eye(N);
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

%% X Right border
A(N-Ny+1 : N, N-Ny+1 : N) = eye(Ny);
F(N-Ny+1 : N, :) = 0;

%% Solution
f_val = f(params.x, params.y);
rhs = F*f_val(:);
u = sparse(A) \ rhs;
u = reshape(u, Ny, Nx);
end
