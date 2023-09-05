function [A, F] = compact_matrix(Nx, Ny, params, scheme)
N = Nx*Ny;
hx = params.hx;

if strcmp(scheme, 'standard')
    a = 1;
    b = 0;
    c = -4;
    p = hx^2;
    q = 0;
    r = 0;
elseif strcmp(scheme, 'compact')
    a = -1/5;
    b = -1/20;
    c = 1;
    r = params.r;
    p = -hx^2/5 + 4*r;
    q = -hx^2/40 - 2*r;
else
    error(['Scheme "', scheme, '" is not supported yet.']);
end

A = sparse(c*eye(N));
F = sparse(p*eye(N));

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
end
