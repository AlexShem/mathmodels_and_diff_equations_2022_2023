function [u, Niter] = seidel_system_poisson_dirichlet(scheme, params, u0, f, tol, method, mg_lvls)
if nargin <= 6
    mg_lvls = 1;
end

hx = params.hx;
% hy = params.hy;
N = numel(params.x);
Nx = params.Nx;
Ny = params.Ny;
Niter = 0;


%% Integration
Niter_vec = 3.^((mg_lvls) : -1 : 1);
dir = 1; % == 1 -- direction up (to the finer grid)

u0 = u0(:);
u = Inf(size(u0));

err = norm(u - u0, 2);

while err > tol
    u = u0(:);
    dir = 1;
    %     flag = false; % Forward-backward interation done?
    
    %% Forward grid
    for lvl = 1 : mg_lvls
        [A, F] = compact_matrix(Nx, Ny, params, scheme);
        x = linspace(0, params.L, Nx);
        y = linspace(0, params.L, Ny);
        [x, y] = meshgrid(x, y);
        f_val = f(x, y);
        rhs = F*f_val(:);
        
        for j = 1 : Niter_vec(lvl)
            if method == 1
                ind_ord = randperm(N);
            else
                ind_ord = 1:N;
            end
            u = u(:);
            for k  = 1 : N
                i = ind_ord(k);
                r_new = ind_ord(1 : (k-1));
                r_old = ind_ord((k+1) : end);
                
                s1 = A(i, r_new) * u(r_new);
                s2 = A(i, r_old) * u(r_old);
                u(i) = (rhs(i) - s1 - s2)/A(i, i);
            end
        end
        
        % Interpolation
        u_mat = reshape(u, Ny, Nx);
        u = NaN(Ny, 2*Nx - 1);
        for row = 1 : Ny
            u(row, :) = interp_compact(u_mat(row, :));
        end
        u_mat = u;
        u = NaN(2*Ny - 1, 2*Nx - 1);
        for col = 1 : (2*Nx-1)
            u(:, col) = interp_compact(u_mat(:, col));
        end
        [Ny, Nx] = size(u);
        N = Ny*Nx;
    end
    %% Backward grid
    for lvl = mg_lvls-1 : -1 : 1
        % Matrix reduction
        u = reshape(u, Ny, Nx);
        u = u(1:2:end, 1:2:end);
        [Ny, Nx] = size(u);
        N = Ny*Nx;
        
        [A, F] = compact_matrix(Nx, Ny, params, scheme);
        x = linspace(0, params.L, Nx);
        y = linspace(0, params.L, Ny);
        [x, y] = meshgrid(x, y);
        f_val = f(x, y);
        rhs = F*f_val(:);
        
        for j = 1 : Niter_vec(lvl)
            if method == 1
                ind_ord = randperm(N);
            else
                ind_ord = 1:N;
            end
            u = u(:);
            for k  = 1 : N
                i = ind_ord(k);
                r_new = ind_ord(1 : (k-1));
                r_old = ind_ord((k+1) : end);
                
                s1 = A(i, r_new) * u(r_new);
                s2 = A(i, r_old) * u(r_old);
                u(i) = (rhs(i) - s1 - s2)/A(i, i);
            end
        end
    end
    
    u = reshape(u, Ny, Nx);
    u = u(1:2:end, 1:2:end);
    [Ny, Nx] = size(u);
    N = Ny*Nx;
    
    err = norm(u(:) - u0(:), 2);
    u0 = u;
    %     Niter = Niter + 1;
    %     if Niter > 500
    %         error('Seidel method did not converge');
    %     end
    
    Niter = Niter + sum(Niter_vec);
end

u = reshape(u, Ny, Nx);
end
