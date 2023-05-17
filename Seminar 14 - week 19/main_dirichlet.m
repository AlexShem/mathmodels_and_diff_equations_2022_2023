%% Parameters

scheme = 'compact'; params.r = 0;
% scheme = 'standard';

L = 2*pi;
Nx_set = 2.^(3:7);
C_norm = nan(length(Nx_set), 1);

show_vis = false;

%% Right-hand side
syms xs ys
Lu = @(u) diff(u, xs, 2) + diff(u, ys, 2);

u_sym = sin(xs)*sin(ys);
% u_sym = sin(xs)*;
% u_sym = sin(xs)*(-(ys - sym(pi))^2 + sym(pi)^2);

f_sym = Lu(u_sym);
u_true = matlabFunction(simplify(u_sym), 'Vars', [xs, ys]);
f = matlabFunction(simplify(f_sym), 'Vars', [xs, ys]);
clear xs ys u_sym f_sym

%% Integration
for n = 1:length(Nx_set)
    Nx = Nx_set(n);
    Ny = Nx;

    h = L / Nx;
    % hy = L / Ny;

    x = linspace(0, L, Nx);
    y = linspace(0, L, Ny);
    [x, y] = meshgrid(x, y);

    params.x = x;
    params.y = y;
    params.h = h;
    params.Nx = Nx;
    params.Ny = Ny;

    % Integration
    u = system_poisson_dirichlet(scheme, params, f);

    u_true_val = u_true(x, y);
    C_norm(n) = max(abs(u - u_true_val), [], "all");

    if show_vis
        % Visualisation
        figure(1)
        % contour(x, y, u, 'ShowText', 'on');
        surf(x, y, u);
        xl = xlabel('$x$'); xl.Interpreter = 'latex'; xl.FontSize = 16;
        yl = ylabel('$y$'); yl.Interpreter = 'latex'; yl.FontSize = 16;

        % Visualisation 2
        figure(2)
        contour(x, y, log10(abs(u - u_true_val)), 'ShowText', 'on');
        % contour(x, y, log10(abs(u - u_true_val)), 18, 'ShowText', 'on');
        xl = xlabel('$x$'); xl.Interpreter = 'latex'; xl.FontSize = 16;
        yl = ylabel('$y$'); yl.Interpreter = 'latex'; yl.FontSize = 16;
        ttl = title(['$\log_{10} |u - u^*|$']); ttl.Interpreter = 'latex'; ttl.FontSize = 16;
    else
        continue;
    end
end

%% Error study
figure(3)
loglog(Nx_set, C_norm, '-o');
hold on;

C_fit = fitlm(log10(Nx_set), log10(C_norm));
disp("[Order] " + scheme + " scheme: " + num2str(-C_fit.Coefficients.Estimate(2)));
