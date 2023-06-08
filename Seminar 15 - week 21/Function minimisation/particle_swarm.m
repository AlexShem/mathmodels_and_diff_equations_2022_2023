function [xmin, fmin] = particle_swarm(f, ff, lb, ub, S, niter)
% lb - vector of lower bounds across all dimensions
% ub - vector of upper bounds across all dimensions
if nargin < 4
    S = 50; % Number of particles
end
if nargin < 5
    niter = 100;
end
% omega = .1;
% phip = 1.8;
% phig = 2.5;
omega = 0.5;
phip = 1.2;
phig = .9;

d = length(lb); % Number of dimensions

X = rand(S, d);
% Scale and shift to the propper domain
for i = 1 : d
    X(:, i) = X(:, i) * (ub(i) - lb(i)) + lb(i);
end

P = X; % Best particle's location

% Best location found by the swarm's
Fp = arrayfun(@(ind) f(P(ind, :)), (1:S).');
[Fg, g_ind] = min(Fp);
G = P(g_ind, :);

% Velosity of particles
V = rand(S, d);
for i = 1 : d
    span = ub(i) - lb(i);
    V(:, i) = 2*span*(V(:, i) -.5);
end

% fig = figure("Name", 'Particle Swarm animation');
fig = figure(3);
fsurf(ff, [lb(1) ub(1) lb(2) ub(2)], 'ShowContours', 'on', 'FaceAlpha', .3)
% fcontour(ff, [lb(1) ub(1) lb(2) ub(2)], '-k', 'Fill', 'on');
hold on;
% scatter(X(:, 1), X(:, 2), '*r');
scatter3(X(:, 1), X(:, 2), Fp, 30, 'or', 'filled');
hold off;
xlabel('x');
ylabel('y');

iter = 0;
while iter <= niter
    rp = rand(S, 1);
    rg = rand(S, 1);

    V = omega*V + ...
        phip*repmat(rp, 1, 2).*(P - X) + ...
        phig*repmat(rg, 1, 2).*(G - X);
    X = X + V;
    % Shift X back to the domain
    for i = 1 : d
        X(:, i) = max(X(:, i), lb(i));
        X(:, i) = min(X(:, i), ub(i));
    end

    Fx = arrayfun(@(ind) f(X(ind, :)), (1:S).');
    F_imp = Fx < Fp;
    P(F_imp, :) = X(F_imp, :);
    Fp(F_imp) = arrayfun(@(ind) f(P(ind, :)), find(F_imp));

    F_imp_global = Fp < Fg;
    if any(F_imp_global)
        [Fg, g_ind] = min(Fp);
        G = P(g_ind, :);
    end
    iter = iter + 1;

%     figure(fig);
    figure(3);
    fsurf(ff, [lb(1) ub(1) lb(2) ub(2)], 'ShowContours', 'on', 'FaceAlpha', .3)
%     fcontour(ff, [lb(1) ub(1) lb(2) ub(2)], '-k', 'Fill', 'on');
    hold on;
%     scatter(X(:, 1), X(:, 2), '*r');
    scatter3(X(:, 1), X(:, 2), Fx, 30, 'or', 'filled');
    hold off;
    xlabel('x');
    ylabel('y');
    drawnow;
end

xmin = G;
fmin = Fg;
end
