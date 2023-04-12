function [U, Nt, T, x, h, t] = B_integration(L, T, D, Nx, tau, nu, u_0, correction)
%% Secondary paramters
x = linspace(0, L, Nx + 1);
h = x(2) - x(1);

Nt = ceil(T/tau) + 1;
T = tau*(Nt - 1);
t = 0 : tau : T;

%% Integration
U = zeros(Nt, Nx);
U(1, :) = u_0(x(1 : end-1));

for k = 2 : Nt
    u = U(k-1, :);
    u_ex = burgers_adams_bashforth(u, D, tau, h);
    if strcmp(correction, 'compact')
        comp_eps = burgers_compact_correction(u, u_ex, D, h, tau);
    else
        comp_eps = zeros(length(u), 1);
    end
    u_k = u_ex + comp_eps.';
    U(k, :) = u_k;
end
U = [U, U(:, 1)];
end
