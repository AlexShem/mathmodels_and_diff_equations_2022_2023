function vareps = burgers_compact_correction(u, u_ex, D, h, tau)
nu = D*tau/h^2;

a1 = -2/3 + 2*nu; b1 = -8/3 - 4*nu; c1 = -2/3 + 2*nu;
a0 = 2/3 + 2*nu; b0 = 8/3 - 4*nu; c0 = 2/3 + 2*nu;

p1 = tau/h; q1 = 0; r1 = -tau/h;
p0 = tau/h; q0 = 0; r0 = -tau/h;

F = -a1*circshift(u_ex, 1) - b1*u_ex - c1*circshift(u_ex, -1) ...
    -a0*circshift(u, 1) - b0*u - c0*circshift(u, -1) ...
    -p0*circshift(u, 1).^2/2 - q0*u.^2/2 - r0*circshift(u, -1).^2/2 ...
    -p1*circshift(u_ex, 1).^2/2 - q1*u_ex.^2/2 - r1*circshift(u_ex, -1).^2/2;
A = diag(b1 + q1*u_ex) + ...
    diag(c1 + r1*u_ex(2:end), 1) + ...
    diag(a1 + p1*u_ex(1:end-1), -1);
A(1, end) = a1 + p1*u_ex(end);
A(end, 1) = c1 + r1*u_ex(1);

vareps = sparse(A) \ F.';
end
