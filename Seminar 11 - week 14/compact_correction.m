function vareps = compact_correction(u, u_ex, c, h, tau)
a1 = 1; b1 = 10;
a0 = -1; b0 = -10;

p1 = 6*tau/h^2; q1 = -12*tau/h^2;
p0 = 6*tau/h^2; q0 = -12*tau/h^2;

F = -a1*circshift(u_ex, 1) - b1*u_ex - a1*circshift(u_ex, -1) ...
    -a0*circshift(u, 1) - b0*u - a0*circshift(u, -1) ...
    +c*(p0*circshift(u, 1).^2/2 + q0*u.^2/2 + p0*circshift(u, -1).^2/2 ...
    +p1*circshift(u_ex, 1).^2/2 + q1*u_ex.^2/2 + p1*circshift(u_ex, -1).^2/2);
A = diag(b1 - q1*u_ex) + ...
    diag(a1 - p1*u_ex(2:end), 1) + ...
    diag(a1 - p1*u_ex(1:end-1), -1);
A(1, end) = a1 - p1*u_ex(end);
A(end, 1) = a1 - p1*u_ex(1);

vareps = sparse(A) \ F.';
end
