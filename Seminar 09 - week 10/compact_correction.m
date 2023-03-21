function vareps = compact_correction(u, u_ex, h, tau)
a1 = 1; b1 = 4; c1 = 1;
a0 = -1; b0 = -4; c0 = -1;
p1 = -3*tau/(2*h); q1 = 0; r1 = 3*tau/(2*h);
p0 = -3*tau/(2*h); q0 = 0; r0 = 3*tau/(2*h);

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
