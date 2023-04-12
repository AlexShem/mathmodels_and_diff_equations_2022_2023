function u_new = explicit_euler(u, f, nu)
    u_new = u + nu*(f(circshift(u, -1)) - 2*f(u) + f(circshift(u, 1)));
end
