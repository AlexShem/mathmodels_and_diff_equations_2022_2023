function u_new = burgers_adams_bashforth(u, D, tau, h)
    nu = D*tau/h^2;
    u_jm = circshift(u, 1);
    u_jp = circshift(u, -1);
    u_p = u + nu*(u_jm - 2*u + u_jp)/2 + tau*(u_jm.^2 - u_jp.^2)/(4*h);
    u_pm = circshift(u_p, 1);
    u_pp = circshift(u_p, -1);
    u_new = u + nu*(u_pm - 2*u_p + u_pp) + tau*(u_pm.^2 - u_pp.^2)/(4*h);
end
