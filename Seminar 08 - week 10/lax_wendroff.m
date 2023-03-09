function u_new = lax_wendroff(u, f, tau, h)
    pred = (u + circshift(u, -1))/2 ...
        - tau/h * (f(circshift(u, -1)) - f(u))/2;
    u_new = u - tau/h * (f(pred) - f(circshift(pred, 1)));
end
