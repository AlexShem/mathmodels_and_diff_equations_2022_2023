function u_new = explicit_euler(u, f, tau, h)
%     u_new = u - tau*(circshift(u, 1).^2 - circshift(u, -1).^2)/(4*h);
    u_new = u - tau*(f(circshift(u, -1)) - f(circshift(u, 1)))/(2*h);
end
