function Pu = P_smoothing(u)
Pu = u - (circshift(u, -1) + circshift(u, 1))/2;
end
