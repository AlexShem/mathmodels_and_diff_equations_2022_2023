function u_new = interp_compact(u)
u = u(:);
a = 2/3;
p = 1/6;
Nx = length(u);
Nf = Nx - 1;

A = a*(diag(ones(Nf, 1)) + diag(ones(Nf - 1, 1), 1));

if mod(Nf, 2) == 0
    A(end, :) = (-1).^(1:Nf);
else
    A(end, 1:2:end) = 2/Nx;
    A(end, 2:2:end) = -2/Nx;
end

F = eye(Nx) + p*(diag(ones(Nx - 1, 1), 1) + diag(ones(Nx - 1, 1), -1));
F([1, end], :) = [];
F(end + 1, :) = 0;

u_int = A \ (F * u);
u_new = [u.'; [u_int.', NaN]];
u_new = u_new(1 : (Nx + Nf));
end
