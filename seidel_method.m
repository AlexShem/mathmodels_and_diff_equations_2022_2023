function [x, Niter] = seidel_method(A, b, x0, tol, method)
N = length(b);
x = Inf(N);
err = norm(x - x0, 2);
Niter = 0;

while err > tol
    if method == 1
        ind_ord = randperm(N);
    else
        ind_ord = 1:N;
    end
    x = x0;
    for k  = 1 : N
        i = ind_ord(k);
        r_new = ind_ord(1 : (k-1));
        r_old = ind_ord((k+1) : end);

        s1 = A(i, r_new) * x(r_new);
        s2 = A(i, r_old) * x(r_old);
        x(i) = (b(i) - s1 - s2)/A(i, i);
    end
    err = norm(x - x0, 2);
    x0 = x;
    Niter = Niter + 1;
    if Niter > 500
        error('Seidel method did not converge');
    end
end
end
