function [xmin, fmin, niter, path] = grad_descent(x0, f, tol, maxiter)
if nargin < 4
    maxiter = 50;
end
if nargin < 3 || (nargin == 4 && isempty(tol))
    tol = 1e-8;
end

x0 = x0(:);
path = zeros(length(x0), maxiter + 1);
path(:, 1) = x0;

lam = 1;
grad_x = grad_f(x0, f);
x = x0 - grad_x;
iter = 1;
path(:, iter + 1) = x;

df = abs(f(x) - f(x0));

while (df > tol) && iter < maxiter
    grad_x = grad_f(x, f);
    
    x0 = x;
    x = x0 - lam * grad_x/norm(grad_x, 2);
    while f(x) > f(x0)
        lam = lam/2;
        x = x0 - lam * grad_x/norm(grad_x, 2);
    end
    
    iter = iter + 1;
    df = abs(f(x) - f(x0));
    path(:, iter + 1) = x;
end

if iter == maxiter && df > tol
    warning('Gradient descent algorithm did not converge properly');
end
xmin = x;
fmin = f(x);
niter = iter;
path = path(:, 1 : niter + 1);

function grad = grad_f(x, f)
d = length(x);
h = 1e-8;
grad = arrayfun(@(ind) (f(x + h*((1:d).' == ind)) - f(x))/h, (1:d).');
