N_set = 2.^(2:10);
N_iter = zeros(size(N_set));
tol = 1e-6;
method = 0;

for k = 1 : length(N_set)
    N = N_set(k);
    a = N;
    A = unifrnd(-1, 1, N, N);

    % A(abs(A) < .5) = 0;
    n1 = norm(A, 1);
    nInf = norm(A, "inf");
    colsum = sum(abs(A), 2);
    rowsum = sum(abs(A), 1);
    maxsum = max(colsum, rowsum.');

    A(logical(eye(size(A)))) = ceil(maxsum);
    A(abs(A) < .8) = 0;

    b = normrnd(0, 1, N, 1);
    x0 = zeros(N, 1);

    [x, Ni] = seidel_method(A, b, x0, tol, method);
    N_iter(k) = Ni;
end

%% Visualisation
figure(1);
if method == 1
    semilogx(N_set, N_iter, '--');
else
    semilogx(N_set, N_iter);
end
hold on;
