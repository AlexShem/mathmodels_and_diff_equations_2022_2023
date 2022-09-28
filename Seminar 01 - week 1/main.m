%% Initial conditions

% Speed
a = 2; b = 1;

% u = @(t, x) a*x + b;
% dudt = @(t, x) zeros(size(x));
% dudx = @(t, x) a*ones(size(x));

u = @(t, x) -a*x/200 + .5*t.*x;
dudt = @(t, x) 2*pi*sin(2*pi*t);
dudx = @(t, x) a*pi*cos(pi*x);

% Integration time
T = 1;
t = linspace(0, T, 11);

% Segment of interest
x = linspace(-1, 1, 11);

% Initial density distribution
rho_fun = @(x) .5*normpdf(x, 0, .5) + ...
    .25*normpdf(x, -.5, .25) + ...
    .25*normpdf(x, .5, .25);
rho = rho_fun(x);


figure(1);
fplot(rho_fun, [min(x), max(x)]);
xlabel('$x$', Interpreter = 'latex', FontSize = 16);
legend('$\rho(0, x)$', Interpreter = 'latex', FontSize = 16);

%% Characteristics

x0 = zeros(size(x));

figure(2)
for k = 1 : length(x)
    [t_path, x_path] = ode45(u, [T 0], x(k));
    plot(x_path, t_path);
    hold on;
    x0(k) = x_path(end);

    % Verification
%     x_path = x0(k) * exp(a*t) + b/a*(exp(a*t) - 1);
%     plot(x_path, t, '*k');
end
hold off;
xlabel('$x$', Interpreter = 'latex', FontSize = 16);
ylabel('$t$', Interpreter = 'latex', FontSize = 16);

%% Interpolation of rho at x0
rho_0 = spline(x, rho, x0);

figure(3)
fplot(rho_fun, [min([x, x0]), max([x, x0])]);
hold on;
plot(x0, rho_0, '*r');
hold off;
