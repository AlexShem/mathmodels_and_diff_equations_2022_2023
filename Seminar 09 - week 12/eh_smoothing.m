function u_s = eh_smoothing(u, param_jump, param_smooth)
% u - function values to smooth
% param_jump - jump parameter smoothing. gamma <?> 10^(+-)param_jump
% param_smooth - internal point smoothing parameter,
%               S[u] ?> 10^param_smooth * u_max
n = length(u);

% Jump detection
gamma = abs(u - circshift(u, 1)) ./ abs(u - circshift(u, -1));

% On the right from the jump
jind_right = find(gamma > 10^param_jump);
if any(jind_right)
    jind_right1 = jind_right + 1;
    jind_right2 = jind_right + 2;
    % Periodic correction
    jind_right1(jind_right1 == n + 1) = 1;
    jind_right2(jind_right2 == n + 1) = 1;
    jind_right2(jind_right2 == n + 2) = 2;

    u(jind_right) = 2*u(jind_right1) - u(jind_right2);
end

jind_left = find(gamma < 10^-param_jump);
if any(jind_left)
    jind_left1 = jind_left - 1;
    jind_left2 = jind_left - 2;
    % Periodic correction
    jind_left1(jind_left1 == 0) = n;
    jind_left2(jind_left2 == 0) = n;
    jind_left2(jind_left2 == -1) = n - 1;

    u(jind_left) = 2*u(jind_left1) - u(jind_left2);
end

% Non-jump point smoothing
Su = P_smoothing(P_smoothing(P_smoothing(u)))/8;
s_ind = find(abs(Su) > 10^param_smooth * max(u));
s_ind = setdiff(s_ind, [jind_right, jind_left]);

if any(s_ind)
    s_ind_left = s_ind - 1;
    s_ind_right = s_ind + 1;

    % Periodic correction
    n = length(u);
    s_ind_left(s_ind_left == 0) = n;
    s_ind_right(s_ind_right == n + 1) = 1;

    u_s = u;
    u_s(s_ind) = .25*(u(s_ind_left) + 2*u(s_ind) + u(s_ind_right));
else
    u_s = u;
end
end
