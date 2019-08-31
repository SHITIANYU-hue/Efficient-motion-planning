function [outputArg1] = ode(~, state)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
width = 1.7;
length = 4.4;
a_max = 6 ;
lane_width = 3.75;
x = state(1);
y = state(2);
v = state(3);
theta = state(4);
x_obstacle = state(5);
v_obstacle = state(6);
x1_obstacle = state(7);
v1_obstacle = state(8);
x2_obstacle = state(9);
v2_obstacle = state(10);
global ps1 ps2 net initial_i
temp = [v; x; y; theta; initial_i'];
input=mapminmax('apply',temp, ps1);
output=net(input);
[control]=mapminmax('reverse',output,ps2);
dx = v * cos(theta);
dy = v * sin(theta);
dv = control(1);
dtheta = control(2);
% IDM for x_obstacle
% the following car is following the leading car
v0 = 30/3.6;
a = a_max;  % maximum acceleration
a_obstacle = a * (1 - (v_obstacle / v0)^4);
d_x_obstacle = v_obstacle;
d_v_obstacle = a_obstacle;
% IDM for x1_obstace
% the following car is following the leading car
v0 = 30/3.6;
a = a_max;  % maximum acceleration
a1_obstacle = a * (1 - (v1_obstacle / v0)^4);
d_x1_obstacle = v1_obstacle;
d_v1_obstacle = a1_obstacle;
% IDM for x2_obstacle
if x + width/2 <= lane_width/2 % ego car is in the ego lane
    % the following car is following the leading car
    v0 = 30/3.6;
    v_head = v1_obstacle;
    x_head = x1_obstacle;
    sn = x_head - x2_obstacle - length;
    s0 = 2; % minimum distance
    T = 1.5;  % safe time headway
    a = a_max;  % maximum acceleration
    b = 1.67;  % comfortable deceleration
    s_ = s0 + v2_obstacle*T + v2_obstacle * (v2_obstacle - v_head) / (2 * sqrt(a * b));
    a2_obstacle = a * (1 - (v2_obstacle / v0)^4 - (s_/sn)^2);
    d_x2_obstacle = v2_obstacle;
    d_v2_obstacle = a2_obstacle;

else
    v0 = 30/3.6;
    v_head = v;
    x_head = x;
    sn = x_head - x2_obstacle - length;
    s0 = 2; % minimum distance
    T = 1.5;  % safe time headway
    a = a_max;  % maximum acceleration
    b = 1.67;  % comfortable deceleration
    s_ = s0 + v2_obstacle*T + v2_obstacle * (v2_obstacle - v_head) / (2 * sqrt(a * b));
    a2_obstacle = a * (1 - (v2_obstacle / v0)^4 - (s_/sn)^2);
    d_x2_obstacle = v2_obstacle;
    d_v2_obstacle = a2_obstacle;

end
outputArg1 = [dx, dy, dv, dtheta, d_x_obstacle, d_v_obstacle,...
    d_x1_obstacle, d_v1_obstacle, d_x2_obstacle, d_v2_obstacle]';
end

