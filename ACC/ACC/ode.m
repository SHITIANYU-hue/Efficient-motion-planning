function [outputArg1] = ode(t, state)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
x = state(1);
y = state(2);
v = state(3);
theta = state(4);
global ps1 ps2 net initial_i acc;
initial_state = initial_i;
input=mapminmax('apply',[v; x; y; theta; initial_state'], ps1);
output=net(input);
[control]=mapminmax('reverse',output,ps2);
dx = v * cos(theta);
dy = v * sin(theta);
dv = control(1);
acc = [acc, control(1)];
dtheta = control(2);
outputArg1 = [dx, dy, dv, dtheta]';
end

