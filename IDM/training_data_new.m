%% settings
N = 50; nx = 4; nu = 2; tf = 5; lane_width = 3.75;

%% load data
load('data1.mat');
data1 = data1;
load('initial_state1.mat');
initial_state1 = initial_state1;
load('data2.mat');
data2 = data;
load('initial_state2.mat');
initial_state2 = initial_state2;
load('data3.mat');
data3 = data;
load('initial_state3.mat');
initial_state3 = initial_state3;
data = [data1(1:720);data2(1:1369);data3(1:358)];
initial_state = [initial_state1(1:720);
    initial_state2(1:1369);initial_state3(1:358)];
temp = data{2,1};
x_pos = [];
y_pos = [];
v_pos = [];
theta_pos = [];
a_pos = [];
omega_pos = [];
temp_initial = [];
if abs(temp(4*51 - 2) - 3.75) < .1
    x_pos = temp(1:4:(N+1)*nx-4);
    y_pos = temp(2:4:4*51-4);
    v_pos = temp(3:4:4*51-4);
    theta_pos = temp(4:4:4*51-4);
    a_pos = temp(4*51+1:2:end);
    omega_pos = temp(4*51+2:2:end);
    temp_initial = repmat(initial_state{2, 1}, 50, 1);
end
for i=3:size(data, 1)
    temp=data{i,1};
    if abs(temp(4*51-2) - 3.75) <= .1
        x_pos = [x_pos;temp(1:4:(N+1)*nx-4)];
        y_pos = [y_pos;temp(2:4:4*51-4)];
        v_pos = [v_pos;temp(3:4:4*51-4)];
        theta_pos = [theta_pos;temp(4:4:4*51-4)];
        omega_pos = [omega_pos;temp(4*51+2:2:end)];
        a_pos = [a_pos;temp(4*51+1:2:end)];
        temp_initial = [temp_initial; repmat(initial_state{i, 1}, 50, 1)];
        if (size(temp_initial, 1) < size(v_pos, 1))
        end
    end
end
global ps1 ps2;
state = [v_pos x_pos y_pos theta_pos temp_initial]';

%% data processing
input = [];
target = [];
for i = 1:size(state, 2)/50
    state_temp = state(:,(i-1)*50+1:i*50);
    ini_temp = state_temp(:, 1);
    x_obstacle = repmat(ini_temp(6), N+1, 1) + ini_temp(7)/3.6 * linspace(0, tf, N+1)' +...
                 0.5 * ini_temp(8) * linspace(0, tf, N+1)'.^2; 
    x_obstacle = x_obstacle(1:50) - state_temp(2, :)';
    x1_obstacle = repmat(ini_temp(9), N+1, 1) + ini_temp(10)/3.6 * linspace(0, tf, N+1)' +...
                 0.5 * ini_temp(11) * linspace(0, tf, N+1)'.^2;
    x1_obstacle = x1_obstacle(1:50) - state_temp(2, :)';
    x2_obstacle = repmat(ini_temp(12), N+1, 1) + ini_temp(13)/3.6 * linspace(0, tf, N+1)' +...
                 0.5 * ini_temp(14) * linspace(0, tf, N+1)'.^2;
    x2_obstacle = x2_obstacle(1:50) - state_temp(2, :)';
    x_obstacle_v = ini_temp(7)/3.6 * ones(50, 1) +...
                 ini_temp(8) * linspace(0, tf, N)'; 
    x1_obstacle_v = ini_temp(10)/3.6 * linspace(0, 1, N)' +...
                 ini_temp(8) * linspace(0, tf, N)'; 
    x2_obstacle_v = ini_temp(13)/3.6 * linspace(0, 1, N)' +...
                 ini_temp(8) * linspace(0, tf, N)'; 
    % IDM for x2_obstacle
    x2_obstacle_temp = x2_obstacle(0);
    for j = 1:1:N
        state_temp_ = state_temp(:, j);
        if state_temp_(2) + 1.7/2 <= lane_width/2 % ego car is in the ego lane
            % the following car is following the leading car
            v0 = 30/3.6;
            v1_obstacle = x1_obstacle_v(j);
            sn = x1_obstacle(j) - x2_obstacle_temp - length;
            s0 = 2; % minimum distance
            T = 1.5;  % safe time headway
            a = a_max;  % maximum acceleration
            b = 1.67;  % comfortable deceleration
            s_ = s0 + v2_obstacle*T + v2_obstacle * (v2_obstacle - x1_obstacle_v(j)) / (2 * sqrt(a * b));
            a2_obstacle = a * (1 - (x2_obstacle_v / v0)^4 - (s_/sn)^2);
            v2_obstacle = a2_obstacle * dt + v2_obstacle;
            x2_obstacle_temp = v2_obstacle * dt + x2_obstacle_temp;
            x2_obstacle = [x2_obstacle; x2_obstacle_temp];
        else
            v0 = 30/3.6;
            v1_obstacle = x((i-1)*nx+3);
            sn = x((i-1)*nx+1) - x2_obstacle - length;
            s0 = 2; % minimum distance
            T = 1.5;  % safe time headway
            a = a_max;  % maximum acceleration
            b = 1.67;  % comfortable deceleration
            s_ = s0 + v2_obstacle*T + v2_obstacle * (v2_obstacle - v1_obstacle) / (2 * sqrt(a * b));
            a2_obstacle = a * (1 - (v2_obstacle / v0)^4 - (s_/sn)^2);
            v2_obstacle = a2_obstacle * dt + v2_obstacle;
            x2_obstacle = v2_obstacle * dt + x2_obstacle;
        end
    end
    temp = [state_temp(1:4, :); x1_obstacle'; x2_obstacle'; x3_obstacle'; ...
            x1_obstacle_v'; x2_obstacle_v'; x3_obstacle_v'; state_temp([8,11,14], :); [1:50]];
    input = [input temp];
end
[input,ps1]=mapminmax(input);
[target,ps2]=mapminmax([a_pos omega_pos]');

%% train the net
% net = net;
net = newff(input,target,6,{'tansig','purelin'},'trainlm');
net.trainParam.epochs=1000000;%最大训练次数
net.trainParam.goal=0.00001;%目标最小误差
LP.lr=0.00000001;%学习速率
net.trainParam.max_fail=100;  
net=train(net,input,target);

%%
global net;
global initial_state data
global initial_i
% dnn_output = net(input);
% prediction1=mapminmax('reverse',dnn_output,ps2);
path_record_size = size(temp_initial, 1)/50;
path_record = cell(path_record_size, 1);
for i = 1:path_record_size
    initial_i = temp_initial((i - 1) * 50 + 1, :);
    [t, record_temp] = ode45(@ode,[0:0.1:5],[0 0 initial_i(1)/3.6 0]);
    path_record{i} = [t record_temp];
    
end
