clc;
clear;
load('data2.mat');
data1 = data;
load('initial_state2.mat');
initial_state1 = initial_state;
load('IDM2.mat');
IDM1 = IDM;
load('data3.mat');
data2 = data;
load('initial_state3.mat');
initial_state2 = initial_state;
load('IDM3.mat');
IDM2 = IDM;
load('data4.mat');
data3 = data;
load('initial_state4.mat');
initial_state3 = initial_state;
load('IDM4.mat');
IDM3 = IDM;
load('data5.mat');
data4 = data;
load('initial_state5.mat');
initial_state4 = initial_state;
load('IDM5.mat');
IDM4 = IDM;
load('data6.mat');
data5 = data;
load('initial_state6.mat');
initial_state5 = initial_state;
load('IDM6.mat');
IDM5 = IDM;
load('data7.mat');
data6 = data;
load('initial_state7.mat');
initial_state6 = initial_state;
load('IDM7.mat');
IDM6 = IDM;
load('data8.mat');
data7 = data;
load('initial_state8.mat');
initial_state7 = initial_state;
load('IDM8.mat');
IDM7 = IDM;
data = [data1(1:5);data2(1:8);data3(1:143);...
        data4(1:12);data5(1:37);data6(1:669);data7(409)];
initial_state = [initial_state1(1:5);
    initial_state2(1:8);initial_state3(1:143);...
    initial_state4(1:12);initial_state5(1:37);...
    initial_state6(1:669);initial_state7(1:409)];
IDM = [IDM1(1:5);IDM2(1:8);IDM3(1:143);IDM4(1:12);...
       IDM5(1:37);IDM6(1:669);IDM7(1:409)];
% good_data = cell(3000, 1);
% good_initial_state = cell(3000, 1);
% bad_data = cell(3000, 1);
% bad_initial_state = cell(3000, 1);
% fail_data = cell(3000, 1);
% fail_initial_state = cell(3000, 1);
good_initial_state = [];
bad_initial_state = [];
fail_initial_state = [];


%% data load and processing
% initial_state = bad_initial_state;
% data = bad_data;
N = 50;
tf = 5;
lane_width = 3.75;
length = 4.4;
width = 1.7;
is_exist_1 = true; is_exist_2 = true; is_exist_3 = true;
good_ = 0;
bad_ = 0;
fail_ = 0;

%% loop and plot
for i = 1:size(data, 1)
%     figure;
%     axis equal;
%     figure_iter = 1;
%     line([-100,100],[-lane_width/2, -lane_width/2]);
%     hold on;
%     line([-100,100],[0, 0], 'LineStyle', '-.');
%     hold on;
%     line([-100,100],[lane_width/2, lane_width/2]);
%     hold on;
%     line([-100,100],[lane_width, lane_width],'LineStyle', '-.');
%     hold on;
%     line([-100,100],[3*lane_width/2, 3*lane_width/2]);
%     hold on;
    x_pos = data{i, 1}(1:4:(N+1)*4);
    y_pos = data{i, 1}(2:4:(N+1)*4);
    i_ = 1;flag_i_ = true;
    if (y_pos(end) ~= 0)
        while(y_pos(i_) < 3.75/2)
            if (y_pos(i_ + 1) < y_pos(i_))
                flad_i_ = false;
                fprintf('bad solution!\n');
                bad_ = bad_ + 1;
%                 bad_initial_state{bad_, 1} = initial_state{i, 1};
                bad_initial_state(bad_, :) = initial_state{i, 1};
                bad_data{bad_, 1} = data{i, 1};
                bad_IDM{bad_, 1} = IDM{i, 1};
                break;
            end
            i_ = i_ + 1;
        end
        if (y_pos(i_) >= 3.75/2)
            fprintf('good solution!\n')
            good_ = good_ + 1;
    %         good_initial_state{good_, 1} = initial_state{i, 1};
            good_initial_state(good_, :) = initial_state{i, 1}; 
            good_data{good_, 1} = data{i, 1};
            good_IDM{good_, 1} = IDM{i, 1};
        end
        %% plot
%         v_pos = data{i, 1}(3:4:(N+1)*4);
%         theta_pos= data{i, 1}(4:4:(N+1)*4);
%         obs = initial_state{i, 1}(2: end);
%         x_obstacle = repmat(obs(1), N+1, 1) + obs(2)/3.6 * linspace(0, tf, N+1)' +...
%                      0.5 * obs(3) * linspace(0, tf, N+1)'.^2;
%         y_obstacle = zeros(N + 1, 1);
%         x1_obstacle = repmat(obs(4), N+1, 1) + obs(5)/3.6 * linspace(0, tf, N+1)' +...
%                      0.5 * obs(6) * linspace(0, tf, N+1)'.^2;
%         y1_obstacle = ones(N + 1, 1) * lane_width;
%         x2_obstacle = repmat(obs(7), N+1, 1) + obs(8)/3.6 * linspace(0, tf, N+1)' +...
%                      0.5 * obs(9) * linspace(0, tf, N+1)'.^2;
%         y2_obstacle = ones(N + 1, 1) * lane_width;
% 
%         for t_num = 1:N+1
%             obstacle_x_temp = x_pos(t_num, figure_iter);
%             obstacle_y_temp = y_pos(t_num, figure_iter);
%             obstacle_theta_temp =theta_pos(t_num, figure_iter);
%             x = [obstacle_x_temp + length/2 * cos(obstacle_theta_temp) + width/2 * sin(obstacle_theta_temp),...
%                 obstacle_x_temp + length/2 * cos(obstacle_theta_temp) - width/2 * sin(obstacle_theta_temp),...
%                 obstacle_x_temp - length/2 * cos(obstacle_theta_temp) - width/2 * sin(obstacle_theta_temp),...
%                 obstacle_x_temp - length/2 * cos(obstacle_theta_temp) + width/2 * sin(obstacle_theta_temp),...
%                 obstacle_x_temp + length/2 * cos(obstacle_theta_temp) + width/2 * sin(obstacle_theta_temp)];
%             y = [obstacle_y_temp + length/2 * sin(obstacle_theta_temp) - width/2 * cos(obstacle_theta_temp),...
%                 obstacle_y_temp + length/2 * sin(obstacle_theta_temp) + width/2 * cos(obstacle_theta_temp),...
%                 obstacle_y_temp - length/2 * sin(obstacle_theta_temp) + width/2 * cos(obstacle_theta_temp),...
%                 obstacle_y_temp - length/2 * sin(obstacle_theta_temp) - width/2 * cos(obstacle_theta_temp),...
%                 obstacle_y_temp + length/2 * sin(obstacle_theta_temp) - width/2 * cos(obstacle_theta_temp)];
%             l1 = line(x,y, 'color', 'r');
%             hold on;
%             if is_exist_1 == true
%                 obstacle_x_temp = x_obstacle(t_num) - 4.4/2;
%                 obstacle_y_temp = y_obstacle(t_num) - 1.7/2;
%                 x = [obstacle_x_temp, obstacle_x_temp+length, obstacle_x_temp+length, obstacle_x_temp, obstacle_x_temp];
%                 y = [obstacle_y_temp, obstacle_y_temp, obstacle_y_temp+width, obstacle_y_temp+width, obstacle_y_temp];
%                 hold on;
%                 l2 = line(x,y, 'color', 'm');
%             end
%             if is_exist_2 == true
%                 obstacle_x_temp = x1_obstacle(t_num) - 4.4/2;
%                 obstacle_y_temp = y1_obstacle(t_num) - 1.7/2;
%                 x = [obstacle_x_temp, obstacle_x_temp+length, obstacle_x_temp+length, obstacle_x_temp, obstacle_x_temp];
%                 y = [obstacle_y_temp, obstacle_y_temp, obstacle_y_temp+width, obstacle_y_temp+width, obstacle_y_temp];
%                 hold on;
%                 l3 = line(x,y, 'color', 'b');
%             end
%             if is_exist_3 == true
%                 obstacle_x_temp = x2_obstacle(t_num) - 4.4/2;
%                 obstacle_y_temp = y2_obstacle(t_num) - 1.7/2;
%                 x = [obstacle_x_temp, obstacle_x_temp+length, obstacle_x_temp+length, obstacle_x_temp, obstacle_x_temp];
%                 y = [obstacle_y_temp, obstacle_y_temp, obstacle_y_temp+width, obstacle_y_temp+width, obstacle_y_temp];
%                 hold on;
%                 l4 = line(x,y, 'color', 'c');
%             end
%             pause(.05);
%             if (t_num < N + 1)
%                 delete(l1);
%                 if is_exist_1 == true
%                     delete(l2);
%                 end
%                 if is_exist_2 == true
%                     delete(l3);
%                 end
%                 if is_exist_3 == true
%                     delete(l4);
%                 end
%             end
%         end
    else
        fprintf('no solution!\n');
        fail_ = fail_ + 1;
%         fail_initial_state{fail_, 1} = initial_state{i, 1};
        if size(initial_state{i, 1}, 1) == 0
            fprintf('error here!');
        end
        fail_initial_state(fail_, :) = initial_state{i, 1};
        fail_data{fail_, 1} = data{i, 1};
        fail_IDM{fail_, 1} = IDM{i, 1};
    end
%     v_pos = data{i, 1}(3:4:(N+1)*4);
%     theta_pos= data{i, 1}(4:4:(N+1)*4);
%     obs = initial_state{i, 1}(2: end);
%     x_obstacle = repmat(obs(1), N+1, 1) + obs(2)/3.6 * linspace(0, tf, N+1)' +...
%                  0.5 * obs(3) * linspace(0, tf, N+1)'.^2;
%     y_obstacle = zeros(N + 1, 1);
%     x1_obstacle = repmat(obs(4), N+1, 1) + obs(5)/3.6 * linspace(0, tf, N+1)' +...
%                  0.5 * obs(6) * linspace(0, tf, N+1)'.^2;
%     y1_obstacle = ones(N + 1, 1) * lane_width;
%     x2_obstacle = repmat(obs(7), N+1, 1) + obs(8)/3.6 * linspace(0, tf, N+1)' +...
%                  0.5 * obs(9) * linspace(0, tf, N+1)'.^2;
%     y2_obstacle = ones(N + 1, 1) * lane_width;
% 
%     for t_num = 1:N+1
%         obstacle_x_temp = x_pos(t_num, figure_iter);
%         obstacle_y_temp = y_pos(t_num, figure_iter);
%         obstacle_theta_temp =theta_pos(t_num, figure_iter);
%         x = [obstacle_x_temp + length/2 * cos(obstacle_theta_temp) + width/2 * sin(obstacle_theta_temp),...
%             obstacle_x_temp + length/2 * cos(obstacle_theta_temp) - width/2 * sin(obstacle_theta_temp),...
%             obstacle_x_temp - length/2 * cos(obstacle_theta_temp) - width/2 * sin(obstacle_theta_temp),...
%             obstacle_x_temp - length/2 * cos(obstacle_theta_temp) + width/2 * sin(obstacle_theta_temp),...
%             obstacle_x_temp + length/2 * cos(obstacle_theta_temp) + width/2 * sin(obstacle_theta_temp)];
%         y = [obstacle_y_temp + length/2 * sin(obstacle_theta_temp) - width/2 * cos(obstacle_theta_temp),...
%             obstacle_y_temp + length/2 * sin(obstacle_theta_temp) + width/2 * cos(obstacle_theta_temp),...
%             obstacle_y_temp - length/2 * sin(obstacle_theta_temp) + width/2 * cos(obstacle_theta_temp),...
%             obstacle_y_temp - length/2 * sin(obstacle_theta_temp) - width/2 * cos(obstacle_theta_temp),...
%             obstacle_y_temp + length/2 * sin(obstacle_theta_temp) - width/2 * cos(obstacle_theta_temp)];
%         l1 = line(x,y, 'color', 'r');
%         hold on;
%         if is_exist_1 == true
%             obstacle_x_temp = x_obstacle(t_num) - 4.4/2;
%             obstacle_y_temp = y_obstacle(t_num) - 1.7/2;
%             x = [obstacle_x_temp, obstacle_x_temp+length, obstacle_x_temp+length, obstacle_x_temp, obstacle_x_temp];
%             y = [obstacle_y_temp, obstacle_y_temp, obstacle_y_temp+width, obstacle_y_temp+width, obstacle_y_temp];
%             hold on;
%             l2 = line(x,y, 'color', 'm');
%         end
%         if is_exist_2 == true
%             obstacle_x_temp = x1_obstacle(t_num) - 4.4/2;
%             obstacle_y_temp = y1_obstacle(t_num) - 1.7/2;
%             x = [obstacle_x_temp, obstacle_x_temp+length, obstacle_x_temp+length, obstacle_x_temp, obstacle_x_temp];
%             y = [obstacle_y_temp, obstacle_y_temp, obstacle_y_temp+width, obstacle_y_temp+width, obstacle_y_temp];
%             hold on;
%             l3 = line(x,y, 'color', 'b');
%         end
%         if is_exist_3 == true
%             obstacle_x_temp = x2_obstacle(t_num) - 4.4/2;
%             obstacle_y_temp = y2_obstacle(t_num) - 1.7/2;
%             x = [obstacle_x_temp, obstacle_x_temp+length, obstacle_x_temp+length, obstacle_x_temp, obstacle_x_temp];
%             y = [obstacle_y_temp, obstacle_y_temp, obstacle_y_temp+width, obstacle_y_temp+width, obstacle_y_temp];
%             hold on;
%             l4 = line(x,y, 'color', 'c');
%         end
%         pause(.05);
%         if (t_num < N + 1)
%             delete(l1);
%             if is_exist_1 == true
%                 delete(l2);
%             end
%             if is_exist_2 == true
%                 delete(l3);
%             end
%             if is_exist_3 == true
%                 delete(l4);
%             end
%         end
%     end
end

fprintf('%d \n', size(good_initial_state, 1));
fprintf('%d \n', size(bad_initial_state, 1));
fprintf('%d \n', size(fail_initial_state, 1));
good_initial_state_all = good_initial_state;
bad_initial_state_all = bad_initial_state;
fail_initial_state_all = fail_initial_state;
number_data_svm = floor(size(fail_initial_state, 1) * 0.8);
good_initial_state = good_initial_state(1:number_data_svm, :);
bad_initial_state = bad_initial_state(1:number_data_svm, :);
fail_initial_state = fail_initial_state(1:number_data_svm, :);
% record_ = zeros(111, 1);
% for j =1:111
%      record_(j) = trainedModel.predictFcn(bad_initial_state_all(j,:))
% end
state_ = [good_initial_state; bad_initial_state; fail_initial_state];
% out_ = [-ones(size(good_initial_state, 1), 1); zeros(size(bad_initial_state, 1), 1);...
%     ones(size(fail_initial_state, 1), 1)];
% out_ = [repmat([1, 0, 0], size(good_initial_state, 1), 1);...
%         repmat([0, 1, 0], size(bad_initial_state, 1), 1);...
%         repmat([0, 0, 1], size(fail_initial_state, 1), 1)];
response_ = [repmat([-1], size(good_initial_state, 1), 1);...
        repmat([0], size(bad_initial_state, 1), 1);...
        repmat([1], size(fail_initial_state, 1), 1)];
index = randperm(size(state_, 1));
state_rand_all = state_(index, :);
% out_rand_all = out_(index, :);
response_rand_all = response_(index, :);
train_size = round(size(state_, 1));
state_rand_ = state_rand_all(1:train_size, :);
% out_rand_ = out_rand_all(1:train_size, :);
response_rand_ = response_rand_all(1:train_size, :);
% [input_,ps3]=mapminmax(state_rand_');
% [target_,ps4]=mapminmax(out_rand_');
svm_data = [state_rand_ response_rand_];
% yfit = trainedModel.predictFcn([25.8940467249387,21.5783722707823,27.0219686482920,0.891405292200046,24.0828490186278,39.8585688260799,0.675047029620226,-96.8803356517333,42.9391073634308,-0.825750851773087]);

%% test the model
load('trainedModel.mat');
test_set = [good_initial_state_all(number_data_svm:end, :); bad_initial_state_all(number_data_svm:end, :);...
            fail_initial_state_all(number_data_svm:end, :)];
test_out = zeros(size(test_set, 1), 1);
for j = 1:size(test_set, 1)
    test_out(j) = trainedModel.predictFcn(test_set(j, :));
end
num_good = size(good_initial_state_all(number_data_svm:end, :), 1);
num_bad = size(bad_initial_state_all(number_data_svm:end, :), 1);
num_fail = size(fail_initial_state_all(number_data_svm:end, :), 1);
num_good_test = size(find(test_out(1:num_good) == -1), 1);
num_bad_test = size(find(test_out(num_good+1:num_good+num_bad) == 0), 1);
num_fail_test = size(find(test_out(num_good+num_bad+1:num_good+num_bad+...
    num_fail) == 1), 1);
fprintf('%d \n', num_good_test/num_good);
fprintf('%d \n', num_bad_test/num_bad);
fprintf('%d \n', num_fail_test/num_fail);
% result:
% 9.506048e-01 
% 9.782609e-01 
% 1 

%% retrain_data
nx = 4; nu = 2; N = 50;
temp = good_data{1, 1};
temp_IDM = good_IDM{1, 1};
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
    a_pos = temp(4*51+1:2:304);
    omega_pos = temp(4*51+2:2:304);
    obstacle_pos = temp_IDM(1:end-1, :);
    % obstacle0 in the whole trajectory
    x_obstacle = obstacle_pos(:, 1);
    y_obstacle = zeros(N + 1, 1);
    v_obstacle = obstacle_pos(:, 2);
    theta_obstacle = ones(51,1) * 0/180*pi;
    a_obstacle = obstacle_pos(:, 3);

    % obstacle1 in the whole trajectory
    x1_obstacle = obstacle_pos(:, 4);
    y1_obstacle = ones(N + 1, 1) * lane_width;
    v1_obstacle = obstacle_pos(:, 5);
    theta1_obstacle = ones(51, 1) * 0/180*pi;
    a1_obstacle = obstacle_pos(:, 6);

    % obstacle2 in the whole trajectory
    x2_obstacle = obstacle_pos(:, 7);
    y2_obstacle = ones(N + 1, 1) * lane_width;
    v2_obstacle = obstacle_pos(:, 8);
    theta2_obstacle = ones(N + 1, 1) * 0/180*pi;
    a2_obstacle = obstacle_pos(:, 9);
    
    % store information for 3 obstacles
    all_obstacle_state = [x_obstacle(1:50),...
                          x1_obstacle(1:50),...
                          x2_obstacle(1:50)];
    % relative states
    all_obstacle_state_r = [x_obstacle(1:50) - x_pos, v_obstacle(1:50) - v_pos,...% a_obstacle(1:50) - a_pos,...
                          x1_obstacle(1:50) - x_pos, v1_obstacle(1:50) - v_pos,...% a1_obstacle(1:50) - a_pos,...
                          x2_obstacle(1:50) - x_pos, v2_obstacle(1:50) - v_pos]; %, a2_obstacle(1:50) - a_pos];
    temp_initial = repmat(good_initial_state_all(1, :), 50, 1);
    
end
for i=2:size(good_data, 1)
    temp=good_data{i,1};
    if abs(temp(4*51-2) - 3.75) <= .1
        x_pos_temp = temp(1:4:(N+1)*nx-4);
        y_pos_temp = temp(2:4:4*51-4);
        v_pos_temp = temp(3:4:4*51-4);
        theta_pos_temp = temp(4:4:4*51-4);
        a_pos_temp = temp(4*51+1:2:304);
        
        x_pos = [x_pos;temp(1:4:(N+1)*nx-4)];
        y_pos = [y_pos;temp(2:4:4*51-4)];
        v_pos = [v_pos;temp(3:4:4*51-4)];
        theta_pos = [theta_pos;temp(4:4:4*51-4)];
        omega_pos = [omega_pos;temp(4*51+2:2:304)];
        a_pos = [a_pos;temp(4*51+1:2:304)];
            % obstacle0 in the whole trajectory
        temp_IDM = good_IDM{i, 1};
        x_obstacle = obstacle_pos(:, 1);
        y_obstacle = zeros(N + 1, 1);
        v_obstacle = obstacle_pos(:, 2);
        theta_obstacle = ones(51,1) * 0/180*pi;
        a_obstacle = obstacle_pos(:, 3);

        % obstacle1 in the whole trajectory
        x1_obstacle = obstacle_pos(:, 4);
        y1_obstacle = ones(N + 1, 1) * lane_width;
        v1_obstacle = obstacle_pos(:, 5);
        theta1_obstacle = ones(51, 1) * 0/180*pi;
        a1_obstacle = obstacle_pos(:, 6);

        % obstacle2 in the whole trajectory
        x2_obstacle = obstacle_pos(:, 7);
        y2_obstacle = ones(N + 1, 1) * lane_width;
        v2_obstacle = obstacle_pos(:, 8);
        theta2_obstacle = ones(N + 1, 1) * 0/180*pi;
        a2_obstacle = obstacle_pos(:, 9);

        % store information for 3 obstacles
        all_obstacle_state = [x_obstacle(1:50),...
                              x1_obstacle(1:50),...
                              x2_obstacle(1:50)];
        % relative states
        all_obstacle_state_r = [x_obstacle(1:50) - x_pos_temp, v_obstacle(1:50) - v_pos_temp,...% a_obstacle(1:50) - a_pos_temp,...
                              x1_obstacle(1:50) - x_pos_temp, v1_obstacle(1:50) - v_pos_temp,...% a1_obstacle(1:50) - a_pos_temp,...
                              x2_obstacle(1:50) - x_pos_temp, v2_obstacle(1:50) - v_pos_temp]; %, a2_obstacle(1:50) - a_pos_temp];
        temp_initial = [temp_initial; repmat(good_initial_state_all(i, :), 50, 1)];
    end
end
global ps1 ps2;
state = [v_pos x_pos y_pos theta_pos temp_initial]';

%% data processing
input = [];
target = [];
% for i = 1:size(state, 2)/50
%     state_temp = state(:,(i-1)*50+1:i*50);
%     ini_temp = state_temp(:, 1);
%     x1_obstacle = repmat(ini_temp(6), N+1, 1) + ini_temp(7)/3.6 * linspace(0, tf, N+1)' +...
%                  0.5 * ini_temp(8) * linspace(0, tf, N+1)'.^2;    
%     x2_obstacle = repmat(ini_temp(9), N+1, 1) + ini_temp(10)/3.6 * linspace(0, tf, N+1)' +...
%                  0.5 * ini_temp(11) * linspace(0, tf, N+1)'.^2;
%     x3_obstacle = repmat(ini_temp(12), N+1, 1) + ini_temp(13)/3.6 * linspace(0, tf, N+1)' +...
%                  0.5 * ini_temp(14) * linspace(0, tf, N+1)'.^2;
% end
[input,ps1]=mapminmax(state);
[target,ps2]=mapminmax([a_pos omega_pos]');

%% train the net
% net = net;
net = newff(input,target,10,{'tansig','purelin'},'trainlm');
net.trainParam.epochs=1000000;%最大训练次数
net.trainParam.goal=0.00001;%目标最小误差
LP.lr=0.00000001;%学习速率
net.trainParam.max_fail=100;  
% net=train(net,input,target);

%%
load('net.mat');
load('ps1.mat');
load('ps2.mat');
global net;
global temp_initial good_data good_IDM good_initial_state_all
global initial_i IDM_i
path_record_size = size(temp_initial, 1)/50;
path_record = cell(path_record_size, 1);
for i = 1:path_record_size
    initial_i = temp_initial((i - 1) * 50 + 1, :);
    IDM_i = good_IDM{i, 1};
    x_obstacle = IDM_i(1, 1);
    v_obstacle = IDM_i(1, 2);
    x1_obstacle = IDM_i(1, 4);
    v1_obstacle = IDM_i(1, 5);
    x2_obstacle = IDM_i(1, 7);
    v2_obstacle = IDM_i(1, 8);
    [t, record_temp] = ode45(@ode, [0:0.1:5],[0 0 initial_i(1)/3.6 0 x_obstacle v_obstacle ...
        x1_obstacle v1_obstacle x2_obstacle v2_obstacle]);
    path_record{i} = [t record_temp]; 
end


%% 
path_index = 3;
record__1 = path_record{path_index};
plot(record__1(:, 2), record__1(:, 3), '-');
hold on;
plot(record__1(:,6), repmat(0, 51, 1));
hold on;
plot(record__1(:,8), repmat(3.75, 51, 1));
hold on;
plot(record__1(:,10), repmat(3.75, 51, 1));
i = path_index;
x_real = x_pos((i - 1) * 50 + 1:i * 50);
y_real = y_pos((i - 1) * 50 + 1:i * 50);
hold on;
plot(x_real, y_real,'*');
theta_real = theta_pos((i - 1) * 50 + 1:i * 50);
plot(theta_real,'-.');
figure;
plot(record__1(:, 4), '-');
hold on;
v_real = v_pos((i - 1) * 50 + 1:i * 50);
plot(v_real,'-.');

