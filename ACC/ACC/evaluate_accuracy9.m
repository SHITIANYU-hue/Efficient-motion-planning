bad_initial = [];
good_initial = [];
path_record_size = 4022;
for i = 1:path_record_size
    record = path_record{i};
    initial = temp_initial((i-1) * 50 + 1, :);
    x = record(:, 2);
    y = record(:, 3);
    v = record(:, 4);
    theta = record(:, 5);
    x_real = x_pos(i * 50);
    y_real = y_pos(i * 50);
    v_real = v_pos(i * 50);
    theta_real = theta_pos(i * 50);
    x_differ_temp(i) = abs(x(end) - x_real);
    y_differ_temp(i) = abs(y(end) - y_real);
    v_differ_temp(i) = abs(v(end) - v_real);
    theta_differ_temp(i) = abs(theta(end) - theta_real);
%     figure;
%     plot(x, y, '*');
%     hold on;
%     plot(x_pos((i-1)*50+1:i*50), y_pos((i-1)*50+1:i*50), '*');
    if (x_differ_temp(i) > 4 || y_differ_temp(i) > 0.33 || v_differ_temp(i) > 1 || ...
        theta_differ_temp(i) > 0.05)
%         figure;
%         plot(x, y, '*');
        bad_initial = [bad_initial; [v(1) initial]];
    else
%         figure;
%         plot(x, y, '*');
        good_initial = [good_initial; [v(1) initial]];
    end
end
x_differ = mean(x_differ_temp);
y_differ = mean(y_differ_temp);
v_differ = mean(v_differ_temp);
theta_differ = mean(theta_differ_temp);




for i = 1:path_record_size
    record = path_record{i};
    x = record(:, 2);x = x(1:end-1);
    y = record(:, 3);y = y(1:end-1);
    v = record(:, 4);v = v(1:end-1);
    theta = record(:, 5);theta = theta(1:end-1);
    x_real = x_pos((i - 1) * 50 + 1:i * 50);
    y_real = y_pos((i - 1) * 50 + 1:i * 50);
    v_real = v_pos((i - 1) * 50 + 1:i * 50);
    theta_real = theta_pos((i - 1) * 50 + 1:i * 50);
    x_differ_temp((i-1)*50 + 1:i*50) = abs(x - x_real);
    y_differ_temp((i-1)*50 + 1:i*50) = abs(y - y_real);
    v_differ_temp((i-1)*50 + 1:i*50) = abs(v - v_real);
    theta_differ_temp((i-1)*50 + 1:i*50) = abs(theta - theta_real);  
end
x_differ = mean(x_differ_temp)
y_differ = mean(y_differ_temp)
v_differ = mean(v_differ_temp)
theta_differ = mean(theta_differ_temp)

%% evluate from the new scenario
% load('ps1.mat'); load('ps2.mat'); 
load('net.mat');
global ps1 ps2 net
global initial_i
data_eva = cell(200, 1);
initial_state_eva = cell(200, 1);
for i = 1:200
    obs1_v = rand() * 20 + 20;
    obs1_a = rand() * (-2) + 1;
    ego_v = obs1_v * (1 + (-rand() * (0.2 ) + 0.1));
    obs2_v = rand() * 20 + 20;
    obs2_a = rand() * (-2) + 1;
    obs3_v = obs2_v * (1 + (-rand() * (0.2) + 0.1));
    obs3_a = rand() * (-2) + 1;
    obs1_y = 3 * ego_v/3.6;
    obs2_y = rand() * 50;
    obs3_y = obs2_y - 3 * obs3_v  / 3.6 - rand() * 100;
    initial_state_eva{i} = [ego_v, obs1_y,obs1_v,obs1_a,obs2_y,obs2_v,obs2_a,obs3_y,obs3_v,obs3_a];
%     initial_state_eva{i} = [33.5762   27.9801   33.9624   -0.4136   10.5493   21.3813    0.3427  -93.9939   20.8299...
%         0.9499];
    initial_state_eva{i} = [33.1704706096487,27.6420588413739,36.2856965213763,0.512950062550021,23.6644424451365, ...
        26.9996753196962,0.606809499137584,-35.1213468636565,28.3438063229919,-0.232089352293278];
    ego_v = initial_state_eva{i}(1);
    [flag, data_temp] = run(ego_v, initial_state_eva{i}(2:end), 3.75);
    data_eva{i} = data_temp;
    trainedModel.predictFcn(initial_state_eva{i})
    figure;
    y_pos = data_temp(2:4:4 * 51);
    x_pos = data_temp(1:4:4 * 51);
    plot(x_pos, y_pos, '*');
    initial_i = initial_state_eva{i};
    [t, record_temp] = ode45(@ode, [0:0.1:5],[0 0 initial_i(1)/3.6 0]);
    y_pos_nn = record_temp(:, 2);
    x_pos_nn = record_temp(:, 1);
    hold on;
    plot(x_pos_nn, y_pos_nn, '*');
    
end