data = cell(100000, 1);
initial_state = cell(100000, 1);
IDM = cell(100000, 1);
i = 1;
% load('trainedModel.mat');
while (1)
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
    initial = [ego_v, obs1_y,obs1_v,obs1_a,obs2_y,obs2_v,obs2_a,obs3_y,obs3_v,obs3_a];
    ego_v = initial(1);
%     if (trainedModel.predictFcn(initial_state{i}) == -1) 
        try
            [flag, data_temp, IDM_output] = run(ego_v, initial(2:end), 3.75);
        catch
            continue
        end
        if flag == true
            data{i} = data_temp;
            IDM{i} = IDM_output;
            initial_state{i} = initial;
            i = i + 1;
        end
%     end

end