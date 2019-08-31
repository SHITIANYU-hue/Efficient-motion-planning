function [flag_output, output_args, IDM_output] = run_old(ego, obs, lane_width )
    output_args = [];
    IDM_output = [];
    flag_output = true;
    close all;
    
    %% initial condition
    is_exist_1 = true;
    is_exist_2 = true;
    is_exist_3 = true;
    if obs(1) == -99999
        is_exist_1 = false;
    end
    if obs(4) == -99999
        is_exist_2 = false;
    end
    if obs(7) == -99999
        is_exist_3 = false;
    end
    width = 1.7;
    length = 4.4;
    omega_max = 30/180*pi/2;
    omega_min = -30/180*pi/2;
    a_max = 6 ;
    a_min = -6;
    v_max = 60/3.6;
    v_min = 0/3.6;
    epsilon = lane_width/15;

    nx = 4;    % number of state
    nu = 2;    % number of control
    tf = 5;    % time span
    N = 50;   % number of discreet point
    dt = tf/N ;
    dt_value = linspace(0,tf,N + 1)';

    x_0 = 0;
    y_0 = 0;
    v_0 = ego/3.6;
    theta_0 = 0;

    yf = lane_width;
    thetaf = 0;

    x_obstacle = repmat(obs(1), N+1, 1) + obs(2)/3.6 * linspace(0, tf, N+1)' +...
                 0.5 * obs(3) * linspace(0, tf, N+1)'.^2;
    y_obstacle = zeros(N + 1, 1);
    theta_obstacle = ones(51,1) * 0/180*pi;

    x1_obstacle = repmat(obs(4), N+1, 1) + obs(5)/3.6 * linspace(0, tf, N+1)' +...
                 0.5 * obs(6) * linspace(0, tf, N+1)'.^2;
    y1_obstacle = ones(N + 1, 1) * lane_width;
    theta_obstacle1 = ones(51, 1) * 0/180*pi;

    x2_obstacle = obs(7);
    x2_obstacle_ = repmat(obs(7), N+1, 1);
    v2_obstacle = obs(8)/3.6;
    y2_obstacle = ones(N + 1, 1) * lane_width;
    theta_obstacle2 = ones(N + 1, 1) * 0/180 * pi;

    %% color of line
    color = 'b';
    color1 = strcat(color,'--');
    color2 = strcat(color,'-');
    lineW = 1.5;  

    %%  defination of something
    ele = nx*(N+1) + nu*N + N - 1;
    A =zeros(nx*N,ele);
    b =zeros(nx*N,1);

    %% solvers setting
    ops = sdpsettings('solver','mosek','verbose',1);

    %% solve
    flag = 1;   
    iter = 1; 
    c1 = zeros(1,ele);
    for i=1:1:N
        c1((N+1)*nx+(i-1)*nu+1) = dt; % a
    end
    c2 = zeros(1,ele);
    for i=1:1:N
        c2((N+1)*nx+(i-1)*nu+nu) = dt; % omega
    end
    c3 = zeros(1,ele);
    for i=1:1:N
        c3((i-1)*nx+2) = dt; % y
    end
    c4 = zeros(1, ele);
    for i=1:1:N - 1
        c4((N+1)*nx + N*nu + i) = dt; % y
    end
     x = zeros(ele,1);
     
%% initilize the state (N + 1) * nx
    ts = 0;
    te = tf;

    T = [ts^5    ts^4    ts^3   ts^2   ts^1 1;
         5*ts^4  4*ts^3  3*ts^2 2*ts^1 1    0;
         20*ts^3 12*ts^2 6*ts^1 2      0    0;
         te^5    te^4    te^3   te^2   te^1 1;
         5*te^4  4*te^3  3*te^2 2*te^1 1    0;
         20*te^3 12*te^2 6*te^1 2      0    0
        ];
    %第一列初始位置，第二列初始速度，第三列初始加速度；第四列末端位置；第四列末端速度；第四列末端加速度
    X = [0 v_0 0 (v_0 + obs(5)/3.6)*tf/2 obs(5)/3.6 0];
    Y = [0 0 0  lane_width 0 0];

    A = linsolve(T, X')
    B = linsolve(T, Y')

    t = ts:(te-ts)/50:te;

    x_ini = A(1)*t.^5 + A(2)*t.^4 + A(3)*t.^3 + A(4)*t.^2 + A(5)*t + A(6)*1;
    y_ini = B(1)*t.^5 + B(2)*t.^4 + B(3)*t.^3 + B(4)*t.^2 + B(5)*t + B(6)*1;
    x(1:4:(N+1)*nx) = ego/3.6 * dt * (0:1:N)';
    x(2:4:(N+1)*nx) = y_ini';
    x(3:4:(N+1)*nx) = linspace(v_0, obs(5)/3.6, N + 1);
    x(4:4:(N+1)*nx) = repmat(90/180 * pi - atan(x_ini(end)/y_ini(end)), N + 1, 1);
%     x(1:4:(N+1)*nx) = ego/3.6 * dt * (0:1:N)';
%     x(2:4:(N+1)*nx) = linspace(0,lane_width,N + 1);
%     x(3:4:(N+1)*nx) = repmat(ego/3.6, N + 1, 1);
%     x(4:4:(N+1)*nx) = repmat(90/180 * pi - atan(30/3.6/0.375), N + 1, 1);
     
%% while loop
     while flag 
        fprintf('================== iter=%d ====================\n',iter);
        if iter > 10
            flag_output = false;
            break;
        end
        % ===== define the optimization variables  ===== %%
         y = sdpvar(ele,1); 

        % ===== include all the constraints in F ==== %%
        A =zeros(nx * (N+1) + N - 1, ele);
        b = zeros(nx * (N+1) + N - 1, 1);
        for i=1:N
           x_i = x((i-1)*nx+1);  
           y_i = x((i-1)*nx+2); 
           v_i = x((i-1)*nx+3);   
           theta_i  = x((i-1)*nx+4);
           x_i_ = x((i)*nx+1);  
           y_i_ = x((i)*nx+2); 
           v_i_ = x((i)*nx+3);   
           theta_i_  = x((i)*nx+4);
           AA1 = [-1 0 -cos(theta_i)*dt/2 v_i*sin(theta_i)*dt/2;
                  0 -1 -sin(theta_i)*dt/2 -v_i*cos(theta_i)*dt/2;
                  0 0 -1 0;
                  0 0 0 -1];
           AA2 = [1 0 -dt/2*cos(theta_i_) dt/2*v_i_*sin(theta_i_);
                  0 1 -dt/2*sin(theta_i_) -dt/2*v_i_*cos(theta_i_);
                  0 0 1 0;
                  0 0 0 1];
           BB1 = [0 0;0 0;-dt 0;0 -dt];
           CC1 = [dt/2*v_i_*sin(theta_i_)*theta_i_; -dt/2*v_i_*cos(theta_i_)*theta_i_;0;0] + ...
                 [dt/2*v_i*sin(theta_i)*theta_i; -dt/2*v_i*cos(theta_i)*theta_i;0;0];
           A(i*nx+1:(i+1)*nx,(i-1)*nx+1:i*nx) = AA1;
           A(i*nx+1:(i+1)*nx,i*nx+1:(i+1)*nx) = AA2;
           A(i*nx+1:(i+1)*nx,(N+1)*nx+(i-1)*nu+1:(N+1)*nx+i*nu) = BB1;
    %        A(i*nx+1:(i+1)*nx,(N+1)*nx+i*nu+1:(N+1)*nx+(i+1)*nu) = BB1;
           b(i*nx+1:(i+1)*nx) = CC1;
           if i < N
               A((N+1) * nx + i, [(N+1)*nx + (i-1)*nu + 1,(N+1)*nx + i*nu + 1, (N+1)*nx+N*nu+i]) = [-1,1,-1];
               b((N+1) * nx + i) = 0;
           end
        end
        F = [ A*y == b ];
        % reinitialize x2_obstacle
        x2_obstacle = obs(7);
        for i = 1:1:N
            F = F + [y((i-1)*nx+3) <= v_max, y((i-1)*nx+3) >= v_min];
            F = F + [y((N+1)*nx+(i-1)*nu+1) <= a_max, y((N+1)*nx+(i-1)*nu+1) >= a_min];
            F = F + [y((N+1)*nx+(i-1)*nu+2) <= omega_max, y((N+1)*nx+(i-1)*nu+2) >= omega_min];
            F = F + [y((i-1)*nx+2) <= lane_width + epsilon, y((i-1)*nx+2) >= 0 - epsilon];
            % IDM for x2_obstacle
            if x((i-1)*nx+2) + 1.7/2 <= lane_width/2 % ego car is in the ego lane
                % the following car is following the leading car
                v0 = 30/3.6;
                v1_obstacle = obs(5)/3.6 + obs(6) * dt * (N - 1);
                sn = x1_obstacle(i) - x2_obstacle - length;
                s0 = 2; % minimum distance
                T = 1.5;  % safe time headway
                a = a_max;  % maximum acceleration
                b = 1.67;  % comfortable deceleration
                s_ = s0 + v2_obstacle*T + v2_obstacle * (v2_obstacle - v1_obstacle) / (2 * sqrt(a * b));
                a2_obstacle = a * (1 - (v2_obstacle / v0)^4 - (s_/sn)^2);
                v2_obstacle = a2_obstacle * dt + v2_obstacle;
                x2_obstacle = v2_obstacle * dt + x2_obstacle;
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
            x2_obstacle_(i) = x2_obstacle;
            if is_exist_1 == true
                F = F + implies(y((i-1)*nx+2) - 1.7/2 <= lane_width/2, y((i-1)*nx+1) + 3 * y((i-1)*nx+3) <= x_obstacle(i));
            end
            if is_exist_2 == true
                F = F + implies(-y((i-1)*nx+2) - 1.7/2 <= -lane_width/2, y((i-1)*nx+1) + 3 * y((i-1)*nx + 3) <= x1_obstacle(i));
            end
            if is_exist_3 == true
                F = F + implies(-y((i-1)*nx+2) - 1.7/2 <= -lane_width/2, y((i-1)*nx+1) >= x2_obstacle);
            end
                %         F = F + implies(-y((i-1)*nx+2) <= -1.5, y((i-1)*nx+1) + 1 <= x3_obstacle(i) | y((i-1)*nx+1) >= 1 + x3_obstacle(i));
    %         F = F + implies(-y((i-1)*nx+2) <= -1.5, y((i-1)*nx+1) + 1 <= x4_obstacle(i) | y((i-1)*nx+1) >= 1 + x4_obstacle(i));
            F = F + implies(y((i-1)*nx+2) - 1.7/2 <= 1.5, y((i-1)*nx+1) <= y((i)*nx+1));
        end
        i = N + 1;
        F = F + [y((i-1)*nx+3) <= v_max, y((i-1)*nx+3) >= v_min];
        F = F + [y((i-1)*nx+2) <= lane_width + epsilon, y((i-1)*nx+2) >= 0 - epsilon];
        % IDM
        if x((i-1)*nx+2) + 1.7/2 <= lane_width/2
            v0 = 30/3.6;
            v1_obstacle = obs(5)/3.6 + obs(6) * dt * (N - 1);
            sn = x1_obstacle(i) - x2_obstacle - length;
            s0 = 2; % minimum distance
            T = 1.5;  % safe time headway
            a = a_max;  % maximum acceleration
            b = 1.67;  % comfortable deceleration
            s_ = s0 + v2_obstacle*T + v2_obstacle * (v2_obstacle - v1_obstacle) / (2 * sqrt(a * b));
            a2_obstacle = a * (1 - (v2_obstacle / v0)^4 - (s_/sn)^2);
            v2_obstacle = a2_obstacle * dt + v2_obstacle;
            x2_obstacle = v2_obstacle * dt + x2_obstacle;
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
        x2_obstacle_(i) = x2_obstacle;
        if is_exist_1 == true
            F = F + implies(y((i-1)*nx+2) - 1.7/2 <= lane_width/2, y((i-1)*nx+1) + 3 * y((i-1)*nx+3) <= x_obstacle(i));
        end
        if is_exist_2 == true
            F = F + implies(-y((i-1)*nx+2) - 1.7/2 <= -lane_width/2, y((i-1)*nx+1) + 3 * y((i-1)*nx+3) <= x1_obstacle(i));
        end
        if is_exist_3 == true
            F = F + implies(-y((i-1)*nx+2) - 1.7/2 <= -lane_width/2, x2_obstacle <= y((i-1)*nx+1));
        end
        F = F + [y(1) == x_0];
        F = F + [y(2) == y_0];
        F = F + [y(3) == ego/3.6];
        F = F + [y(4) == 0];
        % no terminal contraint
        F = F + [y(nx*N + 2) == yf];
        coe1 = 0.5;
        coe2 = 1;
        coe3 = 0.005;
        coe4 = 100;
        if is_exist_2 == false
            coe3 = 0;
        end
        Temp = zeros(size(y,1), 1);
        index = 2:nx:nx*(N + 1);
        Temp(index) = lane_width;
        solvesdp(F, coe1*c1*y.^2 + coe2*c3*(y - Temp).^2 ...
            + coe3 * (y(N*nx + 1) - x1_obstacle(end)).^2 ...
             + coe4 * c4 * y.^2, ops);

        y = double(y);
        output_args = y;
        x_pos = y(1:nx:(N+1)*nx);
        y_pos = y(2:nx:(N+1)*nx);
        v_pos =  y(3:nx:(N+1)*nx);
        theta_pos = y(4:nx:(N+1)*nx);
        u1 = y(nx*(N+1)+1:nu:nx*(N+1)+nu*N-1);
        u2 =  y(nx*(N+1)+2:nu:nx*(N+1)+nu*N);
%         figure(91)
%             plot(x_pos, y_pos, color1, 'LineWidth', lineW); grid on; hold on;
%             xlabel('x position (m)'); ylabel('y position (m)');
%             plot(x_pos, y_pos, '*');
%             hold on;
%             plot(x_obstacle', zeros(N + 1, 1), '*');
%             hold on;
%             plot(x1_obstacle', ones(N + 1, 1) * lane_width, '*');
       if iter >= 2 
            x_differ = zeros(N+1,1);
            y_differ = zeros(N+1,1);
            for i=1:N

               x_differ(i) = ( y((i-1)*nx+1) - x((i-1)*nx+1) );
               y_differ(i) = ( y((i-1)*nx+2) - x((i-1)*nx+2) );
            end



            fprintf('Maximum difference between two successive solutions are:\n');
            fprintf('x: %5f, y:%5f \n', max(abs(x_differ)), max(abs(y_differ)) );

            if max(abs(x_differ)) <= 0.1 && max(abs(y_differ)) <= 0.1 
               flag = 0; 
            else 
               x = y;
               iter = iter + 1; 
            end
        end 
        if iter == 1
           x = y;
           iter = iter + 1;
        end
     end

     %% plot
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
%             obstacle_x_temp = x2_obstacle_(t_num) - 4.4/2;
%             obstacle_y_temp = y2_obstacle(t_num) - 1.7/2;
%             x = [obstacle_x_temp, obstacle_x_temp+length, obstacle_x_temp+length, obstacle_x_temp, obstacle_x_temp];
%             y = [obstacle_y_temp, obstacle_y_temp, obstacle_y_temp+width, obstacle_y_temp+width, obstacle_y_temp];
%             hold on;
%             l4 = line(x,y, 'color', 'c');
%         end
%         pause(.05);
%         if (t_num < N + 1)
% %             delete(l1);
%             if is_exist_1 == true
% %                 delete(l2);
%             end
%             if is_exist_2 == true
% %                 delete(l3);
%             end
%             if is_exist_3 == true
% %                 delete(l4);
%             end
%         end
%     end

    % plot v, theta, a, omega information

    %     figure;
    %     t = linspace(0, tf, 51);
    %     plot(v_pos, '-.');

    %    
    IDM_output = x2_obstacle_;
    
end