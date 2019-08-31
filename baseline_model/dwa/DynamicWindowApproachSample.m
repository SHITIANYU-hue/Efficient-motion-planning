%时间已过 216.508322 秒。
function [] = DynamicWindowApproachSample()  
   
% close all;  
% clear all;  
   
disp('Dynamic Window Approach sample program start!!')  
  
% x1=[0 0 0 0 0]';% 机器人1的初期状态[x(m),y(m),yaw(Rad),v(m/s),w(rad/s)]
x2=[2 0 0 0 0]';% 机器人2的初期状态[x(m),y(m),yaw(Rad),v(m/s),w(rad/s)]
x3=[4 2 0 0 0]';% 机器人3的初期状态[x(m),y(m),yaw(Rad),v(m/s),w(rad/s)]
x4=[0 0 0 2 0]';% 机器人4的初期状态[x(m),y(m),yaw(Rad),v(m/s),w(rad/s)]
% x5=[-6 0 0 0 0]';% 机器人5的初期状态[x(m),y(m),yaw(Rad),v(m/s),w(rad/s)]
% goal1=[3,0, 0,2];% 目标点1位置 [x(m),y(m),yaw,v]  
goal2=[8,0, 0,2];% 目标点2位置 [x(m),y(m)]
goal3=[8,2,0,2];% 目标点2位置 [x(m),y(m)]
goal4=[10,3.75,0,6];% 目标点2位置 [x(m),y(m)]
% goal5=[-5,1];% 目标点2位置 [x(m),y(m)]
% 障碍物位置列表 [x(m) y(m)]  
% obstacle=[0 2;  
%           4 2;  
%           4 4;  
%           5 4;  
%           5 5;  
%           5 6;  
%           5 9  
%           8 8  
%           8 9  
%           7 9];  
%更新状态矩阵 
obstacle=[10 10];  
obstacleR=0.1;% 冲突判定用的障碍物半径  
global dt; dt=0.1;% 时间[s]  
  
% 机器人运动学模型  
% 最高速度m/s],最高旋转速度[rad/s],加速度[m/ss],旋转加速度[rad/ss],  
% 速度分辨率[m/s],转速分辨率[rad/s]]  
Kinematic=[5.0,toRadian(20.0),0.5,toRadian(50.0),0.01,toRadian(1)];  
% 评价函数参数 [heading,dist,velocity,predictDT]  
evalParam=[12,0.2,0.1,3.0];  
area=[0 45 0 4];% 模拟区域范围 [xmin xmax ymin ymax]  
  
% 模拟实验的结果  
result.x1=[];
result.x2=[];
result.x3=[];
result.x4=[];
result.x5=[];
tic;  
% movcount=0;  
% Main loop  
for i=1:5000  
    t1=clock;
    % 机器人1DWA参数输入  
%     [u,traj]=DynamicWindowApproach(x1,Kinematic,goal1,evalParam,obstacle,obstacleR); 
%     x1=f(x1,u);% 机器人1移动到下一个时刻   
    [u,traj]=DynamicWindowApproach(x2,Kinematic,goal2,evalParam,obstacle,obstacleR);
     x2=f(x2,u);% 机器人2移动到下一个时刻  
  [u3,traj]=DynamicWindowApproach(x3,Kinematic,goal3,evalParam,obstacle,obstacleR);
     x3=f(x3,u3);% 机器人3移动到下一个时刻    
      [u,traj]=DynamicWindowApproach(x4,Kinematic,goal4,evalParam,obstacle,obstacleR);
     x4=f(x4,u);% 机器人4移动到下一个时刻  
%     [u,traj]=DynamicWindowApproach(x5,Kinematic,goal5,evalParam,obstacle,obstacleR);
%      x5=f(x5,u);% 机器人5移动到下一个时刻     
    % 模拟结果的保存  
%     result.x1=[result.x1; x1'];  
    result.x2=[result.x2; x2'];
   result.x3=[result.x3; x3'];
    result.x4=[result.x4; x4'];
%    result.x5=[result.x5; x5'];
    % 是否到达目的地  
    if  norm(x2(1:4)-goal2')<0.5&& norm(x3(1:4)-goal3')<0.5&& norm(x4(1:4)-goal4')<0.5
        disp('Arrive Goal!!');break;  
    end   
    %====Animation====  
    hold off;  
    ArrowLength=0.5;%   
    % 机器人  
%     quiver(x1(1),x1(2),ArrowLength*cos(x1(3)),ArrowLength*sin(x1(3)),'ok');hold on; 
%     plot(result.x1(:,1),result.x1(:,2),'-b');hold on;  
%     plot(goal1(1),goal1(2),'*r');hold on;  
%     quiver(x2(1),x2(2),ArrowLength*cos(x2(3)),ArrowLength*sin(x2(3)),'ok');hold on;
%      plot(result.x2(:,1),result.x2(:,2),'-b');hold on;  
% %     plot(goal2(1),goal2(2),'*r');hold on;  
%     quiver(x3(1),x3(2),ArrowLength*cos(x3(3)),ArrowLength*sin(x3(3)),'ok');hold on;
%        plot(result.x3(:,1),result.x3(:,2),'-b');hold on;  
%     plot(goal3(1),goal3(2),'*r');hold on;  
    quiver(x4(1),x4(2),ArrowLength*cos(x4(3)),ArrowLength*sin(x4(3)),'ok');hold on;
       plot(result.x4(:,1),result.x4(:,2),'-b');hold on;  
%     plot(goal4(1),goal4(2),'*r');hold on;  
%     quiver(x5(1),x5(2),ArrowLength*cos(x5(3)),ArrowLength*sin(x5(3)),'ok');hold on;
%       plot(result.x5(:,1),result.x5(:,2),'-b');hold on;  
%     plot(goal5(1),goal5(2),'*r');hold on;  
    plot(obstacle(:,1),obstacle(:,2),'*k');hold on;  
    % 探索轨迹  
    if ~isempty(traj)  
        for it=1:length(traj(:,1))/5  
            ind=1+(it-1)*5;  
            plot(traj(ind,:),traj(ind+1,:),'-g');hold on;  
        end  
    end  
    axis(area);  
    grid on;  
    drawnow;  
    %movcount=movcount+1;  
    %mov(movcount) = getframe(gcf);%   
end  
t1=clock;
estime(t2,t1)
toc  
%movie2avi(mov,'movie.avi');  
   
  
function [u,trajDB]=DynamicWindowApproach(x,model,goal,evalParam,ob,R)  
  
% Dynamic Window [vmin,vmax,wmin,wmax]  
Vr=CalcDynamicWindow(x,model);  
  
% 评价函数的计算  
[evalDB,trajDB]=Evaluation(x,Vr,goal,ob,R,model,evalParam);  
  
if isempty(evalDB)  
    disp('no path to goal!!');  
    u=[0;0];return;  
end  
  
% 各评价函数正则化  
evalDB=NormalizeEval(evalDB);  
  
% 最终评价函数的计算  
feval=[];  
for id=1:length(evalDB(:,1))  
    feval=[feval;evalParam(1:3)*evalDB(id,3:5)'];  
end  
evalDB=[evalDB feval];  
  
[maxv,ind]=max(feval);% 最优评价函数  
u=evalDB(ind,1:2)';%   

function [evalDB,trajDB]=Evaluation(x,Vr,goal,ob,R,model,evalParam)  
%   
evalDB=[];  
trajDB=[];  
for vt=Vr(1):model(5):Vr(2)  
    for ot=Vr(3):model(6):Vr(4)  
        % 轨迹推测; 得到 xt: 机器人向前运动后的预测位姿; traj: 当前时刻 到 预测时刻之间的轨迹  
        [xt,traj]=GenerateTrajectory(x,vt,ot,evalParam(4),model);  %evalParam(4),前向模拟时间;  
        % 各评价函数的计算  
        heading=CalcHeadingEval(xt,goal);  
        dist=CalcDistEval(xt,ob,R);  
        vel=abs(vt);  
        % 制动距离的计算  
        stopDist=CalcBreakingDist(vel,model);  
        if dist>stopDist %   
            evalDB=[evalDB;[vt ot heading dist vel]];  
            trajDB=[trajDB;traj];  
        end  
    end  
end  
  
function EvalDB=NormalizeEval(EvalDB)  
% 评价函数正则化  
if sum(EvalDB(:,3))~=0  
    EvalDB(:,3)=EvalDB(:,3)/sum(EvalDB(:,3));  
end  
if sum(EvalDB(:,4))~=0  
    EvalDB(:,4)=EvalDB(:,4)/sum(EvalDB(:,4));  
end  
if sum(EvalDB(:,5))~=0  
    EvalDB(:,5)=EvalDB(:,5)/sum(EvalDB(:,5));  
end  
  
function [x,traj]=GenerateTrajectory(x,vt,ot,evaldt,model)  
% 轨迹生成函数  
% evaldt：前向模拟时间; vt、ot当前速度和角速度;   
global dt;  
time=0;  
u=[vt;ot];% 输入值  
traj=x;% 机器人轨迹  
while time<=evaldt  
    time=time+dt;% 时间更新  
    x=f(x,u);% 运动更新  
    traj=[traj x];  
end  
  
function stopDist=CalcBreakingDist(vel,model)  
% 根据运动学模型计算制动距离,这个制动距离并没有考虑旋转速度，不精确吧！！！  
global dt;  
stopDist=0;  
while vel>0  
    stopDist=stopDist+vel*dt;% 制动距离的计算  
    vel=vel-model(3)*dt;%   
end  
  
function dist=CalcDistEval(x,ob,R)  
% 障碍物距离评价函数  
  
dist=100;  
for io=1:length(ob(:,1))  
    disttmp=norm(ob(io,:)-x(1:2)')-R;%
    if dist>disttmp% 离障碍物最小的距离  
        dist=disttmp;  
    end  
end  
  
% 障碍物距离评价限定一个最大值，如果不设定，一旦一条轨迹没有障碍物，将太占比重  
if dist>=2*R  
    dist=2*R;  
end  
  
function heading=CalcHeadingEval(x,goal)  
% heading的评价函数计算  
  
theta=toDegree(x(3));% 机器人朝向  
goalTheta=toDegree(atan2(goal(2)-x(2),goal(1)-x(1)));% 目标点的方位  
  
if goalTheta>theta  
    targetTheta=goalTheta-theta;% [deg]  
else  
    targetTheta=theta-goalTheta;% [deg]  
end  
  
heading=180-targetTheta;  
  
function Vr=CalcDynamicWindow(x,model)  
%  
global dt;  
% 车子速度的最大最小范围  
Vs=[0 model(1) -model(2) model(2)];  
  
% 根据当前速度以及加速度限制计算的动态窗口  
Vd=[x(4)-model(3)*dt x(4)+model(3)*dt x(5)-model(4)*dt x(5)+model(4)*dt];  
  
% 最终的Dynamic Window  
Vtmp=[Vs;Vd];  
Vr=[max(Vtmp(:,1)) min(Vtmp(:,2)) max(Vtmp(:,3)) min(Vtmp(:,4))];  
  
function x = f(x, u)  
% Motion Model  
% u = [vt; wt];当前时刻的速度、角速度  
global dt;  
   
F = [1 0 0 0 0  
     0 1 0 0 0  
     0 0 1 0 0  
     0 0 0 0 0  
     0 0 0 0 0];  
   
B = [dt*cos(x(3)) 0  
    dt*sin(x(3)) 0  
    0 dt  
    1 0  
    0 1];  
  
x= F*x+B*u;  
  
function radian = toRadian(degree)  
% degree to radian  
radian = degree/180*pi;  
  
function degree = toDegree(radian)  
% radian to degree  
degree = radian/pi*180; 