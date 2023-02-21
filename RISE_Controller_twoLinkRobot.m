close all; clear all; clc;
% RISE Based Non-Linear controller for 2-Link rigid, revolute, direct drive
% robot Arm.
% 2 - link IMI arm dynamics
p1       = 3.473;   % kg.m^2
p2       = 0.196;   % kg.m^2
p3       = 0.242;   % kg.m^2
f1       = 5.3;     % Nm.sec
f2       = 1.1;     % Nm.sec
% Constant system parameters
theta    = [p1;p2;p3;f1;f2];
% Simulation time
tf       = 45;
% initial coditions vector X0 = [e0;r0;thetahat0];
X0       = 2*[ones(6,1)];
% ode options 
options  = odeset('RelTol',1e-3,'AbsTol',1e-3);

% Integration of Dynamic system model
[t,X]     = ode45(@(t,X) twoLinkdynamics(t,X,theta),[0 tf],X0,options);
% initializing the input vector
u         =   zeros(2,length(t));
for i     =   1:length(t)
    u(:,i)  =   getcontrol(t(i),X(i,:),theta);
end
% Desired trajectory
for i   =       1:length(t)
    t1 = t(i);
    qd(:,i)   = [(30*(pi/180)*sin(1.5*t1)+20*(pi/180))*(1-exp(-0.01*t1^3)); -(20*(pi/180)*sin(0.5*t1)+10*(pi/180))*(1-exp(-0.01*t1^3))];
    Td(:,i)   = [0.5*cos(0.5*t1);sin(t1)];

end
% Computing the error 
e           =       X(:,1:2)';
r           =       X(:,3:4)';
Usign       =       X(:,5:6)';
% Computing the actual link parameter(Trajectory) q
q           =       qd-e;

% Trajectory Plots ( Actual, Desired Trajector V/s Time)
figure(1)
plot(t,qd,'-','LineWidth',1)
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(t,q,':','LineWidth',1)
title('Actual and Desired Trajectories')
xlabel('time') % x-axis label
ylabel('Trajectories') % y-axis label
legend('qd1','qd2','q1','q2');
hold off

% Tracking error plot
figure(2)
plot(t,e,'--','LineWidth',2)
title('Tracking Error')
xlabel('time') % x-axis label
ylabel('Errors') % y-axis label
legend('e(1)','e(2)')

% Control Input Plot
figure(3)
plot(t,u(1,:),'-',t,u(2,:),':','LineWidth',1)
title('Control Input')
xlabel('time') % x-axis label
ylabel('Torques') % y-axis label
legend('u1','u2');

%Disturbance Vector plot
figure(4)
plot(t,Td(1,:),'-',t,Td(2,:),':','LineWidth',1)
title('Disturbance signal')
xlabel('time') % x-axis label
ylabel('distudbance') % y-axis label
legend('\tau_d1','\tau_d2')

% Filtered Tracking error plot
figure(5)
plot(t,r,'-','LineWidth',2)
title('Tracking Error')
xlabel('time') % x-axis label
ylabel('Errors') % y-axis label
legend('r1','r2')

% Filtered Tracking error plot
figure(6)
plot(t,Usign,'-','LineWidth',2)
title('Parse state Usign')
xlabel('time') % x-axis label
ylabel('Errors') % y-axis label
legend('Usign1','Usign2')

function [XDot] = twoLinkdynamics(t,X,theta)
global e2i
% Parse parameter vector
p1 = theta(1);
p2 = theta(2);
p3 = theta(3);
f1 = theta(4);
f2 = theta(5);

% Controller Gains tuning
K        = 50; 
a1       = 2; 
a2       = 2;
B        = 3;

% Desired trajectory and needed derivatives
[qd,qdDot,qdDotDot] = des_trajectory(t);
%  Parse State Computing 
e        = [X(1);X(2)];
e2       = [X(3);X(4)];
Usign    = [X(5);X(6)];

% Compute current q and qDot for convenience
q        = qd-e;
qDot     = -e2 + a1*e + qdDot;

% Compute cos(x2) and sin(x2) for convenience
c2       = cos(q(2));
s2       = sin(q(2));

% Compute current matrices for the dynamics
M        = [p1 + 2*p3*c2 p2 + p3*c2;p2 + p3*c2 p2];
V        = [-p3*s2*qDot(2) -p3*s2*(qDot(1) + qDot(2));p3*s2*qDot(1) 0];
fd       = [f1 0;0 f2];
Td       = [0.5*cos(0.5*t);2*sin(t)];

S = V*qDot+fd*q+M*qdDot+a1*M*e2-a1^2*M*e2+M*a2*e2;

% Design controller
if (t==0)
    e2i = e2;
end
u         =(K+1)*(e2-e2i)+Usign; 
r         = M\(S-u+Td);
%udot=(K+1)*r+B*sign(e2); 
% Compute current closed-loop dynamics
eDot       = e2 - a1*e;
e2Dot      = r - a2*e2;
UsignDot    = (K+1)*r+B*sign(e2);

% Stacked dynamics vector (XDot is the same size and "form" as X)
XDot        = [eDot;e2Dot;UsignDot];
end
% Desired Trajectory function
function [qd,qd_dot,qd_dot_dot] = des_trajectory(t)

qd = [(30*(pi/180)*sin(1.5*t)+20*(pi/180))*(1-exp(-0.01*t^3)); -(20*(pi/180)*sin(0.5*t)+10*(pi/180))*(1-exp(-0.01*t^3))];
qd_dot = [    (3*t^2*((6283*sin((3*t)/2))/12000 + 6283/18000))/(100*exp(t^3/100)) - (6283*cos((3*t)/2)*(1/exp(t^3/100) - 1))/8000;
        (6283*cos(t/2)*(1/exp(t^3/100) - 1))/36000 - (3*t^2*((6283*sin(t/2))/18000 + 6283/36000))/(100*exp(t^3/100)) ];
qd_dot_dot=[(18849*t^2*cos((3*t)/2))/(400000*exp(t^3/100)) + (3*t*(6283/18000+(6283*sin((3*t)/2))/12000))/(50*exp(t^3/100))-(9*t^4*(6283/18000+(6283*sin((3*t)/2))/12000))/(10000*exp(t^3/100))+(18849*(-1+exp(-t^3/100))*sin((3*t)/2))/16000;
    (-6283*t^2*cos(t/2))/(600000*exp(t^3/100))-(3*t*(6283/36000+(6283*sin(t/2))/18000))/(50*exp(t^3/100))+(9*t^4*(6283/36000+(6283*sin(t/2))/18000))/(10000*exp(t^3/100))-(6283*(-1+exp(-t^3/100))*sin(t/2))/72000];
end

function [u] = getcontrol(t,X,theta)
global e2i
% Parse parameter vector
p1 = theta(1);
p2 = theta(2);
p3 = theta(3);
f1 = theta(4);
f2 = theta(5);

% Select gains for controller
K        = 50; 
a1       = 2; 
a2       = 2;
% Desired trajectory and needed derivatives
[qd,qdDot,qdDotDot] = des_trajectory(t);
% Current error [e;r;thetahat])
e           = [X(1);X(2)];
e2          = [X(3);X(4)];
Usign       = [X(5);X(6)];
% Compute current x and xDot for convenience
q           = qd-e;
qDot        = -e2 + a1*e + qdDot;

% Compute cos(q2) and sin(q2) for convenience
c2          = cos(q(2));
s2          = sin(q(2));

% Compute current matrices for the dynamics
M           = [p1 + 2*p3*c2 p2 + p3*c2;p2 + p3*c2 p2];
V           = [-p3*s2*qDot(2) -p3*s2*(qDot(1) + qDot(2));p3*s2*qDot(1) 0];
fd          = [f1 0;0 f2];
Td          = [0.5*cos(0.5*t);sin(t)];

% S = V*qDot+fd*q+M*qdDot+a1*M*e2-a1^2*M*e2+M*a2*e2;
if (t==0)
    e2i = e2;
end
u        = (K+1)*(e2-e2i)+Usign; 
end

