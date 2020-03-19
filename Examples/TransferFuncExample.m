% this function is for transfer function example
% written by Shangyan
% 1/29/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

%% Transfer Func
num = [1 0];
den = [1 14 56 160];
% compute state space
[A,B,C,D] = tf2ss(num,den);
% compute transfer function
[num1,den1] = ss2tf(A,B,C,D);

%% SS propagation
t_init = 0;
t_step = 0.1;
t_final = 30;
tspan = t_init:t_step:t_final;
x0 = [0.1;0.5;0.1];
[t_vec,S_vec] = ode45(@(t,x) SSfunc(t,x,A),tspan,x0);

%% Results
y_vec = C*transpose(S_vec);
figure
plot(t_vec,y_vec,'linewidth',1.5)
xlabel('Time(s)')
ylabel('y')
set(gca,'Fontsize',18,'FontWeight','bold')


























