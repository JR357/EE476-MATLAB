% this function is used to simulate WEC
% 3/3/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

%% Parameters
rho = 1000;
g = 9.81;
R = 2;
c = 5000;
A = 1e4;
omega = 2*pi/7;  % rad/s
% store the para
sys.rho = rho;
sys.g = g;
sys.R = R;
sys.c = c;
sys.A = A;
sys.omega = omega;

%% Propagation
t_init = 0;
t_step = 0.1;
t_final = 100;
t_vec = t_init:t_step:t_final;
sys.t_step = t_step;
% initial condition
x0 = [0;0];  
% time domain propagation
counter = 1;
for t = t_init:t_step:t_final
    
    S_vec(:,counter) = x0;
    % prop
    x_p = RK4(t,x0,sys);
    
    % update
    x0 = x_p;
    counter = counter + 1;
end

%% Ode 45
[~,S_vec1] = ode45(@(t,x) SS_WEC(t,x,sys),t_vec,x0);


%% Plot the results
figure
plot(t_vec,S_vec(1,:),'linewidth',1.5)
hold on
plot(t_vec,S_vec1(:,1),'linewidth',1.5)
xlabel('Time(s)')
ylabel('z(m)')
legend('RK4','ODE45')
set(gca,'Fontsize',18,'FontWeight','bold')
figure
plot(t_vec,S_vec(2,:),'linewidth',1.5)
hold on
plot(t_vec,S_vec1(:,2),'linewidth',1.5)
xlabel('Time(s)')
ylabel('v(m/s)')
legend('RK4','ODE45')
set(gca,'Fontsize',18,'FontWeight','bold')
figure
plot(t_vec,S_vec(1,:) - transpose(S_vec1(:,1)),'linewidth',1.5)
xlabel('Time(s)')
ylabel('e_z(m)')
set(gca,'Fontsize',18,'FontWeight','bold')

















