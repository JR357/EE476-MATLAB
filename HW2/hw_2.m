% Script for modeling a ball and plank syste

T = 0.001; % Sampling period
tspan = 0:T:8; % Time span
y0 = [0 0]; % Initial conditions: [position, velocity]

m = .1; % Ball mass
g = -9.8; % Gravitational acceleration
I = 1; % Ball moment of inertia
R = 1; % Ball radius

% *****************
% * Input Functions
% *****************
theta1 = @(t) 5*sin(2*pi/5*t);
theta2 = @(t) t >= 0;

% ************
% * Analytically obtained solution
% ************
% Calculate output with input theta1
input1 = theta1(tspan);
forcedResponse1 = -m*g/(m+I/R^2)*tspan;
forcedInput1 = conv(forcedResponse1, input1)*T;
xp1 = forcedInput1(1:length(tspan)) + y0(1) + y0(2).*tspan;

% Calculate output with input theta2
input2 = theta2(tspan);
forcedResponse2 = -m*g/(m+I/R^2)*tspan;
forcedInput2 = conv(forcedResponse2, input2)*T;
xp2 = forcedInput2(1:length(tspan)) + y0(1) + y0(2).*tspan;

% *****************
% * Plot response
% *****************
figure(1)
subplot(2, 1, 1)
xp1_plt = plot(tspan, xp1, 'linewidth', 2);

xp1_h = gca;
xp1_h.FontSize = 16;
xlabel('Time (s)', 'fontsize', 18);
legend('Position $x_p$','interpreter','latex', ...
       'fontsize', 16, 'Location', 'northwest')
title('(Analytical) Output of system with input \Theta = 5sin(2\pi/5*t)', ...
    'fontsize', 16)
grid on;

subplot(2, 1, 2)
xp2_plt = plot(tspan, xp2, 'linewidth', 2);

xp2_h = gca;
xp2_h.FontSize = 16;
xlabel('Time (s)', 'fontsize', 18);
legend('Position $x_p$','interpreter','latex', ...
       'fontsize', 16, 'Location', 'northwest')
title('(Analytical) Output of system with input \Theta = u(t)', ...
    'fontsize', 16)
grid on;


% *****************
% * Numerical solution
% *****************
[t1, y1] = ode45(@(t,y) odefun(t,y, theta1, m, g, I, R), tspan, y0);
[t2, y2] = ode45(@(t,y) odefun(t,y, theta2, m, g, I, R), tspan, y0);

% *****************
% * Plot response
% *****************
figure(2)
subplot(2, 1, 1)
plot(t1, y1, 'linewidth', 2);

xpa1_h = gca;
xpa1_h.FontSize = 16;
xlabel('Time (s)', 'fontsize', 18);
legend('Position $x_p$', 'Velocity $\dot{x_p}$','interpreter','latex', ...
       'fontsize', 16, 'Location', 'northwest')
   title('(Numerical) Output of system with input \Theta = 5sin(2\pi/5*t)', ...
       'fontsize', 16)
grid on

subplot(2, 1, 2);
plot(t2, y2, 'linewidth', 2);

xpa2_h = gca;
xpa2_h.FontSize = 16;
xlabel('Time (s)', 'fontsize', 18);
legend('Position $x_p$', 'Velocity $\dot{x_p}$','interpreter','latex', ...
       'fontsize', 16, 'Location', 'northwest')
title('(Numerical) Output of system with input \Theta = u(t)', ...
    'fontsize', 16)
grid on

% *****************
% Function for describing state space model
% *****************
function dydt = odefun(t, y, theta, m, g, I, R)
    
    dydt = zeros(2, 1);
    
    dydt(1) = y(2);
    dydt(2) = -(m*g)/(m+I/R^2)*theta(t);
end