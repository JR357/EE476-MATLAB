function    dx = SS_WEC(t,x,sys)
% function    [] = SS_WEC()
% SS model of the WEC 
% 3/3/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
rho = sys.rho;
g = sys.g;
R = sys.R;
c = sys.c;
A = sys.A;
omega = sys.omega;
M = 2/3*rho*pi*R^3;

%% Differential equations
Fe = A*cos(omega*t);
dx = zeros(size(x));
dx(1) = x(2);
dx(2) = 1/M*(Fe + 1/3*rho*g*pi*(x(1)^3 - 3*R^2*x(1)) - c*x(2));



end






