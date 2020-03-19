function    x_p = RK4(t,x,sys)
% function    [] = RK4()
% this function is used to apply RK propagation
% 3/3/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read the para
t_step = sys.t_step;

%% Propagation
h = t_step;
k1 = h*SS_WEC(t,x,sys);
k2 = h*SS_WEC(t+h/2,x+k1/2,sys);
k3 = h*SS_WEC(t+h/2,x+k2/2,sys);
k4 = h*SS_WEC(t+h,x+k3,sys);
x_p = x + 1/6*(k1 + 2*k2 + 2*k3 + k4);
t_p = t + h;


end









