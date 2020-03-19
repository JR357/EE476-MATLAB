% this function is to explain the linearization 
% written by Shangyan
% 1/29/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
x = 5:0.1:7;
y = 10:0.1:12;

z = x.*y;



%% Results
figure
plot3(x,y,z,'linewidth',1.5)
xlabel('x_1')
ylabel('x_2')
xlim([-6,6])
ylim([-6,6])
set(gca,'Fontsize',18,'FontWeight','bold')


























