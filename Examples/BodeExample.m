% this function is used to compute the bode plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

%%
G = tf([10 30],[1 3 4 4 0]);
[A,B,C,D] = tf2ss([10 30],[1 3 4 4 0]);

figure
bode([10 30],[1 3 4 4 0])
figure
bode(G)
figure
bode(A,B,C,D)





































