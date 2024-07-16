clc
clear
close all;

load('HW5_question2.mat')

Ts=0.5;
u_id=Z1.u;
y_id=Z1.y;

t=0.5:0.5:500;

%% armax on data Z1
na=3;nb=3;nc=3;
model=armax([y_id u_id],[na nb nc 1]);
aa=model.A;
bb=model.B;

G1 = tf(bb, aa, Ts);

gain_dc=dcgain(G1);
%% armax on data Z3
u_id2=Z3.u;
y_id2=Z3.y;

model2=armax([y_id2 u_id2],[na nb nc 1]);
aa2=model2.A;
bb2=model2.B;

G2 = tf(bb2, aa2, Ts);
gain_dc2=dcgain(G2);

K=(gain_dc-gain_dc2)/(gain_dc*gain_dc2);

G3=G1/(1+K*G1);

[y_hat, t_hat1] = lsim(G3, u_id2, t);

%% validation of identification
SSE = sum((y_id2 - y_hat).^2);
MSE = mean((y_id2 - y_hat).^2);
mean_y_real = mean(y_id2);
SSR = sum((y_hat - mean_y_real).^2);
R_squared = 1 - SSE / SSR;

disp(['SSE: ', num2str(SSE)]);
disp(['MSE: ', num2str(MSE)]);
disp(['R-squared: ', num2str(R_squared)]);

figure;
plot( t,y_id2, 'b',t, y_hat, 'r');
legend('real','estimated');
xlabel('time');
title("transfer function with estemated K")