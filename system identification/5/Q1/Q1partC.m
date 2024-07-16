clc
clear
close all;

load('q1_401126125.mat')

u_id=u(1:500);
u_val=u(500:1000);
v_id=v(1:500);
v_val=v(500:1000);
y_id=y(1:500);
y_val=y(500:1000);
z_id=z(1:500);
z_val=z(500:1000);
% w_id=wgn(500,1,0);
% w_val=wgn(501,1,0);

tt=0:0.1:49.9;

%% B&J Function for C(z)

t=0:0.1:50;
na=1;nb=na;nc=na;nd=nc;p=na+nb;
model=bj([z_id y_id],[na nb nc nd 1]);
aa=model.F;
bb=model.B;
cc=model.C;
dd=model.D;


G2 = tf(bb,aa, 0.1);
[z_hat2, t_hat] = lsim(G2, y_val, t);
%G3 = tf(cc,dd, 0.1);

z_hat=z_hat2;

%% validation of identification
SSE = sum((z_val - z_hat).^2);
MSE = mean((z_val - z_hat).^2);
mean_z_real = mean(z_val);
SSR = sum((z_hat - mean_z_real).^2);
R_squared = 1 - SSE / SSR;

disp(['SSE: ', num2str(SSE)]);
disp(['MSE: ', num2str(MSE)]);
disp(['R-squared: ', num2str(R_squared)]);


figure;
plot( t,z_val, 'b',t, z_hat, 'r');
legend('real','estimated');
xlabel('time');
title("using B&J model")
