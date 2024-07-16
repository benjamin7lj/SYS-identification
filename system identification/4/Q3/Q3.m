clc
clear
close all;
%% reading data
load('401126125-q3.mat')

y_id=id.y;
u_id=id.u;
Ts=id.Ts;

y_id0=y_id(1:20);
u_id0=u_id(1:20);
tt=0:Ts:499*Ts;

y_val=val.y;
u_val=val.u;
na=3;nb=3;
%% using arx function
model=arx([y_id u_id],[na nb 1]);
aa=model.A;
bb=model.B;
z=tf('z',Ts);
GZ2=(bb(1)*(z^3)+bb(2)*(z^2)+bb(3)*z+bb(4))/...
 (z^3+aa(2)*z^2+aa(3)*z+aa(4));


[y_hat, t_model2] = lsim(GZ2, u_val, tt);

%validation of identification
SSE = sum((y_val - y_hat).^2);
MSE = mean((y_val - y_hat).^2);
mean_y_real = mean(y_val);
SSR = sum((y_hat - mean_y_real).^2);
R_squared = 1 - SSE / SSR;

disp(['SSE: ', num2str(SSE)]);
disp(['MSE: ', num2str(MSE)]);
disp(['R-squared: ', num2str(R_squared)]);

%% IIV

N=length(y_id);
na=3;nb=na;

U = zeros(N, na+nb);
for i = 1:n
    U(i+1:N, i) = -y_id(1:N-i, 1);
end
for i = 1:nb
    U(i:N, n+i) = u_id(1:N-i+1, 1);
end
Z=U;

Ut=U';
tetha_hat_Ls=(inv(Ut*U))*(Ut*y_id);

tetha_hat_1=tetha_hat_Ls;

for t=1:50
    y_hat_0=U*tetha_hat_1;
    for i = 1:n
        Z(i+1:N, i) = -y_hat_0(1:N-i, 1);
    end
    tetha_hat_IIV=(inv(Z'*U)*Z')*y_id;
    tetha_hat_1=tetha_hat_IIV;
end
tetha_hat_IIV=tetha_hat_IIV';
Gz = tf(tetha_hat_IIV(na+1:na+nb),[1, tetha_hat_IIV(1:na)], Ts);
[y_hat_IIV, t_mod] = lsim(Gz, u_val, tt);


%% validation of identification
SSE = sum((y_val - y_hat_IIV).^2);
MSE = mean((y_val - y_hat_IIV).^2);
mean_y_real = mean(y_val);
SSR = sum((y_hat_IIV - mean_y_real).^2);
R_squared = 1 - SSE / SSR;

disp(['SSE_IIV: ', num2str(SSE)]);
disp(['MSE_IIV: ', num2str(MSE)]);
disp(['R-squared_IIV: ', num2str(R_squared)]);

figure;
plot( tt,y_val, 'b',tt, y_hat, 'r',tt, y_hat_IIV, 'g');
legend('real','ARX','IIV');
xlabel('time');

