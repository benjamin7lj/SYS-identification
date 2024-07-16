clc
close all;

Ts=4.5;
t=0:Ts:2222*Ts;

figure;
plot(t, y_id, 'b', t, u_id, 'r');
legend('output','input');
xlabel('time');
title("input-output-identification")

figure;
plot(t, y_val, 'b', t, u_val, 'r');
legend('output','input');
xlabel('time');
title("input-output-validation")

tt=0:Ts:(length(u_val)-1)*Ts;

na=6;nb=na;

%% LEAST SQR


N=length(y_id);

U = zeros(N, na+nb);
for i = 1:na
    U(i+1:N, i) = -y_id(1:N-i, 1);
end
for i = 1:nb
    U(i:N, na+i) = u_id(1:N-i+1, 1);
end
Ut=U';
tetha_hat_LS=((inv(Ut*U))*(Ut*y_id))';
Gz = tf(tetha_hat_LS(na+1:na+nb),[1, tetha_hat_LS(1:na)], Ts);
[y_hat, t_mod] = lsim(Gz, u_val, tt);

DET=det(Ut*U);
%% using arx function
model=arx([y_id u_id],[na nb 1]);
aa=model.A;
bb=model.B;

GZ2 = tf(bb,aa, Ts);
[y_hat2, t_model2] = lsim(GZ2, u_val, tt);


%validation of identification
SSE = sum((y_val - y_hat).^2);
MSE = mean((y_val - y_hat).^2);
mean_y_real = mean(y_val);
SSR = sum((y_hat - mean_y_real).^2);
R_squared = 1 - SSE / SSR;

disp(['SSE: ', num2str(SSE)]);
disp(['MSE: ', num2str(MSE)]);
disp(['R-squared: ', num2str(R_squared)]);


figure;
plot( tt,y_val, 'b',tt, y_hat, 'r',tt, y_hat2, 'g');
legend('real','estimated-LS','estimated-arx');
xlabel('time');
title(" real & estimated output")

%% Variance of noise & Covariance of estimate

e=y_val-y_hat;
p=2*na;
var_noise=(e'*e)/(N-p);
cov_est=(var_noise)*(inv(Ut*U));

disp(['var: ', num2str(var_noise)]);
disp(['det(Ut*U): ', num2str(DET)]);

