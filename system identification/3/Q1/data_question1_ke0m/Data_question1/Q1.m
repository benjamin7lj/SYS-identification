clc
clear
close all;

load ('y.mat'); load ('u.mat'); load ('t.mat');

figure(1);
subplot(2,1,1)
plot(t,u);
title("input")
subplot(2,1,2)
plot(t,y);
title("output")
T=0.1;

na=4;%order of denum
nb=na-2;%order of num
N=length(y);

for i=1:na
   for j=1:N
        if j>i
           U(j,i)=-y(j-i);  
        else
           U(j,i)=0;
        end 
    end
end

for i=na+1:na+nb
  for j=1:N
        if j>(i-na+3) % for Z^-3 in transfer function there is 3 shift
           U(j,i)=u(j-i+na-2);  
        elseif j==i-na+3
             U(j,i)=u(1);
        else
           U(j,i)=0;
        end 
  end
end

Ut=U';

tetha_hat_LS=(inv(Ut*U))*(Ut*y);
y_hat=U*tetha_hat_LS;

%validation
SSE = sum((y - y_hat).^2);
MSE = mean((y - y_hat).^2);
mean_y_real = mean(y);
SSR = sum((y_hat - mean_y_real).^2);
R_squared = 1 - SSE / SSR;


figure;
plot(t, y, 'b', t, y_hat, 'r');
legend('real','input','estimated');
xlabel('time');

disp(['SSE: ', num2str(SSE)]);
disp(['MSE: ', num2str(MSE)]);
disp(['R-squared: ', num2str(R_squared)]);