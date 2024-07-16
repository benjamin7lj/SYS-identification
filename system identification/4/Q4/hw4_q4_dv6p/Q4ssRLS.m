clc
clear
close all;
%%  readind data and seprating to identification and validation
load('HW4_Q4.mat')

y_id=detrend(Z2.y(1:500));
u_id=detrend(Z2.u(1:500));

y_id0=detrend(Z2.y(1:20));
u_id0=detrend(Z2.u(1:20));
tt=0:1:499;

y_val=detrend(Z2.y(501:1000));
u_val=detrend(Z2.u(501:1000));


%% zero stage setup

N=length(y_id);
na=4;nb=4;
ut=zeros(2*na,1);
Pt0=10000*eye(2*na+1,2*na+1);

%% zero stage theta0 with LS

for i=1:na
   for j=1:length(y_id0)
        if j>i
           U0(j,i)=-y_id0(j-i);  
        else
           U0(j,i)=0;
        end 
    end
end
for i=na+1:na+nb
  for j=1:length(u_id0)
        if j>(i-na)
           U0(j,i)=u_id0(j-i+na+1);  
        elseif j==i-na
             U0(j,i)=u_id0(1);
        else
           U0(j,i)=0;
        end 
  end
end

Ut0=U0';
tetha_hat_Ls=(inv(Ut0*U0))*(Ut0*y_id0);

tetha_hat_Ls0=zeros(na+nb+1,1);
tetha_hat_Ls0(1,1)=1;
for j=1:2*na
    tetha_hat_Ls0(j+1,1)=tetha_hat_Ls(j);
end

%% implimentation of RLS method for ARX model
for t=na+1:N
    for i=1:na
        ut(i,1)=-y_id(t-i);
        ut(i+na,1)=u_id(t-i+1);
    end
    ut(na+nb+1,1)=u_id(t-2);
    utt=ut';
    Kt=(Pt0*ut)/(1+utt*Pt0*ut);
    tetha_hat_RLS=tetha_hat_Ls0+Kt*(y_id(t)-utt*tetha_hat_Ls0);
    Pt0=(eye(2*na+1,2*na+1)-Kt*utt)*Pt0;
    
    if norm(tetha_hat_RLS-tetha_hat_Ls0,2)<=10^-8 && t>=12
        break
    end
    
    tetha_hat_Ls0=tetha_hat_RLS;
end

%% making transfer function
for i=1:na
 a(i)=tetha_hat_RLS(i);
end
for i=1:nb+1
 b(i)=tetha_hat_RLS(i+nb);
end

z=tf('z',1);

GZ=(b(1)*(z^4)+b(2)*(z^3)+b(3)*z^2+b(4)*z+b(5))/...
 (z^4+a(1)*z^3+a(2)*z^2+a(3)*z+a(4));

% GZ=(b(1)*(z^2)+b(2)*(z)+b(3))/...
%  (z^2+a(1)*z+a(2));


[y_hat, t_hat] = lsim(GZ, u_id, tt);

%% validation of identification
SSE = sum((y_id - y_hat).^2);
MSE = mean((y_id - y_hat).^2);
mean_y_real = mean(y_id);
SSR = sum((y_hat - mean_y_real).^2);
R_squared = 1 - SSE / SSR;

disp(['SSE: ', num2str(SSE)]);
disp(['MSE: ', num2str(MSE)]);
disp(['R-squared: ', num2str(R_squared)]);


figure;
plot( tt,y_id, 'b',tt, y_hat, 'r');
legend('real','estimated');
xlabel('time');
title("using arx model")

%% validation

[y_hat2, t_hat2] = lsim(GZ, u_val, tt);

%validation of validation data
SSE = sum((y_val - y_hat2).^2);
MSE = mean((y_val - y_hat2).^2);
mean_y_real = mean(y_val);
SSR = sum((y_hat2 - mean_y_real).^2);
R_squared = 1 - SSE / SSR;

disp(['SSE_val: ', num2str(SSE)]);
disp(['MSE_val: ', num2str(MSE)]);
disp(['R-squared_val: ', num2str(R_squared)]);


figure;
plot( tt,y_val, 'b',tt, y_hat2, 'r');
legend('real','estimated');
xlabel('time');
title("val_using arx model")

