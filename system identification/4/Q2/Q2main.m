clc
clear
close all;
%% reading data
load('401126125-q2.mat')

y_id=id.y;
u_id=id.u;
Ts=id.Ts;

NY=round(length(y_id)*0.1);
y_id2=id.y(1:NY);
u_id2=id.u(1:NY);

y_id0=y_id(1:20);
u_id0=u_id(1:20);
tt=0:Ts:514*Ts;

y_val=val.y;
u_val=val.u;


N=length(y_id);
na=3*n;nb=na;
ut=zeros(2*na,1);
Pt0=100*eye(2*na+1,2*na+1);

%% estimation of tetha of stage zero witha Least sqr method
for i=1:na
   for j=1:length(y_id0)
        if j>i
           U0(j,i)=-y_id0(j-i);  
        else
           U0(j,i)=0;
        end 
    end
end
for i=n+1:na+nb
  for j=1:length(u_id0)
        if j>(i-n)
           U0(j,i)=u_id0(j-i+n+1);  
        elseif j==i-n
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
%% using RLS function to estimate transfer function parameters

[tetha_hat_RLS_total] = RLS_arx(tetha_hat_Ls0,Pt0, N , u_id , y_id);

tetha_hat_RLS_final=tetha_hat_RLS_total(end,:);

%% transfer function
for i=1:na
 a(i)=tetha_hat_RLS_final(i);
end
for i=1:nb+1
 b(i)=tetha_hat_RLS_final(i+nb);
end

z=tf('z',Ts);
GZ=(b(1)*(z^9)+b(2)*(z^8)+b(3)*(z^7)+b(4)*(z^6)+b(5)*(z^5)+...
    b(6)*z^4+b(7)*z^3+b(8)*z^2+b(9)*z+b(10))/...
    (z^9+a(1)*z^8+a(2)*z^7+a(3)*z^6+a(4)*z^5+...
    a(5)*z^4+a(6)*z^3+a(7)*z^2+a(8)*z+a(9));


[y_hat, t_hat] = lsim(GZ, u_val, tt);

%% validation of identification
SSE = sum((y_val - y_hat).^2);
MSE = mean((y_val - y_hat).^2);
mean_y_real = mean(y_val);
SSR = sum((y_hat - mean_y_real).^2);
R_squared = 1 - SSE / SSR;

disp(['SSE: ', num2str(SSE)]);
disp(['MSE: ', num2str(MSE)]);
disp(['R-squared: ', num2str(R_squared)]);


figure;
plot( tt,y_val, 'b',tt, y_hat, 'r');
legend('real','estimated');
xlabel('time');
title("using arx model")



%% matlab function for RLS
function [tetha_hat_RLS_total] = RLS_arx(tetha_hat_Ls0, Pt0 , N , u_id , y_id)
    n=3;
    na=3*n;nb=na;
    
    
    for t=na+1:N
        for i=1:na
            ut(i,1)=-y_id(t-i);
            ut(i+na,1)=u_id(t-i+1);
        end
        ut(na+nb+1,1)=u_id(t-nb);
        utt=ut';
        Kt=(Pt0*ut)/(1+utt*Pt0*ut);
        tetha_hat_RLS=tetha_hat_Ls0+Kt*(y_id(t)-utt*tetha_hat_Ls0);
        Pt0=(eye(2*na+1,2*na+1)-Kt*utt)*Pt0;
    
        if norm(tetha_hat_RLS-tetha_hat_Ls0,2)<=10^-12 && t>=12
            break
        end
    
        tetha_hat_Ls0=tetha_hat_RLS;
        tetha_hat_RLS_total(t,:)=tetha_hat_RLS';
    end
    
end
