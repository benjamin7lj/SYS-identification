clc
clear
close all;

G=tf([5 0 -5],[2 14 36 44 24]);
P=pole(G);
Ts=(1/(max(abs(P))))*0.1;

%identification using WGN signal and LS
u = wgn(100,1,0);
t=0:Ts:Ts*99;

[y_id,t_id]=lsim(G,u,t);

figure(1);
subplot(2,1,1)
plot(t,u);
title("input")
subplot(2,1,2)
plot(t_id,y_id);
title("output")


na=2;%order of denum
nb=na;%order of num


N=length(y_id);

for i=1:na
   for j=1:N
        if j>i
           U(j,i)=-y_id(j-i);  
        else
           U(j,i)=0;
        end 
    end
end
for i=na+1:na+nb
  for j=1:N
        if j>(i-na)
           U(j,i)=u(j-i+na+1);  
        elseif j==i-na
             U(j,i)=u(1);
        else
           U(j,i)=0;
        end 
  end
end

Ut=U';


tetha_hat_LS=(inv(Ut*U))*(Ut*y_id);
y_hat=U*tetha_hat_LS;

%validation
SSE = sum((y_id - y_hat).^2);
MSE = mean((y_id - y_hat).^2);
mean_y_real = mean(y_id);
SSR = sum((y_hat - mean_y_real).^2);
R_squared = 1 - SSE / SSR;

figure;
plot(t, y_id, 'b', t, y_hat, 'r');
legend('real','input','estimated');
xlabel('time');
title("WGN");



disp(['SSE: ', num2str(SSE)]);
disp(['MSE: ', num2str(MSE)]);
disp(['R-squared: ', num2str(R_squared)]);

%identification using PRBS signal with LS
u_prbs=prbs(7,100);

[y_id_prbs,t_id_prbs]=lsim(G,u_prbs,t);

figure(3);
subplot(2,1,1)
plot(t,u_prbs);
title("PRBS_input")
subplot(2,1,2)
plot(t_id_prbs,y_id_prbs);
title("output")


for i=1:na
   for j=1:N
        if j>i
           U1(j,i)=-y_id_prbs(j-i);  
        else
           U1(j,i)=0;
        end 
    end
end
for i=na+1:na+nb
  for j=1:N
        if j>(i-na)
           U1(j,i)=u_prbs(j-i+na+1);  
        elseif j==i-na
             U1(j,i)=u(1);
        else
           U1(j,i)=0;
        end 
  end
end

U1t=U1';


tetha_hat_LS_prbs=(inv(U1t*U1))*(U1t*y_id_prbs);
y_hat_prbs=U1*tetha_hat_LS_prbs;

%validation
SSE_prbs = sum((y_id_prbs - y_hat_prbs).^2);
MSE_prbs = mean((y_id_prbs - y_hat_prbs).^2);
mean_y_real_prbs = mean(y_id_prbs);
SSR_prbs = sum((y_hat_prbs - mean_y_real_prbs).^2);
R_squared_prbs = 1 - SSE_prbs / SSR_prbs;

figure;
plot(t, y_id_prbs, 'b', t, y_hat_prbs, 'r');
legend('real','input','estimated');
xlabel('time');
title("PRBS")



disp(['SSE_prbs: ', num2str(SSE_prbs)]);
disp(['MSE_prbs: ', num2str(MSE_prbs)]);
disp(['R-squared_prbs: ', num2str(R_squared_prbs)]);
