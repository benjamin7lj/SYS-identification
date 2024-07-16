clc
clear
close all;

load ('401126125.mat');

y_id=id.y;
u_id=id.u;

y_val=val.y;
u_val=val.u;
t = 0:0.1:24.9;


figure(1);
subplot(2,1,1)
plot(u_id);
title("input_id")
subplot(2,1,2)
plot(y_id);
title("output_id")
T=0.1;

figure(2);
subplot(2,1,1)
plot(u_val);
title("input_val")
subplot(2,1,2)
plot(y_val);
title("output_val")

na=5;%order of denum
nb=na-1;%order of num
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
           U(j,i)=u_id(j-i+na+1);  
        elseif j==i-na
             U(j,i)=u_id(1);
        else
           U(j,i)=0;
        end 
  end
end

Ut=U';
tetha_hat_LS=(inv(Ut*U))*(Ut*y_id);
y_hat=U*tetha_hat_LS;

%validation of identification
SSE = sum((y_id - y_hat).^2);
MSE = mean((y_id - y_hat).^2);
mean_y_real = mean(y_id);
SSR = sum((y_hat - mean_y_real).^2);
R_squared = 1 - SSE / SSR;


disp(['SSE: ', num2str(SSE)]);
disp(['MSE: ', num2str(MSE)]);
disp(['R-squared: ', num2str(R_squared)]);



for i=1:nb
    a(i)=tetha_hat_LS(i);
end
for i=1:na
    b(i)=tetha_hat_LS(i+nb);
end


z=tf('z',.1);

GZ=(b(1)*(z^4)+b(2)*(z^3)+b(3)*z^2+b(4)*z+b(5))/...
    (z^4+a(1)*z^3+a(2)*z^2+a(3)*z+a(4));

[num,den]=tfdata(GZ);
GZ_inv=filt(num,den,.1);

[y_model, t_model] = lsim(GZ, u_val, t);

figure;
plot( t,y_val, 'b',t_model, y_model, 'r');
legend('real','estimated');
xlabel('time');

%using arx function
model=arx([y_id u_id],[na nb 1]);
aa=model.A;
bb=model.B;
GZ2=(bb(1)*(z^4)+bb(2)*(z^3)+bb(3)*z^2+bb(4)*z+bb(5))/...
    (z^4+aa(2)*z^3+aa(3)*z^2+aa(4)*z+aa(5));
[y_model2, t_model2] = lsim(GZ2, u_val, t);

figure;
plot( t,y_val, 'b',t_model2, y_model2, 'r');
legend('real','estimated');
xlabel('time');
title("using arx function")

%validation of validation part
SSE_val = sum((y_val - y_model).^2);
MSE_val = mean((y_val - y_model).^2);
mean_y_real_val = mean(y_val);
SSR_val = sum((y_model - mean_y_real_val).^2);
R_squared_val = 1 - SSE_val / SSR_val;

disp(['SSE_val: ', num2str(SSE_val)]);
disp(['MSE_val: ', num2str(MSE_val)]);
disp(['R-squared_val: ', num2str(R_squared_val)]);

%validation of validation part with arx function
SSE_val_arx = sum((y_val - y_model2).^2);
MSE_val_arx = mean((y_val - y_model2).^2);
mean_y_real_val_arx = mean(y_val);
SSR_val_arx = sum((y_model2 - mean_y_real_val_arx).^2);
R_squared_val_arx = 1 - SSE_val_arx / SSR_val_arx;

disp(['SSE_val_arx: ', num2str(SSE_val_arx)]);
disp(['MSE_val_arx: ', num2str(MSE_val_arx)]);
disp(['R-squared_val_arx: ', num2str(R_squared_val_arx)]);