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
tt=0:0.1:49.9;

%% OE model

N=length(y_id);
K=10;

for i=1:40
    na=i;
    nb=i;
    p=na+nb;
    model=oe([y_id u_id],[na nb 0]);
    aa=model.F;
    bb=model.B;

    G2 = tf(bb,aa, 0.1);
    [y_hat, t_hat] = lsim(G2, u_id, tt);
    

    
    e=y_id-y_hat;
    S(i)=e'*e;%sum sqr of errors
    sigma_hat_2(i)=(S(i))/(N-p);%noise varianse
    AIC(i)=N*(log10(S(i)))+K*p;%AIC method
    %cov_tetha_i=(sigma_hat_2(i)).*(inv(Ut*U));
end

figure;
bar(S);
xlabel('system order');
title("fitting principle");

figure;
bar(sigma_hat_2);
xlabel('system order');
title("variance of noise");

figure;
bar(AIC);
xlabel('system order');
title("AIC criterion");


%% OE Function
t=0:0.1:50;
na=3;nb=na;p=na+nb;
model=oe([y_id u_id],[na nb 0]);
aa=model.F;
bb=model.B;
% cc=model.C;

G2 = tf(bb,aa, 0.1);
[y_hat2, t_hat] = lsim(G2, u_val, t);
% G3 = tf(cc,aa, 0.1);
% [e_hat2, t_hat2] = lsim(G3, v_val, t);

y_hat=y_hat2 + v_val;
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
plot( t,y_val, 'b',t, y_hat, 'r');
legend('real','estimated');
xlabel('time');
title("using OE model")


%% varince of noise
e=y_val-y_hat;
S=e'*e;%sum sqr of errors
sigma_hat_2=(S)/(N-p);%noise varianse

%% Error correlation function
Le=length(e);
for i=1:Le
   for j=1:Le
        if j>i
           EE(j,i)=e(j-i);  
        else
           EE(j,i)=e(Le+(j-i));
        end 
    end
end

for k=1:Le
    Ree(k)=(1/Le)*((EE(1:Le,k)')*EE(1:Le,1));
    V_ee(k)=2*((1/Le)^0.5);
end

Ree_sigma=(1/sigma_hat_2).*Ree;
tt=1:1:Le;
figure;
plot( tt,Ree_sigma, 'b',tt, V_ee, 'r',tt,-V_ee,'r')
title('error & input cross correlation function')


%% Error & input cross correlation function
Le=length(e);
for i=1:Le
   for j=1:Le
        if j>i
           XX(j,i)=u_val(j-i);  
        else
           XX(j,i)=u_val(Le+(j-i));
        end 
    end
end

for k=1:Le
    Rxx(k)=(1/Le)*((XX(1:Le,k)')*XX(1:Le,1));
end

PP=Ree*Rxx';
for k=1:Le
    Rex(k)=(1/Le)*((EE(1:Le,k)')*XX(1:Le,1));
    V_ex(k)=2*((PP/Le)^0.5);
end


figure;
plot( t,Rex, 'b',t, V_ex, 'r',t,-V_ex,'r')
title('error & input cross correlation function')
