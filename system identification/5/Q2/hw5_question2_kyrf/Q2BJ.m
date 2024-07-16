clc
clear
close all;

load('HW5_question2.mat')

Ts=0.5;
u_id=Z3.u(1:600);
u_val=Z3.u;
y_id=Z3.y(1:600);
y_val=Z3.y;
tt=0.5:0.5:300;
t=0.5:0.5:500;



%% Box_jenkins model

N=length(y_id);
K=1;

for i=1:10
    na=i;
    nb=i;
    nd=i;nf=nd;
    p=na+nb;
    model=bj([y_id u_id],[na nb nd nf 1]);
    aa=model.F;
    bb=model.B;
%     cc=model.C;
%     dd=model.D;

    G2 = tf(bb,aa, Ts);
    [y_hat, t_hat] = lsim(G2, u_id, tt);
    
%     G3 = tf(cc,dd, 0.1);
%     [e_hat, t_hat2] = lsim(G3, v_id, tt);
%     y_hat=y_hat+e_hat;
%     
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

%% Box_jenkins Function
na=4;nb=na;nc=na;nd=nc;p=na+nb;
model=bj([y_id u_id],[na nb nc nd 1]);
aa=model.F;
bb=model.B;
cc=model.C;
dd=model.D;

G2 = tf(bb,aa, Ts);
[y_hat2, t_hat] = lsim(G2, u_val, t);

%% validation of identification
SSE = sum((y_val(600:1000) - y_hat2(600:1000)).^2);
MSE = mean((y_val(600:1000) - y_hat2(600:1000)).^2);
mean_y_real = mean(y_val(600:1000));
SSR = sum((y_hat2(600:1000) - mean_y_real).^2);
R_squared = 1 - SSE / SSR;

disp(['SSE: ', num2str(SSE)]);
disp(['MSE: ', num2str(MSE)]);
disp(['R-squared: ', num2str(R_squared)]);


figure;
plot( t(600:1000),y_val(600:1000), 'b',t(600:1000), y_hat2(600:1000), 'r');
legend('real','estimated');
xlabel('time');
title("using BJ model")

%% varince of noise
e=y_val(600:1000)-y_hat2(600:1000);
S=e'*e;%sum sqr of errors
sigma_hat_2=(S)/(N-p);%noise varianse
figure;
plot(e)
title('error')
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
title('error correlation function')


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
plot( tt,Rex, 'b',tt, V_ex, 'r',tt,-V_ex,'r')
title('error & input cross correlation function')
%% DC gain
gain_dc=dcgain(G2);
disp(['DC_gain: ', num2str(gain_dc)]);

%% frequency peak
[gpeak,fpeak] = getPeakGain(G2);
disp(['gain peak: ', num2str(gpeak)]);
disp(['Frequncy of gain peak: ', num2str(fpeak)]);
