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

figure(1);

subplot(2,1,1)
plot(t,u_val);
title("u(t)");

subplot(2,1,2)
plot(t,y_val);
title("y(t)")
hold on

%% LEAST SQR with ARX model

N=length(y_id);
K=10;
Cov2=zeros(40,20);
for i=1:110
    na=i;
    nb=na;
    p=na+nb;
    nd=1;
    model=arx([y_id u_id],[na nb 1]);
    aa=model.A;
    bb=model.B;

    G1 = tf(bb,aa, Ts);
    [y_hat, t_hat1] = lsim(G1, u_id, tt);
    U = zeros(N, na+nb);
    for j = 1:na
        U(j+1:N, j) = -y_id(1:N-j, 1);
    end
    for j = 1:nb
        U(j+nd:N, na+j) = u_id(1:N-j+1-nd, 1);
    end
    Ut=U';
%     tetha_hat_LS=((inv(Ut*U))*(Ut*y_id))';
    %y_hat=U*tetha_hat_LS';
    %Gz = tf(tetha_hat_LS(na+1:p), tetha_hat_LS(1:na), Ts);
    %[y_hat, ttt]=lsim(Gz, u_id, tt);

    e=y_id-y_hat;
    S(i)=e'*e;%sum sqr of errors
    sigma_hat_2(i)=(S(i))/(N-p);%noise varianse
    AIC(i)=N*(log10(S(i)))+K*p;%AIC method
    cov_tetha_i=(sigma_hat_2(i)).*(inv(Ut*U));
    Cov=diag(cov_tetha_i);
    Cov2(1:length(Cov),i)=Cov;
    Cov3(i)=sum(Cov)/p;
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

figure;
contour(log(Cov2),'Fill','on');
xlabel('system order');
title("covariance matrix of tetha");

figure;
bar(Cov3);
xlabel('system order');
title("covariance matrix of tetha");

%% ARX Function

na=9;nb=9;
model=arx([y_id u_id],[na nb 1]);
aa=model.A;
bb=model.B;

G2 = tf(bb, aa, Ts);
%G3=(G2)/(1-1*G2);
%u_val2=Z3.u(591:1000);
%t=0.5:0.5:205;
[y_hat2, t_hat] = lsim(G2, u_val, t);
%y_hat2(1:12)=y_hat2(1:12)+1*y_id(589:600);
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
title("using ARX model")


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
title('error  correlation function')


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
