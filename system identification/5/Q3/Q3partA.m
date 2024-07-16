clc
clear
close all;

load('q3_401126125.mat')
u_id=u1;
y_id=y1;

%% LEAST SQR with ARX model

N=length(y_id);
K=10;
d_prime=1;

for i=1:80
    na=i;
    nb=na;
    p=na+nb;
    nd=0;
%     model=arx([y_id u_id],[na nb 1]);
%     aa=model.A;
%     bb=model.B;
% 
%     G1 = tf(bb,aa, Ts);
%     [y_hat, t_hat1] = lsim(G1, u_id, tt);
    U = zeros(N, na+nb);
    for j = 1:na
        U(j+1:N, j) = -y_id(1:N-j, 1);
    end
    for j = 1:nb
    U(j+nd:N, na+j) = u_id(1:N-j+1-nd, 1);
    end
    Ut=U';
    tetha_hat_LS=((inv(Ut*U))*(Ut*y2))';
    y_hat=U*tetha_hat_LS';

    e=y_id-y_hat;
    S(i)=e'*e;%sum sqr of errors
    sigma_hat_2(i)=(S(i))/(N-p);%noise varianse
    AIC(i)=N*(log10(S(i)))+K*p;%AIC method
    %cov_tetha_i=(sigma_hat_2(i)).*(inv(Ut*U));
    cov_tetha_i=(sigma_hat_2(i)).*(inv(Ut*U));
    Cov=diag(cov_tetha_i);
    Cov2(1:length(Cov),i)=Cov;
    Cov3(i)=sum(Cov)/p;
    
    Ry(i)=(1/N)*sum((U(1:4096,i)')*U(1:4096,1));%output correation function
    FPE(i)=0.5*(1+(2*d_prime*p)/N)*(1/N)*S(i);%FPE method
    
    Rxy(i)=-(1/N)*sum(U(1:4096,1)'*(U(1:4096,na+i)));%output & input cross correation function
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

figure;
bar(Ry);
xlabel('system order');
title("output correlation Function");

figure;
bar(FPE);
xlabel('system order');
title("FPE method");

figure;
bar(Rxy);
xlabel('system order');
title("output & input cross correlation Function");


[R,tau]=Cross_Correlation(u1,y1,10);

function [R,tau] = Cross_Correlation(x,y,tau_max)
N=length(y);
tau_max=10;
for tou=0:tau_max
    Ss=zeros(1,N);
for i=1:N
    if i+tou>N Ss(i)=x(i)*y(tou);
    else Ss(i)=x(i)*y(i+tou);
    end
end
R(tou+1)=(1/N)*sum(Ss);
end
tau=0:1:tau_max;
figure()
stem(tau,R)
end
