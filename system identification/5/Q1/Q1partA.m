clc
clear
close all;

load('q1_401126125.mat')

u_id=u(1:500);
u_val=u(500:1000);
y_id=y(1:500);
y_val=y(500:1000);
tt=0:0.1:50;
t=0:0.1:49.9;

%% LEAST SQR with ARX model

N=length(y_id);
K=3;
Cov2=zeros(200,100);
for i=1:120
    
    na=i;
    nb=i;
    p=na+nb;
    model=arx([y_id u_id],[na nb 0]);
    aa=model.A;
    bb=model.B;

    G1 = tf(bb,aa, 0.1);
    [y_hat, t_hat1] = lsim(G1, u_id, t);
    U = zeros(N, na+nb);
    for j = 1:na
        U(j+1:N, j) = -y_id(1:N-j, 1);
    end
    for j = 1:nb
    U(j:N, na+j) = u_id(1:N-j+1, 1);
    end
    Ut=U';
%     tetha_hat_LS=((inv(Ut*U))*(Ut*y_id))';
%     y_hat=U*tetha_hat_LS';

    e=y_id-y_hat;
    S(i)=e'*e;%sum sqr of errors
    sigma_hat_2(i)=(S(i))/(N-p);%noise varianse
    AIC(i)=N*(log10(S(i)))+K*p;%AIC method
    cov_tetha_i=(sigma_hat_2(i)).*(inv(Ut*U));
    Cov=diag(cov_tetha_i);
    Cov2(1:length(Cov),i)=Cov;
    Cov3(i)=trace(cov_tetha_i)/p;
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

% %% RLS methot for arx model
% 
% %% zero stage setup
% na=30;nb=30;nd=0;p=na+nb;
% ut=zeros(2*na,1);
% Pt0=10000*eye(2*na,2*na);
% y_id0=y_id(1:20);
% % zero stage theta0 with LS
% tetha_hat_Ls0=zeros(60,1);
% %% implimentation of RLS method for ARX model
% for t=na+1:N
%     for i=1:na
%         ut(i,1)=-y_id(t-i);
%         ut(i+na,1)=u_id(t-i+1-nd);
%     end
%     %ut(na+nb+1,1)=u2(t-2-nd);
%     utt=ut';
%     Kt=(Pt0*ut)/(1+utt*Pt0*ut);
%     tetha_hat_RLS=tetha_hat_Ls0+Kt*(y_id(t)-utt*tetha_hat_Ls0);
%     Pt0=(eye(2*na,2*na)-Kt*utt)*Pt0;
%     
%     if norm(tetha_hat_RLS-tetha_hat_Ls0,2)<=10^-8 && t>=12
%         break
%     end
%     
%     tetha_hat_Ls0=tetha_hat_RLS;
% end
% G=tf(tetha_hat_RLS(na+1:2*na)',[1,tetha_hat_RLS(1:na)'],0.1);
% 
% [y_hat, t_hat2] = lsim(G, u_val, tt);
%% ARX Function

na=25;nb=25;
model=arx([y_id u_id],[na nb 1]);
aa=model.A;
bb=model.B;

G2 = tf(bb,aa, 0.1);


[y_hat2, t_hat] = lsim(G2, u_val, tt);

%% validation of identification
SSE = sum((y_val - y_hat2).^2);
MSE = mean((y_val - y_hat2).^2);
mean_y_real = mean(y_val);
SSR = sum((y_hat2 - mean_y_real).^2);
R_squared = 1 - SSE / SSR;

disp(['SSE: ', num2str(SSE)]);
disp(['MSE: ', num2str(MSE)]);
disp(['R-squared: ', num2str(R_squared)]);


figure;
plot( tt,y_val, 'b',tt,y_hat2,'r');
legend('real','estimated','ARX');
xlabel('time');
title("using ARX model")


%% varince of noise
e=y_val-y_hat2;
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
