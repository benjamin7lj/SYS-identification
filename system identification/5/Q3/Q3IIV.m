clc
clear
close all;

load('q3_401126125.mat')

u_id=u1;
y_id=y1;
%% IIV

N=length(y1);
na=4;nb=na;
nd=1;
p=na+nb;
U = zeros(N, na+nb);
for j = 1:na
    U(j+1:N, j) = -y_id(1:N-j, 1);
end
for j = 1:nb
    U(j+nd:N, na+j) = u_id(1:N-j+1-nd, 1);
end
Z=U;

Ut=U';
tetha_hat_Ls=(inv(Ut*U))*(Ut*y_id);

tetha_hat_1=tetha_hat_Ls;

for t=1:10
    y_hat_0=U*tetha_hat_1;
    for i = 1:na
        Z(i+1:N, i) = -y_hat_0(1:N-i, 1);
    end
    tetha_hat_IIV=(inv(Z'*U)*Z')*y_id;
    tetha_hat_1=tetha_hat_IIV;
end

y_hat_IIV=U*tetha_hat_IIV;
y_hat=y_hat_IIV;
% Gz = tf(tetha_hat_IIV(na+1:na+nb),[1, tetha_hat_IIV(1:na)], Ts);
% [y_hat_IIV, t_mod] = lsim(Gz, u_val, tt);
G=tf(tetha_hat_IIV(na+1:na+nb)',[1 tetha_hat_IIV(1:na)']);
figure;
pzmap(G)
title("pzmap of best G")

e=y_id-y_hat_IIV;
Le=length(e);
S=e'*e;%sum sqr of errors
sigma_hat_2=(S)/(N-p);%noise varianse
%% validation of identification
SSE = sum((y_id - y_hat_IIV).^2);
MSE = mean((y_id - y_hat_IIV).^2);
mean_y_real = mean(y_id);
SSR = sum((y_hat_IIV - mean_y_real).^2);
R_squared = 1 - SSE / SSR;

disp(['SSE_IIV: ', num2str(SSE)]);
disp(['MSE_IIV: ', num2str(MSE)]);
disp(['R-squared_IIV: ', num2str(R_squared)]);
disp(['variance of noise: ', num2str(sigma_hat_2)]);
figure;
plot(e)
title('error')
figure;
plot( 1:Le,y_id, 'b',1:Le, y_hat, 'r')
legend('real','estimated')
title('real and estimated output')



%% Error correlation function
% Le=length(e);
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
plot( 1:Le,Ree_sigma, 'b',1:Le, V_ee, 'r',1:Le,-V_ee,'r')
title('error  correlation function')


%% Error & input cross correlation function
Le=length(e);
for i=1:Le
   for j=1:Le
        if j>i
           XX(j,i)=u_id(j-i);  
        else
           XX(j,i)=u_id(Le+(j-i));
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
plot( 1:Le,Rex, 'b',1:Le, V_ex, 'r',1:Le,-V_ex,'r')
title('error & input cross correlation function')

