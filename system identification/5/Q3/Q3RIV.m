clc
clear
close all;

load('q3_401126125.mat')

u_id=u1;
y_id=y1;


%% zero stage setup

N=length(y1);
na=4;nb=4;nd=1;p=na+nb;
ut=zeros(2*na,1);
zt=zeros(2*na,1);
Pt0=10000*eye(2*na,2*na);
y_id0=y_id(1:20);
%% zero stage theta0 with LS

U0 = zeros(20, na+nb);
for j = 1:na
    U0(j+1:20, j) = -y_id(1:20-j, 1);
end
for j = 1:nb
    U0(j+nd:20, na+j) = u_id(1:20-j+1-nd, 1);
end
Ut0=U0';
tetha_hat_Ls0=(inv(Ut0*U0))*(Ut0*y_id0);
y_hat_RIV=U0*tetha_hat_Ls0;
% tetha_hat_Ls0=zeros(na+nb+1,1);
% tetha_hat_Ls0(1,1)=1;
% for j=1:2*na
%     tetha_hat_Ls0(j+1,1)=tetha_hat_Ls(j);
% end
mean_y_real = mean(y_id);
%% implimentation of RIV method for ARX model
y_hat_RIV(1:na+1)=y_id(1:na+1);
for j=1:10
    for t=na+1:N
        for i=1:na
           ut(i,1)=-y_id(t-i);
          ut(i+na,1)=u_id(t-i+1-nd);
        end
        for i=1:na
          zt(i,1)=-y_hat_RIV(t-i);
          zt(i+na,1)=u_id(t-i+1-nd);
        end
    %ut(na+nb+1,1)=u2(t-2-nd);
        utt=ut';
        ztt=zt';

        Pt0=Pt0-((Pt0*zt*utt*Pt0)/(1+utt*Pt0*zt));
        tetha_hat_RIV=tetha_hat_Ls0+Pt0*zt*(y_id(t)-utt*tetha_hat_Ls0);
    
%     utt=ut';
%     ztt=zt';
%     Kt=(Pt0*zt)/(1+utt*Pt0*zt);
%     tetha_hat_RLS=tetha_hat_Ls0+Kt*(y_id(t)-utt*tetha_hat_Ls0);
%     Pt0=(eye(2*na,2*na)-Kt*utt)*Pt0;
    
%     if norm(tetha_hat_RLS-tetha_hat_Ls0,2)<=10^-8 && t>=12
%         break
%     end

        y_hat_t=utt*tetha_hat_RIV;
        y_hat_RIV(t)=y_hat_t;
        e(t)=y_id(t)-y_hat_t;
        ssr(t)=y_hat_t-mean_y_real;
        tetha_hat_Ls0=tetha_hat_RIV;
     end
   
end

G=tf(tetha_hat_RIV(na+1:na+nb)',[1 tetha_hat_RIV(1:na)']);
figure;
pzmap(G)
title("pzmap of best G")


y_hat=y_hat_RIV;

Le=length(e);
S=e*e';%sum sqr of errors
sigma_hat_2=(S)/(N-p);%noise varianse
%% validation of identification

SSR = ssr*ssr';
R_squared = 1 - S / SSR;
disp(['SSE_RLS: ', num2str(S)]);
disp(['R-squared_RLS: ', num2str(R_squared)]);

disp(['variance of noise: ', num2str(sigma_hat_2)]);
figure;
plot(e)
title('error')
figure;
plot( 1:Le,y_id, 'b',1:Le, y_hat_RIV, 'r')
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
