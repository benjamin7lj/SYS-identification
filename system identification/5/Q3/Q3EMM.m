clc
clear
close all;

load('q3_401126125.mat')

u_id=u1;
y_id=y1;

%% zero stage setup 
N=length(y_id);
nv=1;ne=1;

v_hat=wgn(ne,1,0);
e_hat=wgn(ne,1,0);

na=4;nb=4;nd=1;p=na+nb;
ut=zeros(2*(na+ne),1);
Pt0=10000*eye(2*(na+ne),2*(na+ne));
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
tetha_hat_Ls=(inv(Ut0*U0))*(Ut0*y_id0);
tetha_hat_Ls0=[tetha_hat_Ls',0,0]';

mean_y_real = mean(y_id);
%% implimentation of EMM method for box-jenkins model
for t=na+1:N-1
    for i=1:na
        ut(i,1)=-y_id(t-i);
        ut(i+na,1)=u_id(t-i+1-nd);
    end
    %ut(2*na+1,1)=u_id(t-na);
    
    for i=1:ne
        ut(2*na+i,1)=-e_hat(ne+1-i);
        ut(2*na+ne+i,1)=v_hat(ne+1-i);
    end
%     for i=1:ne
%         ut2(i,1)=-e_hat(ne+1-i);
%         ut2(ne+i,1)=v_hat(ne+1-i);
%     end

    utt=ut';
    %utt2=ut2';
    Kt=(Pt0*ut)/(1+utt*Pt0*ut);
    tetha_hat_EMM=tetha_hat_Ls0+Kt*(y_id(t)-utt*tetha_hat_Ls0);
    Pt0=(eye(2*(na+ne),2*(ne+na))-Kt*utt)*Pt0;
    v_hat(ne)=y_id(t)-utt*tetha_hat_EMM;
    e_hat(ne)=y_id(t)-utt(1:2*na)*tetha_hat_EMM(1:2*na);
    for i=1:ne
        ut2(i,1)=-e_hat(ne+1-i);
        ut2(ne+i,1)=v_hat(ne+1-i);
    end
    utt2=ut2';
    if norm(tetha_hat_EMM-tetha_hat_Ls0,2)<=10^-8 && t>=12
        break
    end
    y_hat_t(t)=utt*tetha_hat_EMM-utt2*tetha_hat_EMM(2*na+1:2*(na+ne));
    e(t+1)=e_hat(ne);
    %e(t)=e_hat(ne-1);
    v(t+1)=v_hat(ne);
    %v(t)=v_hat(ne-1);
    ssr(t)=y_hat_t(t)-mean_y_real;
    tetha_hat_Ls0=tetha_hat_EMM;
end

%% PZplot
G=tf(tetha_hat_EMM(na+1:na+nb)',[1 tetha_hat_EMM(1:na)']);
figure;
pzmap(G)
title("pzmap of best G")
% U = zeros(N, na+nb+ne+nv);
% % sys=tf(tetha_hat_EMM(na+1:2*na)',tetha_hat_EMM(1:na)');
% % h = pzplot(sys);
% % grid on
% % P = pole(sys);
% for j = 1:na
%     U(j+1:N, j) = -y_id(1:N-j, 1);
% end
% for j = 1:nb
%     U(j+nd:N, na+j) = u_id(1:N-j+1-nd, 1);
% end
% for j=1:ne
%     U(j+1:N,na+nb+j)=-e(1:N-j);
% end
% for j=1:nv
%     U(j+nd:N,na+nb+ne+j)=v(1:N-j);
% end
%y_hat_EMM=U*tetha_hat_EMM;
y_hat_EMM=y_hat_t;
y_hat=y_hat_EMM';


S=e*e';%sum sqr of errors
Le=length(e);
sigma_hat_2=(S)/(N-p);%noise varianse
%% validation of identification
SSE = e*e';
SSR = ssr*ssr';
R_squared = 1 - SSE / SSR;
disp(['SSE_RELS: ', num2str(SSE)]);
disp(['R-squared_RELS: ', num2str(R_squared)]);
disp(['variance of noise: ', num2str(sigma_hat_2)]);
figure;
plot(e)
title('error')
figure;
plot( 1:Le-1,y_id(1:4095), 'b',1:Le-1, y_hat, 'r')
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


% %% Error correlation function
% Le=length(e);
% for i=1:Le
%    for j=1:Le
%         if j>i
%            EE(j,i)=e(j-i);  
%         else
%            EE(j,i)=e(Le+(j-i));
%         end 
%     end
% end
% 
% for k=1:Le
%     Ree(k)=(1/Le)*((EE(1:Le,k)')*EE(1:Le,1));
%     V_ee(k)=2*((1/Le)^0.5);
% end
% 
% Ree_sigma=(1/sigma_hat_2).*Ree;
% tt=1:1:Le;
% figure;
% plot( tt,Ree_sigma, 'b',tt, V_ee, 'r')
% title('error  correlation function')
% 
% %% Error & input cross correlation function
% Le=length(e);
% for i=1:Le
%    for j=1:Le
%         if j>i
%            XX(j,i)=u_id(j-i);  
%         else
%            XX(j,i)=u_id(Le+(j-i));
%         end 
%     end
% end
% 
% for k=1:Le
%     Rxx(k)=(1/Le)*((XX(1:Le,k)')*XX(1:Le,1));
% end
% 
% PP=Ree*Rxx';
% for k=1:Le
%     Rex(k)=(1/Le)*((EE(1:Le,k)')*XX(1:Le,1));
%     V_ex(k)=2*((PP/Le)^0.5);
% end
% 
% 
% figure;
% plot( tt,Rex, 'b',tt, V_ex, 'r');
% title('error & input cross correlation function')
