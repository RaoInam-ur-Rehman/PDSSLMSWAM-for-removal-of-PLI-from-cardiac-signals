clear all; close all; clc;
%% Data loading HRECG
fs=1000;T=1/fs;N=10000;
n=0:1:N-1;
load ecg;
data1=ecg_normal_simulator(:,3);
data=data1-mean(data1); data=data(1:N);
data=1.5*data+.074;%for peak to peak normalization to 1 of 3rd clmn of ecg_normal_simulator sig
SP=(1/N)*sum(data.^2)
data=data';
%% fixed freq PLI disturbance
SNR_in=3;
PLIfreq=48.79;
PLI1=sin(2*pi*PLIfreq*T*n);
PLI3=sin(2*pi*3*PLIfreq*T*n); 
PLI5=sin(2*pi*5*PLIfreq*T*n);
PLI7=sin(2*pi*7*PLIfreq*T*n);
PLI9=sin(2*pi*9*PLIfreq*T*n);
% PLIfreq=PLIfreq*ones(1,N);% for graph only in case of fixed freq
% PLI=(PLI1/10+PLI3/50+PLI5/150+PLI7/250+PLI9/450);
% PLI=(PLI1);
NP=SP/(10.^(SNR_in/10))
A1=(sqrt((.89)*2*NP));
A3=(sqrt(2*NP*(.05)));
A5=(sqrt(2*NP*(.03)));
A7=(sqrt(2*NP*(.02)));
A9=(sqrt(2*NP*(.01)));
PLI1=abs(A1).*PLI1;
PLI3=abs(A3).*PLI3;
PLI5=abs(A5).*PLI5;
PLI7=abs(A7).*PLI7;
PLI9=abs(A9).*PLI9;
PLI=(PLI1+PLI3+PLI5+PLI7+PLI9);
PLI=PLI;
NP3=(1/N)*sum(PLI.^2)
AAA=SP/NP3;
SNR=10*log10(AAA)
ecg_con=data+PLI;
ecg=ecg_con;
%%
freqest=PLIfreq*ones(1,N);

C=[1 0 1 0 1 0 1 0 1 0];
G=eye(10);
xpost=zeros(10,N);
K=zeros(10,N);
K(:,1)=[0;0;0;0;0;0;0;0;0;0];
xpri=zeros(10,N);
sci=zeros(10,N);
mu=zeros(1,N);
mu(1,1)=.01;
alpha=0.00001;
y=ecg;
for k=2:N
    fm=freqest(1,k);  
    w=2*pi*fm;
   A=[cos(w*T) sin(w*T) 0 0 0 0 0 0 0 0;  
    -sin(w*T) cos(w*T) 0 0 0 0 0 0 0 0;
    0 0 cos(3*w*T) sin(3*w*T) 0 0 0 0 0 0;
    0 0 -sin(3*w*T) cos(3*w*T) 0 0 0 0 0 0;
    0 0 0 0 cos(5*w*T) sin(5*w*T) 0 0 0 0;
    0 0 0 0 -sin(5*w*T) cos(5*w*T) 0 0 0 0;
    0 0 0 0 0 0 cos(7*w*T) sin(7*w*T) 0 0;
    0 0 0 0 0 0 -sin(7*w*T) cos(7*w*T) 0 0;
    0 0 0 0 0 0 0 0 cos(9*w*T) sin(9*w*T);
    0 0 0 0 0 0 0 0 -sin(9*w*T) cos(9*w*T)];
           if k==2
             ypred(k)=C*xpri(:,k);%predicted output
%              K=inv(phi)*C';    %5.6          
             E(k)=y(k)-ypred(k); 
             AB=(A - K(:,k)*C*A);
             sci(:,k)=AB*sci(:,k-1) + G*C'*E(k);
             mu(k)=mu(k-1)+alpha*sci(:,k-1)'*A'*C'*E(k);
            K(:,k)=mu(k)*G*C';   
%              E(k)=PLI(k)-ypred(k); %% for Hassan raza algo only
             xpost(:,k)=xpri(:,k)+K(:,k)*E(k);
             yest(k)=C*xpost(:,k);%%  estimated o/p
             EE(k)=y(k)-yest(k);
        else        
%              phi=L*inv(A')*phi*inv(A)+C'*C;
             xpri(:,k)=A*xpost(:,k-1);% Initial Pridicted States  
             ypred(k)=C*xpri(:,k);
%              K=inv(phi)*C';
            E(k)=y(k)-ypred(k);        
            AB=(A - K(:,k)*C*A);
         sci(:,k)=AB*sci(:,k-1) + G*C'*E(k);
             mu(k)=mu(k-1)+alpha*sci(:,k-1)'*A'*C'*E(k);
            K(:,k)=mu(k)*G*C'; 
             
%              E(k)=PLI(k)-ypred(k); %% for Hassan raza algo only
             xpost(:,k)=xpri(:,k)+K(:,k)*E(k);
%              
%                   Gamma=1;
% %               alpha=K'*C';
%               alpha=C*K;
%               xpost=xpri+2*Gamma*alpha(1)*C'*E(k);

             yest(k)=C*xpost(:,k);%%  changed
             EE(k)=y(k)-yest(k);
        end
    
end
%%
% figure(10);plot(n*T,PLI,n*T,yest,'Linewidth',1.5);xlabel('time(sec)');ylabel('Amplitude');
% title ('Tracking of SSLMS during initialization');xlim([0 1])
% legend('PLI','Estimated PLI')
figure (2);subplot(211);plot(n*T,ecg,'Linewidth',1.5);
xlabel('time(sec)'); ylabel('Normalized Amplitude'); title('Contaminated ECG')
ylim([-0.5 1.2])
subplot(212);plot(n*T,EE,'Linewidth',1.5);xlabel('time(sec)');ylabel('Normalized Magnitude')
title ('Filtered ECG'); %axis([0 N -0.01 0.01])
ylim([-0.5 1.2])
figure(3);plot(n*T,PLI-yest,'Linewidth',1.5);xlabel('time(sec)');ylabel('Amplitude')
title ('Tracking error with mu=0.05'); %axis([0 N -0.1 0.1])
%%
% figure
% K=N;
% ecgfft=abs(fft(ecg,N));%for FFT
% xaxis= (1:N).*fs/K;% defining x axis in Hz
% plot(xaxis(1:N/2),20*log(ecgfft(1:N/2)),'LineWidth',1.2)% for plotting +ive freq only
% title('Frequency Spectrum of Tracked Signal');
% xlabel('Frequency (Hz)');ylabel('Magnitude (dB)'); %Grid on
% axis([0 500 -80 150])
% hold on
% PLIfft=abs(fft(yest1,N));%for FFT
% K=N;
% xaxis= (1:N).*fs/K;% defining x axis in Hz
% plot(xaxis(1:N/2),20*log(PLIfft(1:N/2)),'r-.','LineWidth',1)% for plotting +ive freq only
% % title('Frequency Spectrum of Simulated HRECG');
% % xlabel('Frequency (Hz)');ylabel('Magnitude (dB)'); %Grid on
% axis([0 500 -80 150])
% legend('Contaminated HRECG','Tracked PLI Signal')
%%
% figure(1111);
% [ha, pos] = tight_subplot(4,1,[.01 .03],[.1 .05],[.1 .03]);
% axes(ha(1));plot(n*T,data);
% % set(gca,'xtick',[],'ytick',[])
% 
% ylabel('Normalized Magnitude')
% legend ('Clean HRA IEGM signal ');
% set(gca,'xtick',[])
% title('Filteration Comparison')
% axes(ha(2)); plot(n*T,ecg);
% 
% % ylabel('Normalized Magnitude')
% % xlabel('Samples');ylabel('Normalized Magnitude')  
% legend ('PLI corrupted HRA IEGM signal');
% set(gca,'xtick',[])
% axes(ha(3)); plot(n*T,EE,n*T,error,'r');
% % xlabel('Samples');ylabel('Normalized Magnitude')
% legend ('Sequentially operated SSLMSWAM','error','Orientation','horizontal'); 
% set(gca,'xtick',[])
% axes(ha(4)); plot(n*T,EE1,n*T,error1,'r');
% % xlabel('Samples');
% legend ('Proposed PD-SSLMSWAM','error','Orientation','horizontal');
% %set(gca,'xtick',[])
% axis([0 10 -.4 0.6]);
% xlabel('Times (Sec)');