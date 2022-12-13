clear all; close all; clc;
%% Data loading HRECG
% fs=1000;T=1/fs;N=10000;
% n=0:1:N-1;
% load ecg;
% data1=ecg_normal_simulator(:,3);
% data=data1-mean(data1); data=data(1:N);
% data=1.5*data+.074;%for peak to peak normalization to 1 of 3rd clmn of ecg_normal_simulator sig
% SP=(1/N)*sum(data.^2)
% data=data';
%% Data loading UHF-ECG
% N=60000; fs=5000;
% T=1/fs; n=0:N-1;
% load ECG3
% data1=ecg(:,3);
% data=data1-mean(data1);
% data=data./max(data);
% data=data(1:N);
% data=movmean(data,117);
% data=data';
% SP=(1/N)*sum(data.^2)
%% Data HRA IEGM
N=20000; fs=2000; 
T=1/fs; n=0:N-1;
load parveen
data1=data(1:N,4)';
data=data1-mean(data1);
data=data./max(data);
data=data(1:N)';
data=movmean(data,15);
data=data';
SP=(1/N)*sum(data.^2)

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
%% Frequency estimation 
% WinLen=fs/4;
% Frange=[45 55];
% [Fest1,FREQMAX1,Iter1, freqest1] = SAIDFT(ecg, fs, WinLen, Frange);
%%
%  freqest=freqest1(1:N);
freqest=PLIfreq*ones(1,N);

% mu1=0.01; %for HRECG HRA-IEGM signal aplha=0.00001
mu1=0.01; %for UHF-ECG signal aplha=0.000001
% mu2=0.09;
SysOrd=10;
[yest,EE] =SSLMSwithvariant(ecg,fs,freqest,mu1,SysOrd,'SSLMSWAM');
error=data-EE;
[yest1,EE1] =SSLMSwithvariant_PD(ecg,fs,freqest,mu1,SysOrd,'SSLMSWAM');
error1=data-EE1;
%%
figure (2);subplot(211);plot(n*T,ecg,'Linewidth',1.5);
xlabel('time(sec)'); ylabel('Normalized Amplitude'); title('Contaminated ECG')
ylim([-0.5 1.2])
subplot(212);plot(n*T,EE1,'Linewidth',1.5);xlabel('time(sec)');ylabel('Normalized Magnitude')
title ('Filtered ECG'); %axis([0 N -0.01 0.01])
ylim([-0.5 1.2])
figure(3);plot(n*T,PLI-yest,'Linewidth',1.5);xlabel('time(sec)');ylabel('Amplitude')
title ('Tracking error with mu=0.05'); %axis([0 N -0.1 0.1])
%%
figure(1111);
[ha, pos] = tight_subplot(4,1,[.01 .03],[.1 .05],[.1 .03]);
axes(ha(1));plot(n*T,data,'Linewidth',1.5);
% set(gca,'xtick',[],'ytick',[])
axis([0 10 -.2 .8]);

ylabel('Normalized Magnitude')
legend ('Clean HRA-IEGM signal ');
set(gca,'xtick',[])
title('Filteration Comparison for HRA-IEGM Signal')
axes(ha(2)); plot(n*T,ecg,'Linewidth',1.5);

% ylabel('Normalized Magnitude')
% xlabel('Samples');ylabel('Normalized Magnitude')  
legend ('PLI corrupted HRA-IEGM signal');
set(gca,'xtick',[])
axis([0 10 -.3 .8]);

axes(ha(3)); plot(n*T,EE,n*T,error,'r','Linewidth',1.5);
% xlabel('Samples');ylabel('Normalized Magnitude')
legend ('Sequentially operated SSLMSWAM','error','Orientation','horizontal'); 
set(gca,'xtick',[])
axis([0 10 -.2 .8]);

axes(ha(4)); plot(n*T,EE1,n*T,error1,'r','Linewidth',1.5);
% xlabel('Samples');
legend ('Proposed PD-SSLMSWAM','error','Orientation','horizontal');
%set(gca,'xtick',[])
axis([0 10 -.2 .8]);
xlabel('Times (Sec)');
%%
figure
K=N;
ecgfft=abs(fft(ecg,N));%for FFT
xaxis= (1:N).*fs/K;% defining x axis in Hz
plot(xaxis(1:N/2),20*log(ecgfft(1:N/2)),'LineWidth',1.2)% for plotting +ive freq only
title('Frequency Spectrum of Tracked Signal');
xlabel('Frequency (Hz)');ylabel('Magnitude (dB)'); %Grid on
axis([0 500 -80 150])
hold on
PLIfft=abs(fft(yest1,N));%for FFT
K=N;
xaxis= (1:N).*fs/K;% defining x axis in Hz
plot(xaxis(1:N/2),20*log(PLIfft(1:N/2)),'r-.','LineWidth',1)% for plotting +ive freq only
% title('Frequency Spectrum of Simulated HRECG');
% xlabel('Frequency (Hz)');ylabel('Magnitude (dB)'); %Grid on
axis([0 500 -100 150])
legend('Contaminated HRECG','Tracked PLI Signal')