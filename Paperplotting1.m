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
%% HRECG test signal
figure(3)
plot(n*T,data,'LineWidth',2)% to plot graph in sec
title('Pure HRECG Signal'); xlabel('Time (seconds)');ylabel('Normalized Amplitude');%Grid on
axis([0 2 -0.2 1])
datafft=abs(fft(data,N));%for FFT
figure(4)
K=N;
xaxis= (1:N).*fs/K;% defining x axis in Hz
plot(xaxis(1:N/2),20*log(datafft(1:N/2)),'LineWidth',1.2)% for plotting +ive freq only
title('Frequency Spectrum of Pure HRECG');
xlabel('Frequency (Hz)');ylabel('Magnitude (dB)'); %Grid on
axis([0 500 -80 150])
%
%% 
figure(5)
plot(n*T,ecg,'LineWidth',1.2)% to plot graph in sec
title('Corrupted HRECG Signal'); xlabel('Time (seconds)');ylabel('Normalized Amplitude');%Grid on
axis([0 2 -0.3 1.1])
ecgfft=abs(fft(ecg,N));%for FFT
figure(6)
K=N;
xaxis= (1:N).*fs/K;% defining x axis in Hz
plot(xaxis(1:N/2),20*log(ecgfft(1:N/2)),'LineWidth',1.1)% for plotting +ive freq only
title('Frequency Spectrum of Corrupted HRECG');
xlabel('Frequency (Hz)');ylabel('Magnitude (dB)'); %Grid on
axis([0 500 -80 150])
%%
figure(7)
% plot(n,PLI1,'') 
% hold on;
% plot(n,PLI1,n,PLI3,n,PLI5,n,PLI7,n,PLI9) 
% hold on;
% plot(n,PLI1,n,PLI3,n,PLI5,n,PLI7,n,PLI9) 
% hold on;
plot(n,PLI1,'d-',n,PLI3,'^-',n,PLI5,'v-',n,PLI7,'p-',n,PLI9,'s-') 
hold on;
plot(n,PLI,'r','LineWidth',3)
title('PLI Harmonic Components and Composite Signal'); xlabel('Time (miliseconds)');ylabel('Normalized Amplitude');%Grid on
legend('1st Harmonic','3rd Harmonic','5th Harmonic','7th Harmonic','9th Harmonic','Composite PLI','Orientation','vertical')
% axis([0 45 -0.1 0.2])
xlim([0 55])
ylim([-0.2 0.35])
hold off