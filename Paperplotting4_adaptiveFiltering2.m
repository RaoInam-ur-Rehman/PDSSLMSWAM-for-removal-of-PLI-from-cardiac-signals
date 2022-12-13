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
% PLI=(PLI1+PLI3+PLI5);
PLI=PLI;
NP3=(1/N)*sum(PLI.^2)
AAA=SP/NP3;
SNR=10*log10(AAA)
ecg_con=data+PLI;
ecg=ecg_con;
%%
freqest=PLIfreq*ones(1,N);
% mu1=0.01; %for HRECG signal aplha=0.00001
mu1=0.01; %for UHF-ECG signal aplha=0.000001
% mu2=0.09;
SysOrd=10;
[yest,EE] =SSLMSwithvariant(ecg,fs,freqest,mu1,SysOrd,'SSLMSWAM');
error=data-EE;
[yest1,EE1] =SSLMSwithvariant_PD(ecg,fs,freqest,mu1,SysOrd,'SSLMSWAM');
error1=data-EE1;
yfilter=ecg-yest;
%% adaptive notch filter
data=data';
ecg=ecg';
q=100;% quality factor
dftwindow=fs/4;
K=N;% plotting FFT
datafft=fft(data,N);
mdata=20*log(abs(datafft));% magnitude response of FFT
xaxis= (1:N).*fs/K;% defining x axis in Hz

y1=zeros(N,1);y3=zeros(N,1);y5=zeros(N,1);y7=zeros(N,1);y9=zeros(N,1);
for j=3:N
    if j<dftwindow
        fo=50;
    else
        fo=freqest(j);
    end
omega1=2*pi*fo/fs;C1=tan(0.5*omega1/q);beta1=cos(omega1);%for 1st harmonic
omega3=2*pi*3*fo/fs;C3=tan(0.5*omega3/q);beta3=cos(omega3);%for 3rd harmonic
omega5=2*pi*5*fo/fs;C5=tan(0.5*omega5/q);beta5=cos(omega5);%for 5th harmonic
omega7=2*pi*7*fo/fs;C7=tan(0.5*omega7/q);beta7=cos(omega7);%for 7th harmonic
omega9=2*pi*9*fo/fs;C9=tan(0.5*omega9/q);beta9=cos(omega9);%for 9th harmonic

y1(j)=1/(1+C1)*(2*beta1*y1(j-1)-(1-C1)*y1(j-2)+ecg(j)-2*beta1*ecg(j-1)+ecg(j-2));
y3(j)=1/(1+C3)*(2*beta3*y3(j-1)-(1-C3)*y3(j-2)+y1(j)-2*beta3*y1(j-1)+y1(j-2));
y5(j)=1/(1+C5)*(2*beta5*y5(j-1)-(1-C5)*y5(j-2)+y3(j)-2*beta5*y3(j-1)+y3(j-2));
y7(j)=1/(1+C7)*(2*beta7*y7(j-1)-(1-C7)*y7(j-2)+y5(j)-2*beta7*y5(j-1)+y5(j-2));
y9(j)=1/(1+C9)*(2*beta9*y9(j-1)-(1-C9)*y9(j-2)+y7(j)-2*beta9*y7(j-1)+y7(j-2));
end
ynotch=y9;
figure(23)
plot(n*T,ynotch,n*T,data-ynotch)
legend('Notch filtered HRECG signal','Error')
title('Notch Filtered Signal'); xlabel('Time (seconds)');ylabel('Normalized amplitude')
axis([0 2.7 -0.2 1.2])
% Magn & Ph of FFT of notch filtered sig %
ynotchfft = fft(ynotch,N);
mynotch=20*log(abs(ynotchfft));% FFT magnitude in dBs
pynotch= unwrap(angle(ynotchfft)); % FFT Phase
figure(24)
plot(xaxis(1:N/2),mynotch(1:N/2),xaxis(1:N/2),mdata(1:N/2))% for plotting +ive freq only
legend('Frequency spectrum of Notch filter output','Frequency spectrumof clean signal')
title('Frequency Spectrum of Notch Filtered Signal');
xlabel('Frequency (Hz)');ylabel('Magnitude (dB)'); %Grid on
axis([0 500 -60 110])
%% Time domain comparison
figure(25)
%  plot(n*T,yfilter,n*T,yest-PLI)
 plot(n*T,yest-PLI)

legend('SSRLS filter output', 'Error')
title('Filtered Signal Comparison')
xlabel('Time (seconds)');ylabel('Normalized Amplitude');%Grid on
% axis([0 5 -0.2 1.2])
axis([0 5 -0.1 0.1])
figure(31)% Output 
subplot(4,1,1),plot(n*T,data),axis([0.5 4.5 -0.2 1.5]);legend('Clean HRECG');title('Filtration comparison')%xlabel('Time (seconds)');ylabel('Normalized Amplitude');title('Pure HRECG signal')
subplot(4,1,2),plot(n*T,ecg),axis([0.5 4.5 -0.2 1.5]);legend('PLI corrupted signal');%xlabel('Time (seconds)');ylabel('Normalized Amplitude');title('Pure HRECG signal')
subplot(4,1,3),plot(n*T,yfilter),axis([0.5 4.5 -0.2 1.5]);ylabel('Normalized Amplitude');legend('SSRLS filter output');%xlabel('Time (seconds)');title('SSRLS Filter output');
subplot(4,1,4),plot(n*T,ynotch),axis([0.5 4.5 -0.2 1.5]);xlabel('Time (seconds)');legend('Notch filter output');%ylabel('Normalized Amplitude');title('Notch filter output');
% figure(31)% Output 
% subplot(4,1,1),plot(n*T,data),axis([0 10 -0.3 1.3]);legend('Clean HRECG');title('Filtration comparison')%xlabel('Time (seconds)');ylabel('Normalized Amplitude');title('Pure HRECG signal')
% subplot(4,1,2),plot(n*T,ecg),axis([0 10 -0.3 1.3]);legend('PLI corrupted signal');%xlabel('Time (seconds)');ylabel('Normalized Amplitude');title('Pure HRECG signal')
% subplot(4,1,3),plot(n*T,yfilter),axis([0 10 -0.3 1.3]);ylabel('Normalized Amplitude');legend('SSRLS filter output');%xlabel('Time (seconds)');title('SSRLS Filter output');
% subplot(4,1,4),plot(n*T,ynotch),axis([0 10 -0.3 1.3]);xlabel('Time (seconds)');legend('Notch filter output');%ylabel('Normalized Amplitude');title('Notch filter output');
figure(32)
plot(n*T,data,n*T,yfilter,n*T,ynotch)
legend('Clean HRECG signal','SSRLS filter output', 'Notch filter output')
title('Filtered Signal Comparison')
xlabel('Time (seconds)');ylabel('Normalized Amplitude');%Grid on
axis([5.02 5.12 -0.03 0.04])
%%
figure(1111);
[ha, pos] = tight_subplot(4,1,[.01 .03],[.1 .05],[.1 .03]);
axes(ha(4));plot(n*T,ynotch,n*T,data-ynotch,'r','Linewidth',1.5);
%  set(gca,'xtick',[],'ytick',[])
axis([0 10 -.2 .8]);
xlabel('Times (Sec)');
% ylabel('Normalized Magnitude')
legend ('Adaptive notch filter output','error','Orientation','horizontal');
set(gca,'xtick',[])
title('Filteration Comparison for HRA-IEGM Signal')
axes(ha(1)); plot(n*T,ecg,'Linewidth',1.5);

% ylabel('Normalized Magnitude')
% xlabel('Samples');ylabel('Normalized Magnitude')  
legend ('PLI corrupted HRA-IEGM signal');
set(gca,'xtick',[])
axis([0 10 -.3 .8]);
ylabel('Normalized Magnitude')
axes(ha(3)); plot(n*T,EE,n*T,error,'r','Linewidth',1.5);
% xlabel('Samples');ylabel('Normalized Magnitude')
legend ('Sequentially operated SSLMSWAM','error','Orientation','horizontal'); 
set(gca,'xtick',[])
axis([0 10 -.2 .8]);

axes(ha(2)); plot(n*T,EE1,n*T,error1,'r','Linewidth',1.5);
% xlabel('Samples');
legend ('Proposed PD-SSLMSWAM','error','Orientation','horizontal');
%set(gca,'xtick',[])
axis([0 10 -.2 .8]);
% xlabel('Times (Sec)');
