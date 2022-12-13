clear all; close all; clc
%% load ECG signal %%%
Q=1000;
load ecg
% data1=ecg% Q=1000;
data1=ecg_normal_simulator(:,3);
data=data1-mean(data1);
x=1.5*data+.074;%for peak to peak normalization to 1 of 3rd clmn of ecg_normal_simulator sig
LX=length(x); if LX>5*Q; LX=5*Q; end
x=(x(1:LX)'); % Signal length
%% Data loading UHF-ECG
% N=50000; fs=5000;
% Q=5000;
% T=1/fs; n=0:N-1;
% load ECG3
% data1=ecg(:,3);
% data=data1-mean(data1);
% data=data./max(data);
% % data=data(1:N);
% data=movmean(data,257);
% LX=length(data); if LX>10*Q; LX=10*Q; end
% x=(data(1:LX)'); % Signal length
%% Data HRA IEGM
%  N=10000; fs=2000; 
% % T=1/fs; n=0:N-1;
% Q=2000;
% load parveen
% daata=data;
% data1=daata(1:N,4);
% data2=data1-mean(data1);
% data2=data2./max(data2);
% data3=movmean(data2,10);
% LX=length(data2); if LX>10*Q; LX=10*Q; end
% x=(data3(1:LX)'); % Signal length
%%
t = 1/Q;
 a=zeros(1,LX); am=0.35; %%% Amplitude deviation %%%
 for i=1: 1: LX
 a(i) = am; % constant
%  a(i) = 2*am/LX*i; % linear low
%  a(i) = abs(am*cos(2*pi*0.16*i*t)-am)/2; % sinusoidal low
%  if i<LX/2; a(i)=am/2; else a(i)=2*am; end % abrupt change
 end
 dPL=zeros(1,LX); %%% Frequency deviation %%%
 F=50; dF=F*0.015;
for i=1: 1: LX
 dPL(i)=dF; % constant
%  dPL(i)=i*t*dF/8-dF; % linear low
%  dPL(i)=dF*sin(2*pi*0.0625*i*t); % sinusoidal low
%  if i<LX/2; dPL(i)=+dF; else dPL(i)=-dF; end % abrupt change
 end
 %%% Interference composing %%%
 PL=x*0; iF=0;
 for i=1: 1: LX
 iF=iF+dPL(i)*t;
%  PL(i) = a(i)*sin(2*pi*(F+iF/i/t)*i*t);
PL1(i) = a(i)*sin(2*pi*(F+iF/i/t)*i*t);
PL3(i) = (1/5)*a(i)*sin(2*pi*3*(F+iF/i/t)*i*t);
PL5(i) =(1/10)*a(i)*sin(2*pi*5*(F+iF/i/t)*i*t);
PL7(i) = (1/50)*a(i)*sin(2*pi*7*(F+iF/i/t)*i*t);
PL9(i) = (1/100)*a(i)*sin(2*pi*9*(F+iF/i/t)*i*t);
 end
% PL=(PL1+PL3);
PL=(PL1+PL3+PL5+PL7+PL9);
X = x+PL; % Add mains interference
%% Frequency estimation 
% WinLen=Q/4;
% Frange=[45 55];
% Fbins=10;
% res=5;
% [Fest1,FREQMAX1,Iter1, freqest1] = SAIDFT(X, Q, WinLen, Frange);
%  [Fest,FREQMAX,freqest1] = AdapDFT(X, Q, WinLen, Frange, Fbins, res);

%%
% freqest=PLIfreq*ones(1,N);
freqest=F+dPL;
fs=Q;
% mu1=0.01; %for HRECG signal aplha=0.00001
mu1=0.01; %for UHF-ECG signal aplha=0.000001
% mu2=0.09;
SysOrd=10;
[yest,EE] =SSLMSwithvariant(X,fs,freqest,mu1,SysOrd,'SSLMSWAM');
 error=x-EE;
[yest1,EE1] =SSLMSwithvariant_PD(X,fs,freqest,mu1,SysOrd,'SSLMSWAM');
 error1=x-EE1;
%% adaptive notch filter
data=x';
ecg=X';
N=LX;
ynotch = notchfilter(ecg,fs, N,freqest,100);
%%
% q=100;% quality factor
% dftwindow=fs/2;
% % K=N;% plotting FFT
% % datafft=fft(data,N);
% % mdata=20*log(abs(datafft));% magnitude response of FFT
% % xaxis= (1:N).*fs/K;% defining x axis in Hz
% 
% y1=zeros(N,1);y3=zeros(N,1);y5=zeros(N,1);y7=zeros(N,1);y9=zeros(N,1);
% for j=3:N
%     if j<dftwindow
%         fo=50;
%     else
%         fo=freqest(j);
%     end
% omega1=2*pi*fo/fs;C1=tan(0.5*omega1/q);beta1=cos(omega1);%for 1st harmonic
% omega3=2*pi*3*fo/fs;C3=tan(0.5*omega3/q);beta3=cos(omega3);%for 3rd harmonic
% omega5=2*pi*5*fo/fs;C5=tan(0.5*omega5/q);beta5=cos(omega5);%for 5th harmonic
% omega7=2*pi*7*fo/fs;C7=tan(0.5*omega7/q);beta7=cos(omega7);%for 7th harmonic
% omega9=2*pi*9*fo/fs;C9=tan(0.5*omega9/q);beta9=cos(omega9);%for 9th harmonic
% 
% y1(j)=1/(1+C1)*(2*beta1*y1(j-1)-(1-C1)*y1(j-2)+ecg(j)-2*beta1*ecg(j-1)+ecg(j-2));
% y3(j)=1/(1+C3)*(2*beta3*y3(j-1)-(1-C3)*y3(j-2)+y1(j)-2*beta3*y1(j-1)+y1(j-2));
% y5(j)=1/(1+C5)*(2*beta5*y5(j-1)-(1-C5)*y5(j-2)+y3(j)-2*beta5*y3(j-1)+y3(j-2));
% y7(j)=1/(1+C7)*(2*beta7*y7(j-1)-(1-C7)*y7(j-2)+y5(j)-2*beta7*y5(j-1)+y5(j-2));
% y9(j)=1/(1+C9)*(2*beta9*y9(j-1)-(1-C9)*y9(j-2)+y7(j)-2*beta9*y7(j-1)+y7(j-2));
% end
% ynotch=y9;
% ynotch=ynotch';
 %%
tt = 0: 1/Q: LX/Q;
axe=1:LX;
%  Hh1=[0 LX*t min(x)-max(a) max(x)+max(a)];
%  Hh1=[0 LX*t -0.6 1.5]; % For HRECG input
%  Hh2=[0 2 -.4 .6];      % For HRECG input
%  Hh1=[0 LX*t -0.6 1.7]; % For UHFECG input
%  Hh2=[0 2 -.7 .9];      % For UHFECG input
 
 Hh1=[0 LX*t -0.6 1.2]; % For HRA-IEGM input  
 Hh2=[0 2 -.4 .5];  % For HRA-IEGM input
figure(11);
[ha, pos] = tight_subplot(5,1,[.04 .03],[.06 .05],[.06 .03]);
axes(ha(1)); plot(tt(axe),X(axe),'LineWidth',1.5);axis(Hh1);
legend ('PLI corrupted cardiac signal');
% axes(ha(1)); plot(tt(axe),X(axe),tt(axe),a(axe),'g','LineWidth',1.5);axis(Hh1);
% legend ('PLI corrupted cardiac signal','Amplitude deviation','Orientation','horizontal');
set(gca,'xtick',[])
% title('Filtration Comparison with Abrupt Amplitude Deviation')
title('Filtration Comparison for HRECG Signal')
axes(ha(2)); plot(tt(axe),EE1(axe),'LineWidth',1.5);axis(Hh1);
legend ('Filtered Output (Proposed PD-SSLMSWAM)'); 
set(gca,'xtick',[])
axes(ha(3)); plot(tt(axe),EE(axe),'LineWidth',1.5);axis(Hh1);
legend ('Filtered Output (Sequential SSLMSWAM)'); 
set(gca,'xtick',[])
ylabel('Normalized Magnitude')
axes(ha(4)); plot(tt(axe),ynotch(axe),'LineWidth',1.5);axis(Hh1);
legend ('Filtered Output (Notch Filter)'); 
% set(gca,'xtick',[])
axes(ha(5)); plot(tt(axe),error(axe),'g--',tt(axe),error1(axe),'b--',tt(axe),x-ynotch(axe),'r','LineWidth',1);axis(Hh2);
legend ('Zoomed error (sequential)','Zoomed error (Proposed)','Zoomed error (Notch Filter)','Orientation','horizontal');
axis(Hh2);
xlabel('Times (seconds)');