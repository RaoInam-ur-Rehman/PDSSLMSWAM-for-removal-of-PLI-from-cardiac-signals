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
%% Frequency estimation 
WinLen=fs/4;
Frange=[45 55];
[Fest1,FREQMAX1,Iter1, freqest1] = SAIDFT(ecg, fs, WinLen, Frange);
% mn=fix(WinLen/2); % defining move step of sliding window in each step (with round off)
% for m=1:mn:N-WinLen
%  ecg_N=ecg(1,m:m-1+WinLen);
% [Fest,FREQMAX,Iter] = SAsDFTwithvariant(ecg_N, fs, WinLen, Frange,'rSDFT');
% % [Fest1,FREQMAX1] = AdapDFTSingle(PL_N, fs, WinLen, Frange, Fbins, res);
% freqest(m+WinLen:m+mn+WinLen)=Fest;
% % freqest1(m+WinLen:m+mn+WinLen)=Fest1;
% end
%% for plotting
figure; stairs(1:Iter1,FREQMAX1,'-.or','LineWidth',2,'Marker','d','MarkerFaceColor','c'); 
ylim([45 55]);xlim([1 Iter1]); grid on;
xlabel('No. of Iterations');ylabel('Selected Frequency Band'); title('Convergence of Proposed SADFT')

PLI_freq=PLIfreq*ones(1,N);
figure(111)
plot(n*T,PLI_freq,'r',n*T,freqest1(1:N),'b','LineWidth',1.5)
legend('Actual fundamental frequency of PLI','Ongoing Frequency estimation by Proposed SADFT')
xlabel('Time (seconds)'), ylabel('Frequency (Hz)')
ylim([PLIfreq-1 PLIfreq+1])