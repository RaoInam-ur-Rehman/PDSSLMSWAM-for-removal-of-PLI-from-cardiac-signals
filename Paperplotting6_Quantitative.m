clear all; close all; clc;
%% Data loading UHF-ECG
% N=50000; fs=5000;
% T=1/fs; n=0:N-1;
% load ECG3
% data1=ecg(:,3)';
% data=data1-mean(data1);
% data=data./max(data);
% data=data(1:N);
% data=movmean(data,257);
%% Data loading HRECG
fs=1000;T=1/fs;N=10000;
n=0:1:N-1;
load ecg;
data1=ecg_normal_simulator(:,3);
data=data1-mean(data1); data=data(1:N);
data=1.5*data+.074;%for peak to peak normalization to 1 of 3rd clmn of ecg_normal_simulator sig
data=data';
%% fixed freq PLI disturbance
PLIfreq=50;
%% Varying Input Signal SNR
SIR=-25:5:25;
L=length(SIR);
SIRin=zeros(1,L);
for i=1:L
    SP=(1/N)*sum(data.^2);  % Signal Power
    NP=SP/(10.^(SIR(i)/10));
    PLI1=sin(2*pi*PLIfreq*T*n);PLI3=sin(2*pi*3*PLIfreq*T*n); 
    PLI5=sin(2*pi*5*PLIfreq*T*n);PLI7=sin(2*pi*7*PLIfreq*T*n);
    PLI9=sin(2*pi*9*PLIfreq*T*n);

    A1=(sqrt((.89)*2*NP));A3=(sqrt(2*NP*(.05)));
    A5=(sqrt(2*NP*(.03)));A7=(sqrt(2*NP*(.02))); A9=(sqrt(2*NP*(.01)));
    PLI1=abs(A1).*PLI1;PLI3=abs(A3).*PLI3;
    PLI5=abs(A5).*PLI5;PLI7=abs(A7).*PLI7; PLI9=abs(A9).*PLI9;
    PLI=PLI1+PLI3+PLI5+PLI7+PLI9;
%     PLI=sqrt(2*SP/(10.^(SIR(i)/10))).*PLI;
    ecg_con=data+PLI;
%     ecg=ecg(100:end);
    %%%% Tocheck and verify SIR
    NP1=(1/N)*sum(PLI.^2);
    AAA=(SP)/NP1;
    SIRin(i)=10*log10(AAA);
%% Call Functions
    freqest=PLIfreq*ones(1,N);
    mu1=0.01; %for UHF-ECG signal aplha=0.000001
    SysOrd=10;
    [yest1,EE1] =SSLMSwithvariant(ecg_con,fs,freqest,mu1,SysOrd,'SSLMSWAM');
    [yest2,EE2] =SSLMSwithvariant_PD(ecg_con,fs,freqest,mu1,SysOrd,'SSLMSWAM');
    ynotch = notchfilter(ecg_con',fs, N,freqest,100);
    CleanECG=data;
    error1=CleanECG-EE1;
    error2=CleanECG-EE2;
    error3=CleanECG-ynotch;
    %% Performance Parameters
    %%%% LUI paper prformance parameters
    gamma_lui1(i)=20*log10(norm(ecg_con)/norm(EE1));
    gamma_lui2(i)=20*log10(norm(ecg_con)/norm(EE2));
    gamma_lui3(i)=20*log10(norm(ecg_con)/norm(ynotch));
    MSE_lui1(i)=10*log10(mean(error1.^2));
    MSE_lui2(i)=10*log10(mean(error2.^2));
    MSE_lui3(i)=10*log10(mean(error3.^2));
    %%%% bahaz paper prformance parameters
    SNR_bahaz1(i)=10*log10(var(CleanECG)/var(error1));
    SNR_bahaz2(i)=10*log10(var(CleanECG)/var(error2));
    SNR_bahaz3(i)=10*log10(var(CleanECG)/var(error3));
    MSE_bahaz1(i)=mean(error1.^2);
    MSE_bahaz2(i)=mean(error2.^2);
    MSE_bahaz3(i)=mean(error3.^2);
    %%%%harrach paper prformance parameters
%     RMSE_harrach1(i)=rms(error1/CleanECG);
%     RMSE_harrach2(i)=rms(error2/CleanECG);
    %%%% SEdgo paper prformance parameters 
%     R_SEdgo1=EE1-CleanECG;
%     R_SEdgo2=EE2-CleanECG;
%     SIR_SEdgo1(i)=10*log10(mean(CleanECG.^2)/mean(R_SEdgo1.^2));
%     SIR_SEdgo2(i)=10*log10(mean(CleanECG.^2)/mean(R_SEdgo2.^2));
    PCorrC_sedgo1=corrcoef(CleanECG,EE1);
    PCorrC_sedgo2=corrcoef(CleanECG,EE2);
    PCorrC_sedgo3=corrcoef(CleanECG,ynotch);
    Pcorr1(i)=PCorrC_sedgo1(1,2);
    Pcorr2(i)=PCorrC_sedgo2(1,2);
    Pcorr3(i)=PCorrC_sedgo3(1,2);
    %%%% YQin paper prformance parameters
%     PRD_YQin1(i)=sqrt(sum(error1.^2)/sum(CleanECG.^2));
%     PRD_YQin2(i)=sqrt(sum(error2.^2)/sum(CleanECG.^2));
end
%% Plotting
figure;
plot(SIRin,gamma_lui1,'o-','Linewidth',1.5);
hold on;
plot(SIRin,gamma_lui2,'g--','Linewidth',1.5);
hold on;
plot(SIRin,gamma_lui3,'r--','Linewidth',1.5);
xlabel('SNR input [dB]');ylabel('Suppresion Ratio [dB]')
title ('Suppresion Ratio ');
legend('Sequntially operated SSLMSWAM','Proposed PD-SSLMSWAM','Notch Filter');
grid on;
figure;
plot(SIRin,SNR_bahaz1,'o-','Linewidth',1.5);
hold on;
plot(SIRin,SNR_bahaz2,'g--','Linewidth',1.5);
hold on;
plot(SIRin,SNR_bahaz3,'r--','Linewidth',1.5);
xlabel('SNR input [dB]');ylabel('SNR output [dB]')
title ('Output Signal to Noise Ratio ');
legend('Sequntially operated SSLMSWAM','Proposed PD-SSLMSWAM','Notch Filter');
grid on;
figure;
semilogy(SIRin,MSE_bahaz1,'o-','Linewidth',1.5);
hold on;
semilogy(SIRin,MSE_bahaz2,'g--','Linewidth',1.5);
hold on;
semilogy(SIRin,MSE_bahaz3,'r--','Linewidth',1.5);
xlabel('SNR input [dB]');ylabel('MSE')
title ('Mean Square Error ');
legend('Sequntially operated SSLMSWAM','Proposed PD-SSLMSWAM','Notch Filter');
grid on;
figure;
plot(SIRin,Pcorr1,'o-','Linewidth',1.5);
hold on;
plot(SIRin,Pcorr2,'g--','Linewidth',1.5)
hold on;
plot(SIRin,Pcorr3,'r--','Linewidth',1.5);
xlabel('SNR input [dB]');ylabel('Pcorr [dB]')
title ('Correlation Coefficient');
legend('Sequntially operated SSLMSWAM','Proposed PD-SSLMSWAM','Notch Filter');
grid on;