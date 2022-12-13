clear all; close all; clc
%% load ECG signal %%%
% Q=1000;
% load ecg;
% for l=1:1:9
%      data11=eval(['ecg',int2str(l)]);
% for k=1:1:9
% %     data1=ecg_normal_simulator(:,k);
%     data1=data11(:,k);
%     data=data1-mean(data1);
%     x=1.5*data+.074;%for peak to peak normalization to 1 of 3rd clmn of ecg_normal_simulator sig
%     LX=length(x); if LX>10*Q; LX=10*Q; end
%     x=(x(1:LX)); % Signal length
%     t = 1/Q;
    %%
Q=5000;
load ECG3;
M=size(ecg,2);
Lx=10*Q;
for l=1:M
    for k=1:30
        data11(:,k)=ecg(Lx*(k-1)+1:k*Lx,l);
         data1=data11(:,k);
    data=data1-mean(data1);
    data=data./max(data);
    x=movmean(data,257);
%     x=1.5*data+.074;%for peak to peak normalization to 1 of 3rd clmn of ecg_normal_simulator sig
    LX=length(x); if LX>10*Q; LX=10*Q; end
    x=(x(1:LX)); % Signal length
    t = 1/Q;
%%
% Q=2000;
% load parveen;
% M=size(data,2);
% ecg=data;
% Lx=10*Q;
% for l=1:M
%     for k=1:2
%         data11(:,k)=ecg(Lx*(k-1)+1:k*Lx,l);
%          data1=data11(:,k);
%     data=data1-mean(data1);
%     data=data./max(data);
%     x=movmean(data,257);
% %     x=1.5*data+.074;%for peak to peak normalization to 1 of 3rd clmn of ecg_normal_simulator sig
%     LX=length(x); if LX>10*Q; LX=10*Q; end
%     x=(x(1:LX)); % Signal length
%     t = 1/Q;
    %%
a=zeros(1,LX); am=0.15; %%% Amplitude deviation %%%
 for i=1: 1: LX
 a(i) = am; % constant
%  a(i) = am/LX*i; % linear low
%  a(i) = abs(am*cos(2*pi*0.06*i*t)-am)/2; % sinusoidal low
%  if i<LX/2; a(i)=am; else a(i)=0; end % abrupt change
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
%   PL(i) = a(i)*sin(2*pi*(F+iF/i/t)*i*t);
PL1(i) = a(i)*sin(2*pi*(F+iF/i/t)*i*t);
PL3(i) = (1/5)*a(i)*sin(2*pi*3*(F+iF/i/t)*i*t);
PL5(i) =(1/10)*a(i)*sin(2*pi*5*(F+iF/i/t)*i*t);
PL7(i) = (1/50)*a(i)*sin(2*pi*7*(F+iF/i/t)*i*t);
PL9(i) = (1/100)*a(i)*sin(2*pi*9*(F+iF/i/t)*i*t);
 end
% PL=(PL1+PL3);
PL=(PL1+PL3+PL5+PL7+PL9);
X = x+PL'; % Add mains interference
%% Call Algorithms
%     freqest=PLIfreq*ones(1,N);
    freqest=F+dPL;
    mu1=0.01; %for UHF-ECG signal aplha=0.000001
    SysOrd=10;
%      [yest,EE] =SSLMSwithvariant(X',Q  ,freqest,mu1,SysOrd,'SSLMSWAM');
%     [yest,EE] =SSLMSwithvariant_PD(X',Q,freqest,mu1,SysOrd,'SSLMSWAM');
    EE = notchfilter(X,Q, LX,freqest,100);
    error=x-EE';
%%% Call Subtraction procedure for Power-Line Interference removing
%   M=0.1; res=0.002; [Y]=PLinterference_removing_Vd2(X,res,Q,F,M);
 % M=.1; res=0.002; [Y]=PLinterference_removing_Vd21111(X,res,Q,F,M);
%%%%%%% Ptroposed Algorithm %%%%%%%%%%5
% lemda=0.9999;
% freqest=F+dPL;
% %     [yest1,EE1] =ssrlsSeq(ecg',fs,lemda,freqest);
%     [yest2,EE2] =ssrlsPD(X',Q,lemda,freqest);
% %      error1=x-Y;
%     error2=x-EE2';
 %%% Errors calculating %%%
 %%%% gamma    
% gamma_lui1(l,k)=20*log10(norm(X)/norm(Y));
gamma_lui2(l,k)=20*log10(norm(X)/norm(EE));
%%% SNR output & MSE bahaz
% SNR_bahaz1(l,k)=10*log10(var(x)/var(error1));
SNR_bahaz2(l,k)=10*log10(var(x)/var(error));
% MSE_bahaz1(l,k)=mean(error1.^2);
MSE_bahaz2(l,k)=mean(error.^2);
%%% Pcorr
% PCorrC_sedgo1=corrcoef(x,Y);
PCorrC_sedgo2=corrcoef(x,EE);
% Pcorr1(l,k)=PCorrC_sedgo1(1,2);
Pcorr2(l,k)=PCorrC_sedgo2(1,2);
end
end
%%
MSE=reshape(MSE_bahaz2,[],1);
% MSE_PDSSRLS=reshape(MSE_bahaz2,[],1);
SNR=reshape(SNR_bahaz2,[],1);
% SNR_PDSSRLS=reshape(SNR_bahaz2,[],1);
gamma=reshape(gamma_lui2,[],1);
% gamma_PDSSRLS=reshape(gamma_lui2,[],1);
Pcorr=reshape(Pcorr2,[],1);
% Pcorr_PDSSRLS=reshape(Pcorr2,[],1);
save('UHFECG_NotchF','MSE','SNR','gamma','Pcorr')
% save('UHFECG_seq_SSLMSWAM','MSE','SNR','gamma','Pcorr')
% save('UHFECG_PD_SSLMSWAM','MSE','SNR','gamma','Pcorr')
% save('IEGM.mat','MSE_SP','MSE_PDSSRLS','SNR_SP','SNR_PDSSRLS','gamma_SP','gamma_PDSSRLS','Pcorr_SP','Pcorr_PDSSRLS')
% save('HRECG.mat','MSE_SP','MSE_PDSSRLS','SNR_SP','SNR_PDSSRLS','gamma_SP','gamma_PDSSRLS','Pcorr_SP','Pcorr_PDSSRLS')
% save('UHFECG.mat','MSE_SP','MSE_PDSSRLS','SNR_SP','SNR_PDSSRLS','gamma_SP','gamma_PDSSRLS','Pcorr_SP','Pcorr_PDSSRLS')
