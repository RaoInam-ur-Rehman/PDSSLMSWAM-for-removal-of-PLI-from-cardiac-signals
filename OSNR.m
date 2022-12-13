clear all;close all;clc
%%
load HRECG_seq_SSLMSWAM.mat 'SNR';HR1= SNR;
load HRECG_PD_SSLMSWAM.mat 'SNR';HR2= SNR;
load HRECG_NotchF.mat 'SNR';HR3= SNR;
L1=length(HR1);
load UHFECG_seq_SSLMSWAM.mat 'SNR';UHF1= SNR;
load UHFECG_PD_SSLMSWAM.mat 'SNR';UHF2= SNR;
load UHFECG_NotchF.mat 'SNR';UHF3= SNR;
L2=length(UHF1);
load IEGM_seq_SSLMSWAM.mat 'SNR';IGM1= SNR;
load IEGM_PD_SSLMSWAM.mat 'SNR';IGM2= SNR;
load IEGM_NotchF.mat 'SNR';IGM3= SNR;
L3=length(IGM1);
para=[HR1;HR2;HR3;UHF1;UHF2;UHF3;IGM1;IGM2;IGM3];
G=[ones(L1,1);2*ones(L1,1);3*ones(L1,1);4*ones(L2,1);5*ones(L2,1);6*ones(L2,1);7*ones(L3,1);8*ones(L3,1);9*ones(L3,1)];
%%
var_Seq = var(HR1,'omitnan');
var_P = var(HR2,'omitnan');
var_N = var(HR3,'omitnan');
var_HRECG = table(var_Seq,var_P,var_N)

var_Seq = var(UHF1,'omitnan');
var_P = var(UHF2,'omitnan');
var_N = var(UHF3,'omitnan');
var_UHFECG = table(var_Seq,var_P,var_N)

var_Seq = var(IGM1,'omitnan');
var_P = var(IGM2,'omitnan');
var_N = var(IGM3,'omitnan');
var_IGM = table(var_Seq,var_P,var_N)
% %%
% positions = [1 1.25 1.5 2 2.25 2.5 3 3.25 3.5];
% figHandler = figure;
% hold on
% boxplot(para,G, 'positions', positions);
% set(gca,'xtick',[mean(positions(1:3))  mean(positions(4:6)) mean(positions(7:9))])
% set(gca,'xticklabel',{'HRECG','UHFECG','IEGM'})
% title('Mean Square Error for All Cardiac Signals')
% xlabel('Types of Cardiac Signals')
% ylabel('Mean Square Error')
% ylim([-250 40])
% color = ['c', 'y', 'g', 'c', 'y', 'g','c', 'y', 'g'];
% h = findobj(gca,'Tag','Box');
% for j=1:length(h)
%    patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
% end
% c = get(gca, 'Children');
% hleg1 = legend(c(1:3), 'Sequentially operated SSLMSWAM', 'Proposed PD-SSLMSWAM','Notch Filter' );
% hold off; 
% magnifyOnFigure(...
%         figHandler,...
%         'units', 'pixels',...
%         'magnifierShape', 'rectangle',...
%         'initialPositionSecondaryAxes', [230.933 185.189 264.941 100.65],...
%         'initialPositionMagnifier',     [55.56 55.2 250.6823 40.519],...    
%         'mode', 'interactive',...    
%         'displayLinkStyle', 'straight',...        
%         'edgeWidth', .5,...
%         'edgeColor', 'black',...
%         'secondaryAxesFaceColor', [1 1 1]... 
%             );   
% 
