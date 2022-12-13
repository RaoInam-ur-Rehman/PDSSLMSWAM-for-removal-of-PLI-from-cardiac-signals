clear all;close all;clc
%%
load HRECG_seq_SSLMSWAM.mat 'MSE';HR1= MSE;
load HRECG_PD_SSLMSWAM.mat 'MSE';HR2= MSE;
load HRECG_NotchF.mat 'MSE';HR3= MSE;
L1=length(HR1);
load UHFECG_seq_SSLMSWAM.mat 'MSE';UHF1= MSE;
load UHFECG_PD_SSLMSWAM.mat 'MSE';UHF2= MSE;
load UHFECG_NotchF.mat 'MSE';UHF3= MSE;
L2=length(UHF1);
load IEGM_seq_SSLMSWAM.mat 'MSE';IGM1= MSE;
load IEGM_PD_SSLMSWAM.mat 'MSE';IGM2= MSE;
load IEGM_NotchF.mat 'MSE';IGM3= MSE;
L3=length(IGM1);
para=[HR1;HR2;HR3;UHF1;UHF2;UHF3;IGM1;IGM2;IGM3];
G=[ones(L1,1);2*ones(L1,1);3*ones(L1,1);4*ones(L2,1);5*ones(L2,1);6*ones(L2,1);7*ones(L3,1);8*ones(L3,1);9*ones(L3,1)];
%%
mean_Seq = mean(HR1,'omitnan');
std_Seq = std(HR1,'omitnan');
mean_P = mean(HR2,'omitnan');
std_P = std(HR2,'omitnan');
mean_N = mean(HR3,'omitnan');
std_N = std(HR3,'omitnan');
HRECG = table(mean_Seq,std_Seq,mean_P,std_P,mean_N,std_N)

mean_Seq = mean(UHF1,'omitnan');
std_Seq = std(UHF1,'omitnan');
mean_P = mean(UHF2,'omitnan');
std_P = std(UHF2,'omitnan');
mean_N = mean(UHF3,'omitnan');
std_N = std(UHF3,'omitnan');
UHFECG = table(mean_Seq,std_Seq,mean_P,std_P,mean_N,std_N)


mean_Seq = mean(IGM1,'omitnan');
std_Seq = std(IGM1,'omitnan');
mean_P = mean(IGM2,'omitnan');
std_P = std(IGM2,'omitnan');
mean_N = mean(IGM3,'omitnan');
std_N = std(IGM3,'omitnan');
IGM = table(mean_Seq,std_Seq,mean_P,std_P,mean_N,std_N)
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
%%
% positions = [1 1.25 1.5 2 2.25 2.5 3 3.25 3.5];
% figHandler = figure;
% hold on
% boxplot(para,G, 'positions', positions);
% set(gca,'xtick',[mean(positions(1:3))  mean(positions(4:6)) mean(positions(7:9))])
% set(gca,'xticklabel',{'HRECG','UHFECG','IEGM'})
% title('Mean Square Error for All Cardiac Signals')
% xlabel('Types of Cardiac Signals')
% ylabel('Mean Square Error')
% ylim([-1e-4 200e-3])
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
