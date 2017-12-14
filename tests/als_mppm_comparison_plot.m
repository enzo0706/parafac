% Post-analysis plot for the ALS-MPPM algorithm of different configurations.
%
%  ----------    --      --    --------          --        --      --
% |          |  |  |    |  |  |        \        /  \      |  \    |  |
% |   -------   |  |    |  |  |   ---\  \      /    \     |   \   |  |
% |  |          |  |    |  |  |  |    \  \    /      \    |    \  |  |
% |   -------   |  |    |  |  |  |    |  |   /   /\   \   |     \ |  |
% |          |  |  |    |  |  |  |    |  |  |   /  \   |  |      \|  |
% |   -------   |  |    |  |  |  |    |  |  |   ----   |  |  |\      |
% |  |          |  |    |  |  |  |    |  |  |          |  |  | \     |
% |  |          |  |    |  |  |  |    /  /  |   ----   |  |  |  \    |
% |  |          \  \----/  /  |   ---/  /   |  |    |  |  |  |   \   |
% |  |           \        /   |        /    |  |    |  |  |  |    \  |
%  --             --------     --------      --      --    --      --
%  
% |    | |\    |  ---  |    |  -----  ----   ----   ---  ----- \    /
% |    | | \   |   |   |    | |      |    | |        |     |    \  /
% |    | |  \  |   |   \    / |----- |----   ----    |     |     \/
% |    | |   \ |   |    \  /  |      |  \        |   |     |     /
%  ----  |    \|  ---    \/    ----- |   \   ----   ---    |    /
% 
% ----------------------------------------------------------------------
% 
%  ----------    --      --    ----------   |           /\     ------
% |          |  |  \    /  |  |          |  |          /  \   |      \
% |   -------   |   \  /   |   ---    ---   |         /    \  |      |
% |  |          |    \/    |      |  |      |        /      \ |      /
% |   -------   |          |      |  |      |        |      | |------
% |          |  |  |\  /|  |      |  |      |        |------| |      \
% |   -------   |  | \/ |  |      |  |      |        |      | |      |
% |  |          |  |    |  |      |  |      |        |      | |      /
% |   -------   |  |    |  |   ---    ---    ------- |      |  ------
% |          |  |  |    |  |  |          |  
%  ----------    --      --    ----------   
%
% History:
%   2017.12.24 initial committed by Ming Shi
%

clear all
close all
clc

[folder,~,~] = fileparts(which(mfilename));
load([folder,'/results/ALS-MPPM_COMPARISON_RESULT.mat'])

h = figure(1);
h.Units = 'inches';
h.Position = [0.05,0.05,3,2.5];
p(1) = semilogy(SNR,Iter(:,1,1),'LineWidth',1.5,'LineStyle','-','Color',[0.7,0,0],'Marker','None'); hold on, grid on
p(2) = semilogy(SNR,Iter(:,2,1),'LineWidth',1.5,'LineStyle','-','Color',[0,0.7,0],'Marker','x');
p(3) = semilogy(SNR,Iter(:,3,1),'LineWidth',1.5,'LineStyle','--','Color',[0,0,0.7],'Marker','None');
p(4) = semilogy(SNR,Iter(:,4,1),'LineWidth',1.5,'LineStyle','--','Color',[0.7,0.7,0],'Marker','o');
p(5) = semilogy(SNR,Iter(:,5,1),'LineWidth',1.5,'LineStyle','-','Color',[0.7,0,0.7],'Marker','s');
axis tight
xlabel('SNR (dB)','FontName','Times New Roman','FontSize',7)
ylabel('Iterations','FontName','Times New Roman','FontSize',7)
set(gca,'FontName','Times New Roman','FontSize',7)
gridLegend(p,3,Alg,'FontSize',7,'Box','off');

h = figure(2);
h.Units = 'inches';
h.Position = [0.05,0.05,3,2.5];
p(1) = semilogy(SNR,Err(:,1,1),'LineWidth',1.5,'LineStyle','-','Color',[0.7,0,0],'Marker','None'); hold on, grid on
p(2) = semilogy(SNR,Err(:,2,1),'LineWidth',1.5,'LineStyle','-','Color',[0,0.7,0],'Marker','x');
p(3) = semilogy(SNR,Err(:,3,1),'LineWidth',1.5,'LineStyle','--','Color',[0,0,0.7],'Marker','None');
p(4) = semilogy(SNR,Err(:,4,1),'LineWidth',1.5,'LineStyle','--','Color',[0.7,0.7,0],'Marker','o');
p(5) = semilogy(SNR,Err(:,5,1),'LineWidth',1.5,'LineStyle','-','Color',[0.7,0,0.7],'Marker','s');
axis tight
xlabel('SNR (dB)','FontName','Times New Roman','FontSize',7)
ylabel('NMSE','FontName','Times New Roman','FontSize',7)
set(gca,'FontName','Times New Roman','FontSize',7)
gridLegend(p,3,Alg,'FontSize',7,'Box','off');

h = figure(3);
h.Units = 'inches';
h.Position = [0.05,0.05,3,2.5];
p(1) = semilogy(SNR,Time(:,1,1),'LineWidth',1.5,'LineStyle','-','Color',[0.7,0,0],'Marker','None'); hold on, grid on
p(2) = semilogy(SNR,Time(:,2,1),'LineWidth',1.5,'LineStyle','-','Color',[0,0.7,0],'Marker','x');
p(3) = semilogy(SNR,Time(:,3,1),'LineWidth',1.5,'LineStyle','--','Color',[0,0,0.7],'Marker','None');
p(4) = semilogy(SNR,Time(:,4,1),'LineWidth',1.5,'LineStyle','--','Color',[0.7,0.7,0],'Marker','o');
p(5) = semilogy(SNR,Time(:,5,1),'LineWidth',1.5,'LineStyle','-','Color',[0.7,0,0.7],'Marker','s');
axis tight
xlabel('SNR (dB)','FontName','Times New Roman','FontSize',7)
ylabel('CPU Time','FontName','Times New Roman','FontSize',7)
set(gca,'FontName','Times New Roman','FontSize',7)
gridLegend(p,3,Alg,'FontSize',7,'Box','off');

h = figure(4);
h.Units = 'inches';
h.Position = [0.05,0.05,3,2.5];
p(1) = plot(SNR,Fail(:,1,1)*100,'LineWidth',1.5,'LineStyle','-','Color',[0.7,0,0],'Marker','None'); hold on, grid on
p(2) = plot(SNR,Fail(:,2,1)*100,'LineWidth',1.5,'LineStyle','-','Color',[0,0.7,0],'Marker','x');
p(3) = plot(SNR,Fail(:,3,1)*100,'LineWidth',1.5,'LineStyle','--','Color',[0,0,0.7],'Marker','None');
p(4) = plot(SNR,Fail(:,4,1)*100,'LineWidth',1.5,'LineStyle','--','Color',[0.7,0.7,0],'Marker','o');
p(5) = plot(SNR,Fail(:,5,1)*100,'LineWidth',1.5,'LineStyle','-','Color',[0.7,0,0.7],'Marker','s');
axis tight
xlabel('SNR (dB)','FontName','Times New Roman','FontSize',7)
ylabel('Prediction Failure (%)','FontName','Times New Roman','FontSize',7)
set(gca,'FontName','Times New Roman','FontSize',7)
gridLegend(p,3,Alg,'FontSize',7,'Box','off');

p = [];
h = figure(5);
h.Units = 'inches';
h.Position = [0.05,0.05,3,2.5];
p(1) = plot(SNR,Acce_ratio(:,2,1),'LineWidth',1.5,'LineStyle','-','Color',[0,0.7,0],'Marker','x');hold on, grid on
p(2) = plot(SNR,Acce_ratio(:,3,1),'LineWidth',1.5,'LineStyle','--','Color',[0,0,0.7],'Marker','None');
p(3) = plot(SNR,Acce_ratio(:,4,1),'LineWidth',1.5,'LineStyle','--','Color',[0.7,0.7,0],'Marker','o');
p(4) = plot(SNR,Acce_ratio(:,5,1),'LineWidth',1.5,'LineStyle','-','Color',[0.7,0,0.7],'Marker','s');
axis tight
xlabel('SNR (dB)','FontName','Times New Roman','FontSize',7)
ylabel('Accleration Ratio','FontName','Times New Roman','FontSize',7)
set(gca,'FontName','Times New Roman','FontSize',7)
gridLegend(p,2,Alg(2:end),'FontSize',7,'Box','off');