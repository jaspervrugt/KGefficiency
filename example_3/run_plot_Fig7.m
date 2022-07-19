%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% Vrugt, J.A. and de Oliveira, D.Y., 2022. Confidence intervals of the
%%  Kling-Gupta efficiency, Journal of Hydrology, 612, 
%%  10.1016/j.jhydrol.2022.127968
%% https://doi.org/10.1016/j.jhydrol.2022.127968
%%
%% Section 3.3.2. Application: a conceptual watershed model (Figure 7)
%% 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc, close all

% Add folders to path
addpath(genpath([pwd '\results_obs']))
addpath(genpath([pwd '\results_bootstrap']))

%% Initialize figure
fig = figure('Position',[100 100 600 300]);
fontsize = 12;

%% Leaf River near Collins, MS (USGS 02472000)

% Load results
load('SCE_obs_KG_187.mat'), KG_obs = KG;
load('SCE_KG_187.mat')

subplot(1,2,1)
h = histogram(KG(2:1001),'NumBins',15,'normalization','pdf');
h.FaceColor = 'w';    
hold on
plot(KG_obs(1),0,'rx','MarkerSize',10,'LineWidth',2)
text(0.01,0.95,'(A)','Units','normalized','interpreter','latex','fontsize',fontsize,'fontweight','normal')
drawnow

% Adjust figure properties
set(gca,'TickDir','out')
xlabel('KG efficiency','interpreter','latex','fontsize',fontsize);
ylabel('Density','interpreter','latex','fontsize',fontsize);
set(gca,'TickLabelInterpreter','latex','fontsize',fontsize)
    pos = get(gca,'Position');
    set(gca,'Position',[pos(1) 2*pos(2) 1.1*pos(3) 0.9*pos(4)])
ax = gca;
ax.TickLength = [0.025 0.025];
% set box property to off and remove background color
set(ax,'box','off','color','none')
drawnow

%% Kinchafoonee Creek near Dawson, GA (USGS 02350900)

% Load results
load('SCE_obs_KG_165.mat'), KG_obs = KG;
load('SCE_KG_165.mat')

subplot(1,2,2)
h = histogram(KG(2:1001),'NumBins',15,'normalization','pdf');
h.FaceColor = 'w';
hold on
plot(KG_obs(1),0,'rx','MarkerSize',10,'LineWidth',2)
text(0.01,0.95,'(B)','Units','normalized','interpreter','latex','fontsize',fontsize,'fontweight','normal')
drawnow

% Adjust figure properties
set(gca,'TickDir','out')
xlabel('KG efficiency','interpreter','latex','fontsize',fontsize);
set(gca,'TickLabelInterpreter','latex','fontsize',fontsize)
pos = get(gca,'Position');
set(gca,'Position',[pos(1) 2*pos(2) 1.1*pos(3) 0.9*pos(4)])
ax = gca;
ax.TickLength = [0.025 0.025];
% set box property to off and remove background color
set(ax,'box','off','color','none')
drawnow