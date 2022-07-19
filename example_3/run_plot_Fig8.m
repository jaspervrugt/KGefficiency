%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% Vrugt, J.A. and de Oliveira, D.Y., 2022. Confidence intervals of the
%%  Kling-Gupta efficiency, Journal of Hydrology, 612, 
%%  10.1016/j.jhydrol.2022.127968
%% https://doi.org/10.1016/j.jhydrol.2022.127968
%%
%% Section 3.3.2. Application: a conceptual watershed model (Figure 8)
%% 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc, close all

% Add folders to path
addpath(genpath([pwd '\results_obs']))
addpath(genpath([pwd '\results_bootstrap']))

%% Figure settings

% Font size
fontsize = 12;

% Labels
p = 0;
p = p + 1; parlabel{p} = '$I_\mathrm{max}$ (mm)';
p = p + 1; parlabel{p} = '$S_\mathrm{max}$ (mm)';
p = p + 1; parlabel{p} = '$Q_\mathrm{max}$ (mm/d)';
p = p + 1; parlabel{p} = '$\alpha_\mathrm{e}$ (-)';
p = p + 1; parlabel{p} = '$\alpha_\mathrm{f}$ (-)';
p = p + 1; parlabel{p} = '$K_\mathrm{f}$ (d)';
p = p + 1; parlabel{p} = '$K_\mathrm{s}$ (d)';
p = p + 1; parlabel{p} = 'KG statistic';

% 99 95 90 50
color_in = [ 0.3 0.3 0.3; 0.5 0.5 0.5 ; 0.7 0.7 0.7; 0.9 0.9 0.9]; 

% Initialize figure
fig = figure('Position',[100 100 900 450]);

%% Leaf River near Collins, MS (USGS 02472000)

% Load data
load('SCE_obs_xopt_187.mat'), x_MAP = xopt(1,:);
load('SCE_xopt_187.mat'), xopt = xopt(2:1001,:);
load('SCE_KG_187.mat'), KG(1) = [];

% Select parameters to plot
ip2 = [2 2 2 2];
ip1 = [1 5 6 7];

count = 1;
letters = {'(A1)';'(B1)';'(C1)';'(D1)';'(E1)'};
for ip=1:4
    
    subplot(2,4,count)
    idx = KG>prctile(KG,1)&KG<=prctile(KG,5);
    scatter(xopt(idx,ip1(ip)),xopt(idx,ip2(ip)),10,'markerfacecolor',color_in(1,:),...
        'markeredgecolor','none');
    hold on

    idx = KG>prctile(KG,5)&KG<=prctile(KG,10);
    scatter(xopt(idx,ip1(ip)),xopt(idx,ip2(ip)),10,'markerfacecolor',color_in(2,:),...
        'markeredgecolor','none');

    idx = KG>prctile(KG,10)&KG<=prctile(KG,50);
    scatter(xopt(idx,ip1(ip)),xopt(idx,ip2(ip)),10,'markerfacecolor',color_in(3,:),...
        'markeredgecolor','none');
    
    idx = KG>prctile(KG,50);
    scatter(xopt(idx,ip1(ip)),xopt(idx,ip2(ip)),10,'markerfacecolor',color_in(4,:),...
        'markeredgecolor','none');
    
    plot(x_MAP(ip1(ip)),x_MAP(ip2(ip)),'rx','MarkerSize',10,'LineWidth',2)
    
    % Adjust figure properties
    set(gca,'TickDir','out')
    if count == 1
        ylabel(parlabel(ip2(ip)),'interpreter','latex','fontsize',fontsize);
        yh1 = get(gca,'ylabel');            % handle to the label object
        py1 = get(yh1,'position');          % get the current position property
        xl = xlim;
        py1(1) = xl(1) - 0.3*(xl(2)-xl(1)); % double the distance, 
                                            % negative values put the label below the axis
        set(yh1,'position',py1);            % set the new position
    else
        set(gca,'YTickLabel',[])
    end
    set(gca,'TickLabelInterpreter','latex','fontsize',fontsize)
    pos = get(gca,'Position');
    set(gca,'Position',[pos(1) pos(2) 1.1*pos(3) pos(4)])
    text(0.02,0.10,letters(ip),'Units','normalized','interpreter','latex','fontsize',fontsize,'fontweight','normal')
    ax = gca;
    ax.TickLength = [0.035 0.035];
    % set box property to off and remove background color
    set(ax,'box','off','color','none')
    drawnow

    count = count + 1;

end

%% Kinchafoonee Creek near Dawson, GA (USGS 02350900)

% Load data
load('SCE_obs_xopt_165.mat'), x_MAP = xopt(1,:);
load('SCE_xopt_165.mat'), xopt = xopt(2:1001,:);
load('SCE_KG_165.mat'), KG(1) = [];

% Select parameters to plot
ip2 = [2 2 2 2];
ip1 = [1 5 6 7];

count = 1;
letters = {'(A2)';'(B2)';'(C2)';'(D2)';'(E2)'};
for ip=1:4

    subplot(2,4,4+count)
    idx = KG>prctile(KG,1)&KG<=prctile(KG,5);
    scatter(xopt(idx,ip1(ip)),xopt(idx,ip2(ip)),10,'markerfacecolor',color_in(1,:),...
        'markeredgecolor','none');
    hold on

    idx = KG>prctile(KG,5)&KG<=prctile(KG,10);
    scatter(xopt(idx,ip1(ip)),xopt(idx,ip2(ip)),10,'markerfacecolor',color_in(2,:),...
        'markeredgecolor','none');

    idx = KG>prctile(KG,10)&KG<=prctile(KG,50);
    scatter(xopt(idx,ip1(ip)),xopt(idx,ip2(ip)),10,'markerfacecolor',color_in(3,:),...
        'markeredgecolor','none');
    
    idx = KG>prctile(KG,50);
    scatter(xopt(idx,ip1(ip)),xopt(idx,ip2(ip)),10,'markerfacecolor',color_in(4,:),...
        'markeredgecolor','none');
    
    plot(x_MAP(ip1(ip)),x_MAP(ip2(ip)),'rx','MarkerSize',10,'LineWidth',2)
    
    % Adjust figure properties
    set(gca,'TickDir','out')
    xlabel(parlabel(ip1(ip)),'interpreter','latex','fontsize',fontsize);
    if count == 1
        ylabel(parlabel(ip2(ip)),'interpreter','latex','fontsize',fontsize);
        yh1 = get(gca,'ylabel');            % handle to the label object
        py1 = get(yh1,'position');          % get the current position property
        xl = xlim;
        py1(1) = xl(1) - 0.3*(xl(2)-xl(1)); % double the distance, 
                                            % negative values put the label below the axis
        set(yh1,'position',py1);            % set the new position
    else
        set(gca,'YTickLabel',[])
    end    
    set(gca,'TickLabelInterpreter','latex','fontsize',fontsize)
    pos = get(gca,'Position');
    set(gca,'Position',[pos(1) 1.2*pos(2) 1.1*pos(3) pos(4)])
    text(0.02,0.10,letters(ip),'Units','normalized','interpreter','latex','fontsize',fontsize,'fontweight','normal')
    ax = gca;
    ax.TickLength = [0.035 0.035];
    % set box property to off and remove background color
    set(ax,'box','off','color','none')
    drawnow

    count = count + 1;

end