%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% Vrugt, J.A. and de Oliveira, D.Y., 2022. Confidence intervals of the
%%  Kling-Gupta efficiency, Journal of Hydrology, 612, 
%%  10.1016/j.jhydrol.2022.127968
%% https://doi.org/10.1016/j.jhydrol.2022.127968
%%
%% Section 3.2.2. Application: a hydrologic toy model (Figure 3)
%% 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc, close all

%% Plot parameter pairs obtained using the bootstrap method

% Add folders to path
addpath(genpath([pwd '\results_bootstrap']))

% Load results
load('bootstrap_xopt.mat')
load('bootstrap_logL.mat')

xopt_obs = xopt(1,:);
xopt(1,:) = []; logL(1) = [];

fontsize = 24;

p = 0;
p = p + 1; parlabel{p} = '$D$';
p = p + 1; parlabel{p} = '$K_{\mathrm{f}}$';
p = p + 1; parlabel{p} = '$K_{\mathrm{s}}$';

% 99 95 90 50
color_in = [ 0.3 0.3 0.3; 0.5 0.5 0.5 ; 0.7 0.7 0.7; 0.9 0.9 0.9]; % Try reversed

ip1 = [1 2 3];
ip2 = [2 3 1];

xmin_all = [0.635 0.66 0.034]; xmax_all = [0.665 0.755 0.0375];
ymin_all = [0.66 0.034 0.635]; ymax_all = [0.75 0.0375 0.665];

% Initialize figure
fig = figure('Position',[100 100 1920 2*480]);

count = 1;
letters = {'(A1)';'(B1)';'(C1)'};
for ip=1:3
    ax1 = subplot(2,3,count);
    idx = logL>prctile(logL,1)&logL<=prctile(logL,5);
    scatter(xopt(idx,ip1(ip)),xopt(idx,ip2(ip)),40,'markerfacecolor',color_in(1,:),...
        'markeredgecolor','none');
    hold on

    idx = logL>prctile(logL,5)&logL<=prctile(logL,10);
    scatter(xopt(idx,ip1(ip)),xopt(idx,ip2(ip)),40,'markerfacecolor',color_in(2,:),...
        'markeredgecolor','none');

    idx = logL>prctile(logL,10)&logL<=prctile(logL,50);
    scatter(xopt(idx,ip1(ip)),xopt(idx,ip2(ip)),40,'markerfacecolor',color_in(3,:),...
        'markeredgecolor','none');

    idx = logL>prctile(logL,50);
    scatter(xopt(idx,ip1(ip)),xopt(idx,ip2(ip)),40,'markerfacecolor',color_in(4,:),...
        'markeredgecolor','none');

    plot(xopt_obs(ip1(ip)),xopt_obs(ip2(ip)),'rx','MarkerSize',20,'LineWidth',4)
    xlim([xmin_all(ip1(ip)) xmax_all(ip1(ip))])
    ylim([ymin_all(ip1(ip)) ymax_all(ip1(ip))])

    pos = get(gca,'Position');
    set(gca,'Position',[pos(1) pos(2) pos(3) 1.1*pos(4)])

    xl = xlim;
    yl = ylim;

    c_ls = [xopt_obs(ip1(ip)),xopt_obs(ip2(ip))];
    ls_95(1,:) = prctile(xopt(idx,ip1(ip)),[2.5 97.5]);
    ls_95(2,:) = prctile(xopt(idx,ip2(ip)),[2.5 97.5]);
    factor = 1;

    xmin = xl(1); xmax = xl(2); ymin = yl(1); ymax = yl(2);
    unif_intervals = [16 64 207]/255;
    linewdth_intervals = 0.5; linewdth_intervals2 = 2;
    linewdth_type = '--';
    fontsize_labels = factor*25;
    fontsize_axis = factor*22;
    fontsize_percentiles = factor*20;
    fontsize_ctlabels = factor*16;
    fontsize_title = factor*24;

    % horizontal lines: start
    line(ax1,[xmin xmin + 0.03*(xmax-xmin)],[ls_95(2,1) ls_95(2,1)],'linestyle','-','linewidth',linewdth_intervals2,'color',unif_intervals);
    line(ax1,[xmin xmin + 0.03*(xmax-xmin)],[ls_95(2,2) ls_95(2,2)],'linestyle','-','linewidth',linewdth_intervals2,'color',unif_intervals);
    % horizontal lines: end
    line(ax1,[ls_95(1,1) ls_95(1,2)],[ls_95(2,1) ls_95(2,1)],'linestyle',linewdth_type,'linewidth',linewdth_intervals,'color',unif_intervals);
    line(ax1,[ls_95(1,1) ls_95(1,2)],[ls_95(2,2) ls_95(2,2)],'linestyle',linewdth_type,'linewidth',linewdth_intervals,'color',unif_intervals);
    % vertical lines: start
    line(ax1,[ls_95(1,1) ls_95(1,1)],[ymin ymin + 0.03*(ymax-ymin) ],'linestyle','-','linewidth',linewdth_intervals2,'color',unif_intervals);
    line(ax1,[ls_95(1,2) ls_95(1,2)],[ymin ymin + 0.03*(ymax-ymin) ],'linestyle','-','linewidth',linewdth_intervals2,'color',unif_intervals);
    % vertical lines: end
    line(ax1,[ls_95(1,1) ls_95(1,1)],[ls_95(2,1) ls_95(2,2)],'linestyle',linewdth_type,'linewidth',linewdth_intervals,'color',unif_intervals);
    line(ax1,[ls_95(1,2) ls_95(1,2)],[ls_95(2,1) ls_95(2,2)],'linestyle',linewdth_type,'linewidth',linewdth_intervals,'color',unif_intervals);
    % optimum solution
    line(ax1,[xmin c_ls(1)],[c_ls(2) c_ls(2)],'linestyle','--','linewidth',1,'color','r');
    line(ax1,[c_ls(1) c_ls(1)],[ymin c_ls(2)],'linestyle','--','linewidth',1,'color','r');

    set(gca,'TickDir','out')
    set(gca,'XTickLabel',[])
    ylabel(parlabel(ip2(ip)),'interpreter','latex','fontsize',fontsize);
    if count == 1
        set(gca,'YTickLabel',{'0.66';'0.68';'0.70';'0.72';'0.74'})
    else
    end
    xlim([xmin xmax])
    ylim([ymin ymax])

    set(gca,'TickLabelInterpreter','latex','fontsize',fontsize)
    text(0.02,0.92,letters(ip),'Units','normalized','interpreter','latex','fontsize',fontsize,'fontweight','normal')
    ax = gca;
    ax.TickLength = [0.035 0.035];
    % set box property to off and remove background color
    set(ax,'box','off','color','none')
    drawnow
    count = count + 1;

end

%% Plot parameter pairs obtained using the DREAM algorithm

% Add folders to path
addpath(genpath([pwd '\results_DREAM']))

% Load results
load('DREAM_x.mat')
load('DREAM_xpost.mat')
load('DREAM_logL.mat')

[~,iMAP] = max(x(:,end));
x = x(:,1:end-2);
xopt_obs = x(iMAP,:);
x_DREAM = xpost;

count = 1;
letters = {'(A2)';'(B2)';'(C2)'};
for ip=1:3
    ax1 = subplot(2,3,count+3);
    idx = logL>prctile(logL,1)&logL<=prctile(logL,5);
    scatter(x_DREAM(idx,ip1(ip)),x_DREAM(idx,ip2(ip)),40,'markerfacecolor',color_in(1,:),...
        'markeredgecolor','none');
    hold on

    idx = logL>prctile(logL,5)&logL<=prctile(logL,10);
    scatter(x_DREAM(idx,ip1(ip)),x_DREAM(idx,ip2(ip)),40,'markerfacecolor',color_in(2,:),...
        'markeredgecolor','none');

    idx = logL>prctile(logL,10)&logL<=prctile(logL,50);
    scatter(x_DREAM(idx,ip1(ip)),x_DREAM(idx,ip2(ip)),40,'markerfacecolor',color_in(3,:),...
        'markeredgecolor','none');

    idx = logL>prctile(logL,50);
    scatter(x_DREAM(idx,ip1(ip)),x_DREAM(idx,ip2(ip)),40,'markerfacecolor',color_in(4,:),...
        'markeredgecolor','none');

    plot(xopt_obs(ip1(ip)),xopt_obs(ip2(ip)),'rx','MarkerSize',20,'LineWidth',4)
    xlim([xmin_all(ip1(ip)) xmax_all(ip1(ip))])
    ylim([ymin_all(ip1(ip)) ymax_all(ip1(ip))])

    pos = get(gca,'Position');
    set(gca,'Position',[pos(1) 1.2*pos(2) pos(3) 1.1*pos(4)])

    xl = xlim;
    yl = ylim;

    c_ls = [xopt_obs(ip1(ip)),xopt_obs(ip2(ip))];
    ls_95(1,:) = prctile(x_DREAM(idx,ip1(ip)),[2.5 97.5]);
    ls_95(2,:) = prctile(x_DREAM(idx,ip2(ip)),[2.5 97.5]);
    factor = 1;

    xmin = xl(1); xmax = xl(2); ymin = yl(1); ymax = yl(2);
    unif_intervals = [16 64 207]/255;
    linewdth_intervals = 0.5; linewdth_intervals2 = 2;
    linewdth_type = '--';
    fontsize_labels = factor*25;
    fontsize_axis = factor*22;
    fontsize_percentiles = factor*20;
    fontsize_ctlabels = factor*16;
    fontsize_title = factor*24;

    % horizontal lines: start
    line(ax1,[xmin xmin + 0.03*(xmax-xmin)],[ls_95(2,1) ls_95(2,1)],'linestyle','-','linewidth',linewdth_intervals2,'color',unif_intervals);
    line(ax1,[xmin xmin + 0.03*(xmax-xmin)],[ls_95(2,2) ls_95(2,2)],'linestyle','-','linewidth',linewdth_intervals2,'color',unif_intervals);
    % horizontal lines: end
    line(ax1,[ls_95(1,1) ls_95(1,2)],[ls_95(2,1) ls_95(2,1)],'linestyle',linewdth_type,'linewidth',linewdth_intervals,'color',unif_intervals);
    line(ax1,[ls_95(1,1) ls_95(1,2)],[ls_95(2,2) ls_95(2,2)],'linestyle',linewdth_type,'linewidth',linewdth_intervals,'color',unif_intervals);
    % vertical lines: start
    line(ax1,[ls_95(1,1) ls_95(1,1)],[ymin ymin + 0.03*(ymax-ymin) ],'linestyle','-','linewidth',linewdth_intervals2,'color',unif_intervals);
    line(ax1,[ls_95(1,2) ls_95(1,2)],[ymin ymin + 0.03*(ymax-ymin) ],'linestyle','-','linewidth',linewdth_intervals2,'color',unif_intervals);
    % vertical lines: end
    line(ax1,[ls_95(1,1) ls_95(1,1)],[ls_95(2,1) ls_95(2,2)],'linestyle',linewdth_type,'linewidth',linewdth_intervals,'color',unif_intervals);
    line(ax1,[ls_95(1,2) ls_95(1,2)],[ls_95(2,1) ls_95(2,2)],'linestyle',linewdth_type,'linewidth',linewdth_intervals,'color',unif_intervals);
    % optimum solution
    line(ax1,[xmin c_ls(1)],[c_ls(2) c_ls(2)],'linestyle','--','linewidth',1,'color','r');
    line(ax1,[c_ls(1) c_ls(1)],[ymin c_ls(2)],'linestyle','--','linewidth',1,'color','r');

    set(gca,'TickDir','out')
    xlabel(parlabel(ip1(ip)),'interpreter','latex','fontsize',fontsize);
    ylabel(parlabel(ip2(ip)),'interpreter','latex','fontsize',fontsize);
    if count == 1
        set(gca,'YTickLabel',{'0.66';'0.68';'0.70';'0.72';'0.74'})
    else
    end
    xlim([xmin xmax])
    ylim([ymin ymax])

    if count == 2
        set(gca,'XTickLabel',{'0.66';'0.68';'0.70';'0.72';'0.74'})
    end

    set(gca,'TickLabelInterpreter','latex','fontsize',fontsize)
    text(0.02,0.92,letters(ip),'Units','normalized','interpreter','latex','fontsize',fontsize,'fontweight','normal')
    ax = gca;
    ax.TickLength = [0.035 0.035];
    % set box property to off and remove background color
    set(ax,'box','off','color','none')
    drawnow
    count = count + 1;
end