%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% Vrugt, J.A. and de Oliveira, D.Y., 2022. Confidence intervals of the
%%  Kling-Gupta efficiency, Journal of Hydrology, 612, 
%%  10.1016/j.jhydrol.2022.127968
%% https://doi.org/10.1016/j.jhydrol.2022.127968
%%
%% Section 3.3.2. Application: a conceptual watershed model 
%%  This code will generate results stored in '\results_bootstrap'
%% 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc, close all
rng('default')

% Add folders to path
addpath(genpath([pwd '\data']))
addpath(genpath([pwd '\model']))
addpath(genpath([pwd '\SCE-UA']))

iq = 187; % Leaf River near Collins, MS (USGS 02472000)
% iq = 165; % Kinchafoonee Creek near Dawson, GA (USGS 02350900)

% Load watershed area
load(['data_' num2str(iq) '.mat'])

%% Define name of function (.m file) for posterior exploration
Func_name = 'hmodel';

%% Define the structure plugin -- as second argument to model function
plugin.data.wmp = 5*365;
% First five years are warm-up
plugin.data.idx = [1+plugin.data.wmp:size(data,1)]';

%% Then read the boundary conditions
plugin.data.P = data(:,1);
plugin.data.Ep = data(:,2);

%% Initial conditions
plugin.y0 = 1e-6.*ones(5,1);

%% Discharge time series
Y = data(1+plugin.data.wmp:end,3);

%%
% Give the parameter ranges (minimum and maximum values)
if iq==165
    ParRange.minn(1) = 8;    ParRange.maxn(1) = 10;
    ParRange.minn(2) = 400;  ParRange.maxn(2) = 600;
    ParRange.minn(3) = 0.5;  ParRange.maxn(3) = 1.2;
    ParRange.minn(4) = 2;    ParRange.maxn(4) = 4;
    ParRange.minn(5) = -2.5; ParRange.maxn(5) = -0.5;
    ParRange.minn(6) = 5;    ParRange.maxn(6) = 7;
    ParRange.minn(7) = 0;    ParRange.maxn(7) = 100;
elseif iq==187
    ParRange.minn(1) = 9;   ParRange.maxn(1) = 10;
    ParRange.minn(2) = 200; ParRange.maxn(2) = 300;
    ParRange.minn(3) = 0;   ParRange.maxn(3) = 0.0001;
    ParRange.minn(4) = 0;   ParRange.maxn(4) = 0.01;
    ParRange.minn(5) = -6;  ParRange.maxn(5) = -5;
    ParRange.minn(6) = 3;   ParRange.maxn(6) = 4;
    ParRange.minn(7) = 0;   ParRange.maxn(7) = 500;
end

%% Integration options
plugin.options.InitialStep = 1;                 % initial time-step (d)
plugin.options.MaxStep     = 1;                 % maximum time-step (d)
plugin.options.MinStep     = 1e-5;              % minimum time-step (d)
plugin.options.RelTol      = 1e-4;              % relative tolerance
plugin.options.AbsTol      = 1e-4*ones(5,1); % absolute tolerances (mm)
plugin.options.Order       = 2;                 % 2nd order accurate method (Heun)
%% Running time
plugin.tout = [ 0 : size(plugin.data.P,1) ];

% Generate mex file (need to do this only once)
% mex crr_hmodel.c

%% Define SCE parameters
SCEPar.n = 7;               % Dimension of the problem (Nr. parameters to be optimized in the model)
SCEPar.ndraw = 1500;        % Maximum number of function evaluations
SCEPar.p = 10;              % Number of complexes
SCEPar.alpha = 1;	 		% Number of Simplexes
SCEPar.Gamma = 0;			% Kurtosis parameter Bayesian Inference Scheme

% Define the option - KGE
option = 5;

Extra.calPeriod = [1 length(plugin.data.P)-plugin.data.wmp];
Extra.MaxT = length(plugin.data.P);
Extra.plotYN = false;

%% Load replicates
load(['Y_err_' num2str(iq) '.mat'])
N = 1000; % Number of replicates
     
%%
xopt = nan(N+1,SCEPar.n);
Yopt = nan(length(Y),N+1);
KG = nan(N+1,1);

Measurement = Y; Yobs = Y;
for ii=1
       
    % Run shuffled complex evolution optimization algorithm
    [D,x,ParSet] = sce_ua(SCEPar,Func_name,ParRange,Measurement,Extra,option,plugin);
    
    rowSelection = ParSet(:,SCEPar.n+1)==max(ParSet(:,SCEPar.n+1));
    uniqueSets = unique(ParSet(rowSelection,1:SCEPar.n),'rows');
    if size(uniqueSets,1)>1
        warning('Multiple best sets. Taking the first.')
    end

    xopt(ii,:) = uniqueSets(1,:);
    Yopt(:,ii) = hmodel(xopt(ii,:),plugin);
    KG(ii,1) = parKGE(Yopt(:,ii),Yobs);
    
end
% profile report
save(['SCE_xopt_' num2str(iq) '.mat'],'xopt')
save(['SCE_Yopt_' num2str(iq) '.mat'],'Yopt')
save(['SCE_KG_' num2str(iq) '.mat'],'KG')

%%
parfor ii=1:N  
    
    Measurement = Y_err(:,ii);

    % Run shuffled complex evolution optimization algorithm
    [D,x,ParSet] = sce_ua(SCEPar,Func_name,ParRange,Measurement,Extra,option,plugin);
    
    rowSelection = ParSet(:,SCEPar.n+1)==max(ParSet(:,SCEPar.n+1));
    uniqueSets = unique(ParSet(rowSelection,1:SCEPar.n),'rows');
    if size(uniqueSets,1)>1
        warning('Multiple best sets. Taking the first.')
    end
    
    xopt(ii+1,:) = uniqueSets(1,:);
    Yopt(:,ii+1) = hmodel(xopt(ii+1,:),plugin);
    KG(ii+1,1) = parKGE(Yopt(:,ii+1),Yobs);

    disp([num2str(ii) ' of ' num2str(N)])
    
end
save(['SCE_xopt_' num2str(iq) '.mat'],'xopt')
save(['SCE_Yopt_' num2str(iq) '.mat'],'Yopt')
save(['SCE_KG_' num2str(iq) '.mat'],'KG')