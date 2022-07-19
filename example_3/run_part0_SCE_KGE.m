%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% Vrugt, J.A. and de Oliveira, D.Y., 2022. Confidence intervals of the
%%  Kling-Gupta efficiency, Journal of Hydrology, 612, 
%%  10.1016/j.jhydrol.2022.127968
%% https://doi.org/10.1016/j.jhydrol.2022.127968
%%
%% Section 3.3.2. Application: a conceptual watershed model 
%%  This code will generate results stored in '\results_obs'
%% 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc, close all
rng('default')

% Add folders to path
addpath(genpath([pwd '\data']))
addpath(genpath([pwd '\model']))
addpath(genpath([pwd '\SCE-UA']))

% iq = 187; % Leaf Riv er near Collins, MS (USGS 02472000)
iq = 165; %  Kinchafoonee Creek near Dawson, GA (USGS 02350900)

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
ParRange.minn(1) = 0.5;     ParRange.maxn(1) = 10;
ParRange.minn(2) = 10;      ParRange.maxn(2) = 1000;
ParRange.minn(3) = 0;       ParRange.maxn(3) = 1000;
ParRange.minn(4) = 1e-6;    ParRange.maxn(4) = 100;
ParRange.minn(5) = -10;     ParRange.maxn(5) = 10;
ParRange.minn(6) = 0;       ParRange.maxn(6) = 10;
ParRange.minn(7) = 0;       ParRange.maxn(7) = 500;

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
SCEPar.ndraw = 20000;        % Maximum number of function evaluations
SCEPar.p = 10;              % Number of complexes
SCEPar.alpha = 1;	 		% Number of Simplexes
SCEPar.Gamma = 0;			% Kurtosis parameter Bayesian Inference Scheme

% Define the option - KGE
option = 5;

Extra.calPeriod = [1 length(plugin.data.P)-plugin.data.wmp];
Extra.MaxT = length(plugin.data.P);
Extra.plotYN = false;

%%
Measurement = Y;
for ii=1
       
    % Run shuffled complex evolution optimization algorithm
    [D,x,ParSet] = sce_ua(SCEPar,Func_name,ParRange,Measurement,Extra,option,plugin);
    
    rowSelection = ParSet(:,SCEPar.n+1)==max(ParSet(:,SCEPar.n+1));
    uniqueSets = unique(ParSet(rowSelection,1:SCEPar.n),'rows');
    if size(uniqueSets,1)>1
        warning('Multiple best sets. Taking the first.')
    end
    
    xopt = uniqueSets(1,:);   
    Yopt = hmodel(xopt,plugin);
    Yobs = Measurement;
    
    KG = parKGE(Yopt,Yobs);
    
end
save(['SCE_obs_xopt_' num2str(iq) '.mat'],'xopt')
save(['SCE_obs_Yopt_' num2str(iq) '.mat'],'Yopt')
save(['SCE_obs_KG_' num2str(iq) '.mat'],'KG')