%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% Vrugt, J.A. and de Oliveira, D.Y., 2022. Confidence intervals of the
%%  Kling-Gupta efficiency, Journal of Hydrology, 612, 
%%  10.1016/j.jhydrol.2022.127968
%% https://doi.org/10.1016/j.jhydrol.2022.127968
%%
%% Section 3.2.2. Application: a hydrologic toy model 
%%  This code will generate the file 'bootstrap_xopt.mat' stored in
%%  '\results_obs'
%% 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, close all, clc
rng('default')

% Add folders to path
addpath(genpath([pwd '\data']))
addpath(genpath([pwd '\model']))
addpath(genpath([pwd '\SCE-UA']))

% Load data
load('data.mat','data')

%% Define the structure plugin -- as second argument to model function
plugin.data.wmp = 1*365; % First year is warm-up
plugin.data.idx = [1+plugin.data.wmp:size(data,1)]';

%% Then read the boundary conditions
plugin.data.Pt = data(:,1);
plugin.data.Ep = data(:,2);

%% Initial condition
plugin.y0 = 50; 

%% Run model
tmax = length(plugin.data.Pt);
t = 1:tmax-365;
x = [0.65 0.7 0.035]; % model parameters
[y,yf,ys] = HydroModel(x,plugin);

%% Create pseudo-observed discharge record (with heteroscedastic measurement error)
sigma = 0.10*y;
y_obs = y + sigma.*randn(length(y),1);
hold on
plot(t,y_obs,'r.')
drawnow

%% Create replicates
N = 1000;
y_err = nan(length(y),N);
for ii=1:N
    y_err(:,ii) = y_obs + sigma.*randn(length(y),1);
end

%% Define SCE parameters
SCEPar.n = 3;               % Dimension of the problem (Nr. parameters to be optimized in the model)
SCEPar.ndraw = 500;         % Maximum number of function evaluations
SCEPar.p = 5;               % Number of complexes
SCEPar.alpha = 1;	 		% Number of Simplexes
SCEPar.Gamma = 0;			% Kurtosis parameter Bayesian Inference Scheme

%% Define the option - GLS
option = 6;

Extra.calPeriod = [1 length(plugin.data.Pt)-plugin.data.wmp];
Extra.MaxT = length(plugin.data.Pt);
Extra.plotYN = false;

% Give the parameter ranges (minimum and maximum values)
ParRange.minn = [ 0.0 0.0 0.00];
ParRange.maxn = [ 1.0 1.0 0.15];

Func_name = 'HydroModel';

xopt = nan(N+1,3);
Measurement = y_obs;

C = diag(sigma.^2);
invC = inv(C);
Extra.invC = invC;

for ii=1
    
    % Run shuffled complex evolution optimization algorithm
    [D,x,ParSet] = sce_ua(SCEPar,Func_name,ParRange,Measurement,Extra,option,plugin);
    
    rowSelection = ParSet(:,SCEPar.n+1)==max(ParSet(:,SCEPar.n+1));
    uniqueSets = unique(ParSet(rowSelection,1:SCEPar.n),'rows');
    if size(uniqueSets,1)>1
        warning('Multiple best sets. Taking the first.')
    end
    
    % Save parameter set
    xopt(ii,:) = uniqueSets(1,:);
    
end
save('bootstrap_xopt.mat','xopt')

Measurement = y_obs;
for ii=1:N
    
    Measurement = y_err(:,ii);
    
    % Run shuffled complex evolution optimization algorithm
    [D,x,ParSet] = sce_ua(SCEPar,Func_name,ParRange,Measurement,Extra,option,plugin);
    
    rowSelection = ParSet(:,SCEPar.n+1)==max(ParSet(:,SCEPar.n+1));
    uniqueSets = unique(ParSet(rowSelection,1:SCEPar.n),'rows');
    if size(uniqueSets,1)>1
        warning('Multiple best sets. Taking the first.')
    end
    
    % Save parameter set
    xopt(ii+1,:) = uniqueSets(1,:);

    % Display progress
    disp([num2str(ii) ' of ' num2str(N)])
    
end
save('bootstrap_xopt.mat','xopt')