%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% Vrugt, J.A. and de Oliveira, D.Y., 2022. Confidence intervals of the
%%  Kling-Gupta efficiency, Journal of Hydrology, 612, 
%%  10.1016/j.jhydrol.2022.127968
%% https://doi.org/10.1016/j.jhydrol.2022.127968
%%
%% Section 3.2.2. Application: a hydrologic toy model 
%%  This code will generate the file 'bootstrap_logL.mat' stored in
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

%% Compute log-likehood

C = diag(sigma.^2);
invC = inv(C);

load('bootstrap_xopt.mat')

logL = nan(size(xopt,1),1);
for ii=1:N+1
    y = HydroModel(xopt(ii,:),plugin);
    res = y_obs - y;
    logL(ii) = - 1/2 * sum(res.*(invC*res));
end
save('bootstrap_logL.mat','logL')