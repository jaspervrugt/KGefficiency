function [KGE,alpha,beta,r] = parKGE(sim,obs)

sim = sim(:);
obs = obs(:);

% Calculate mean sim and obs
mean_sim = mean(sim);
mean_obs = mean(obs);

% Calculate alpha component
alpha = std(sim) / std(obs);

% Calculate beta component
beta = mean_sim / mean_obs;

% Calculate r component
r = corr(sim, obs, 'type', 'pearson');

% Return Non-Parametric Efficiency value
KGE = (1 - sqrt((alpha - 1)^2 + (beta - 1)^2 + (r - 1)^2));