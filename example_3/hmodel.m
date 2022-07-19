function [ Y ] = hmodel(x,plugin)

%% Assign parameters
plugin.data.Imax  = x(1);      % interception storage capacity (mm)
plugin.data.Sumax = x(2);      % unsaturated zone storage capacity (mm)
plugin.data.Qsmax = x(3);      % maximum percolation rate (mm/d)
plugin.data.aE    = x(4);      % evaporation coefficient
plugin.data.aF    = x(5);      % runoff coefficient
plugin.data.aS    = 1e-6;      % percolation coefficient
plugin.data.Kf    = x(6);      % fast-flow response time (d)
plugin.data.Ks    = x(7);      % slow-flow response time (d)

%% Run model C
y = crr_hmodel(plugin.tout,plugin.y0,plugin.data,plugin.options);

% Now compute discharge
Y = diff(y(5,1:end))'; Y = Y(1+plugin.data.wmp:end);