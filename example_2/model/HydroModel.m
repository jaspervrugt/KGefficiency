function [Q,Qf,Qs] = HydroModel(x,plugin)

% Forcing data
P = plugin.data.Pt;

% Initial state
S0 = plugin.y0;

% Model parameters
D = x(1);
Kf = x(2);
Ks = x(3);

Dt = 1;
tmax = length(P);

% Initialize reservoirs
Qf = nan(tmax,1); S = nan(tmax,1); S(1) = S0; 

% Run analitical solution
for t=1:tmax
    [Qf(t),S(t)] = LinRes(S(t),(1-D)*P(t),Kf,Dt);
    if t<tmax
        S(t+1) = S(t);
    end
end

% iInitialize reservoirs
Qs = nan(tmax,1); S = nan(tmax,1); S(1) = S0;  

% Run analitical solution
for t=1:tmax
    [Qs(t),S(t)] = LinRes(S(t),D*P(t),Ks,Dt);
    if t<tmax
        S(t+1) = S(t);
    end
end
Qf = Qf(1+plugin.data.wmp:end);
Qs = Qs(1+plugin.data.wmp:end);
Q = Qf + Qs;