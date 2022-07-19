function [Q,S_Dt] = LinRes_AS(S0,P,K,Dt)
% Analytical solution of the linear reservoir
S_Dt = P/K + (S0-P/K)*exp(-K*Dt); % Calculate the new storage at time t+Dt
Q = P - (S_Dt-S0)/Dt; % Calculate the average discharge over the time step Dt
end

