function [c_opt,KG] = getoptKGE(f,t,yr)
% Returns optimum parameters of f given replicate yr 
% maximizes KG statistic

sr = std(yr); mr = mean(yr);

c0 = [2 1];

[c_opt,KGE] = fminsearch(@(c) KGE_func(c,f,t,yr,sr,mr),c0); KG = -KGE; 

function KG = KGE_func(c,f,t,yr,sr,mr)

% evaluate model at x
y = f(t,c); 
ss = std(y); ms = mean(y); 
alfa = ss/sr; gamma = ms/mr;
r = corrcoef(y,yr); r = r(1,2);

% Compute KG efficiency
KGE = 1 - sqrt((r-1)^2 + (alfa-1)^2 + (gamma-1)^2);
% Nelder-Mead wants to minimize
KG = -KGE; 


