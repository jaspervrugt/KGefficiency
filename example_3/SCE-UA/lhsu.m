function ParGen = lhsu(Xmin,Xmax,s)
% This function performs Latin Hypercube sampling
% 
% ------------------------------------------------------------------------- 
% Input:    Xmin    Minimum possible value of weigths and variances
%           Xmax    Maximum possible value of weigths and variances
%           s       Number of samples 
%
% Output:   ParGen  Samples generated with LHSU sampling
% -------------------------------------------------------------------------

% Calculate the number of variables
nvar = length(Xmin);
% Generate random vector
ran  = rand(s,nvar);
% Initialize ParGen with zeros
ParGen = zeros(s,nvar);
% Do LHSU sampling
for j=1: nvar
   idx = randperm(s);
   P = (idx' - ran(:,j))/s;
   ParGen(:,j) = Xmin(j) + P.* (Xmax(j)-Xmin(j));
end