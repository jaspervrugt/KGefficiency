function [Wb,Cb] = calcCbWb(Beta);
% This function calculates the exponential power density
% Equation [19 and 20] paper by Thiemann et al
A1 = gamma(3*(1+Beta)/2); A2 = gamma((1+Beta)/2); 
Cb = (A1/A2)^(1/(1+Beta));
Wb = sqrt(A1)/((1+Beta)*(A2^(1.5)));


