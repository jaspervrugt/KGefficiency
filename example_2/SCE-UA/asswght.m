function [rho] = asswght(m);
% Assigns the weights to the individual points in the complex
for i=m:-1:1,
   rho(i,1) = ((2*(m+1-i))/(m*(m+1)));
end;

