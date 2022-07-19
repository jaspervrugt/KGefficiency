function [D] = sortje(OF);
% Sorts the s points in order of increasing function value
D = -sortrows(-OF,[1]);