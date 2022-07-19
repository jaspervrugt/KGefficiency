function [Ir,exitflag] = secant(Ii,maxiter,tolfun,b,beta,s2_eps,Q,Z,alfa,p,nu,method)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The Secant method for root finding of infiltration equation of Haverkamp           %%
%%  SYNOPSIS: Ir = secant(fun,Ii);                                                    %%
%%            [Ir,exitflag] = secant(fun,Ii);                                         %%
%%            [Ir,exitflag] = secant(fun,Ii,maxiter);                                 %%
%%            [Ir,exitflag] = secant(fun,Ii,maxiter,tolfun);                          %%
%% where                                                                              %%
%%       fun:      [input] anonymous function handle of right-hand-side Haverkamp     %%
%%       Ii:       [input] starting point (cumulative infiltration in cm)             %%
%%       maxiter:  [optional input] maximum number of iterations (default: 20)        %%
%%       tolfun:   [optional input] tolerance on root function value (default: 1e-10) %%
%%       Ir:       [output] cumulative infiltration (in cm) at specified time         %%
%%       exitflag: [output] exit condition. Possible conditions are:                  %%
%%                          1: secant found a zero point ( = root)                    %%
%%                          2: secant terminated with infinite root                   %%
%%                          3: secant terminated with infinite function value root    %%
%%                          4: secant terminated with imaginary function value root   %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization part: Define settings of Secant method
if nargin < 12, s2_eps = []; end
if nargin < 12, Z = []; end
if nargin < 11, Q = []; end
if nargin < 4, tolfun = 1e-12; end
if nargin < 3, maxiter = 20; end
exitflag = 1; n = 3;                                          % Exitflag and iteration
switch method 
    case 1 % Draper and Guttman
        % https://www.jstor.org/stable/pdf/2348711.pdf?refreqid=excelsior%3A56d9087072e4bd5e64c55efdcb39afbf
        delta_F = s2_eps * p * finv(alfa,p,nu); 
        % Note for homoscedastic errors matrix Q will equal matrix D
        % method 3 and method 1 are identical --> method 1 for all type of errors
        fun = @(b2) (beta - [b;b2])' * (Q'*Q) * (beta - [b;b2]) - delta_F;
    case 2 % Press et al (1992) --> uncorrelated and homoscedastic errors only
       % delta_WSSR = chi2inv(alfa,p);
        delta_F = p * finv(alfa,p,nu); 
        % Q --> from function call in main program, Q = inv(C) (thus includes sigma2 multiplier in cov matrix)
        fun = @(b2) (beta - [b;b2])' * Q * (beta - [b;b2]) - delta_F; %delta_WSSR;
    case 3 % Test - heteroscedastic/correlated error
        var_res = (Z'*Z - beta'*Q'*Z)/nu; % --> variance of residuals (= s2_eps)
        delta_F = var_res * p * finv(alfa,p,nu); % allowable increment
        fun = @(b2) (beta - [b;b2])' * (Q'*Q) * (beta - [b;b2]) - delta_F; 
        % Equivalent to case 1 above
    case 4 % JAV method: as used in paper
        Vi = Q; D = Z;
        delta_F = s2_eps * p * finv(alfa,p,nu);
        fun = @(b2) (beta - [b;b2])' * (D'*Vi*D) * (beta - [b;b2]) - delta_F;        
end

I(n-2) = Ii; I(n-1) = Ii + 0.2;                               % First pair of triple
fun(I(2));
I(n) = (I(1)*fun(I(2))-I(2)*fun(I(1))) ... 
    /(fun(I(2))-fun(I(1))); % Last point of triple

%% Dynamic part: Iteratively refine estimate of root 
while (abs(fun(I(n))) > tolfun) && ((n-2) < maxiter)        % Secant method
    n = n + 1;                                                % Increment iteration
    I(n) = (I(n-2)*fun(I(n-1)) - I(n-1)*fun(I(n-2))) ...        % Update guess of root
        / (fun(I(n-1)) - fun(I(n-2)));
end 
%% End of Dynamic Part

Ir = I(n); fIr = fun(Ir);                                     % root guess + fun(root)   
if abs(fIr) > tolfun 
    exitflag = 5;                                             % no root
elseif ~isfinite(I)                                           % determine exitflag
    exitflag = 2;
elseif ~isfinite(fIr)       
    exitflag = 3; 
elseif ~isreal(fIr)
    exitflag = 4;
end