%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% Vrugt, J.A. and de Oliveira, D.Y., 2022. Confidence intervals of the
%%  Kling-Gupta efficiency, Journal of Hydrology, 612, 
%%  10.1016/j.jhydrol.2022.127968
%% https://doi.org/10.1016/j.jhydrol.2022.127968
%%
%% Test different methods for estimation of confidence region/intervals of 
%% regression coefficients: linear regression function using methods that
%% built on foundation of least squares
%% 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% close all active figures, clear memory and screen
close all hidden; clear; clc

%% Synthetic case study: Create measurement data using one of four methods
approach = 1;   % [1] homoscedastic error only
                % [2] heteroscedastic error only
                % [3] homoscedastic + correlated error
                % [4] heteroscedastic + correlated error
load_data = 2;  % [1] remake data 
                % [2] load data of trial (of paper)
fig1 = 1;       % [0] do not create Figure 1 / [1] create figure 1 
fig11 = 1;      % [0] do not create Figure 11 / [1] create figure 11

%% Define function and basic settings
f = @(t,c) c(1)*t + c(2);       % function
c_dgp = [1 2]';                 % coefficients of dgp
alfa = [0.01 0.05 0.10 0.50];   % significance levels
N = 1e5;                        % number of bootstrap samples
M = 200;                        % number of samples for gridding WSSR

%% Define model output, data, measurement error cov, D, V, P, Q, etc.
switch load_data
    case 1
        t = (0:1:50)';                  % independent variable
        p = numel(c_dgp);               % number of parameters
        y = f(t,c_dgp);                 % dgp: true
        n = numel(t);                   % number of data points
        D = [ t  ones(n,1) ];           % define design matrix
        % define measurement error
        switch approach
            case 1
                sigma2_eps = 0.5; V = eye(n); 
            case 2
                a = 0.1; b = 0.01;
                sigma2_eps_all = (a*abs(y) + b).^2;
                sigma2_eps = min(sigma2_eps_all);
                V = diag(sigma2_eps_all/sigma2_eps);
                clear sigma2_eps_all; 
            otherwise
                phi = 0.8; R = nan(n,n);
                % Now create first order correlation matrix 
                for i = 1:n
                    R(i,i) = 1; idl = i-1:-1:1; idr = i+1:n;  
                    R(i,1:i-1) = phi.^idl;              % lower triangle
                    R(i,idr) = phi.^(1:numel(idr));     % upper triangle
                end
                V = R;
                switch approach
                    case 3 % V matrix is done
                        sigma2_eps = 0.5;
                    case 4 % Update diagonal of V for heteroscedasticity
                        a = 0.1; b = 0.01;
                        sigma2_eps_all = (a*abs(y) + b).^2;
                        sigma2_eps = min(sigma2_eps_all);
                        for i = 1:n, V(i,i) = sigma2_eps_all(i)/sigma2_eps; end
                        clear sigma2_eps_all; 
                end
        end
        Sigma2_eps = sigma2_eps * V;    % Measurement error covariance matrix
        Vi = inv(V);                    % Get inverse of matrix V 
       % P = chol(V);                   % Now get P'*P=V or P*P=V or P^2 = V
        meas_err = mvnrnd(zeros(n,1),...    
            Sigma2_eps,1)';             % Realization of measurement errors
        y_meas = y + meas_err;          % Measurements of dgp
       % Z = P\y_meas; Q = P\D;         % Get nx1 vector Z and nxp matrix Q
       % Save data for DREAM - to compare results with bootstrap
       % evalstr = strcat('save',{' '},'data_',num2str(approach),...
       % {' '},'t n p meas_err D c_dgp Sigma2_eps sigma2_eps V y_meas'); 
       % eval(char(evalstr));
    case 2 % Load data
        evalstr = strcat('load',{' '},'data_',num2str(approach)); 
        eval(char(evalstr)); Vi = inv(V);
end

% Compute confidence levels
gam = 1 - alfa; 


%% Different procedures for constructing confidence regions and intervals

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1: Least squares
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_ls = (D'*Vi*D)\(D'*Vi*y_meas);    % Least squares solution - sigma2_eps_max multiplier cancels
res = y_meas - f(t,c_ls);           % residual
s2_eps = res'*Vi*res/(n-p);         % sample variance of residual (should be sigma2_eps for n --> infty)
S2_eps = s2_eps * V;                % sample covariance matrix measurement errors
C = s2_eps * inv(D'*Vi*D);          % variance-covariance matrix of parameters
   % = inv(D'*inv(S2_eps)*D)
   % = s2_eps*inv(Q'*Q)
dc = tinv(0.975,n-p)*sqrt(diag(C)); % coefficient ranges at 95%
ls_95 = [ c_ls-dc , c_ls+dc ];      % 95% parameter ranges
SSR_opt = res'*inv(S2_eps)*res;     % weigted SSR with GLS parameter values
    % = res'/S2_eps * res;
% We could also get the decorrelated normalized residuals via the weight matrix
% W = chol(inv(Sigma2_eps));        % W'*W = inv(Sigma2_eps); so that W has units of residuals
W = chol(inv(S2_eps));              % W'*W = inv(S2_eps);
vareps = (res'*Vi*res)/s2_eps;      % shows that you can extract sigma2
                                    % denominator normalizes vareps -->
                                    % var(vareps) = 1; 
vareps = W*res;                     % Thus, vareps is equal to line above
% Now look at original residuals
rho2 = autocorr(res,10);            % Indeed, if approach 3 or 4
% Indeed rho(1) = 0.8 = phi1; rho(2) = phi^2, etc.
rho1 = autocorr(vareps,10);         % Indeed, no correlation
% -----------------------------

%%% Additional notes
% 1. homoscedastic/heteroscedastic errors only
% sum((res./sqrt(diag(S2_eps))).^2) % --> via sum operator (weighted residuals)
% (W*res)'*(W*res)                  % --> vector form

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2: Bootstrap method
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bopt = nan(p+1,N); resb = nan(n,N);               
R = y_meas + mvnrnd(zeros(n,1),S2_eps,N)';              % Create ensemble of N measured records
for i = 1:N
    Bopt(1:p,i) = (D'*Vi*D)\(D'*Vi*R(1:n,i));           % Optimum
    resb(1:n,i) = y_meas - f(t,Bopt(1:p,i));            % residual - store for later
    %Bopt(p+1,i) = resb'/S2_eps * resb;                 % weigted SSR with GLS parameter values
    Bopt(p+1,i) = resb(1:n,i)'*Vi*resb(1:n,i)/s2_eps;   % --> much faster
end
bs_95 = nan(p,2);                   % 95% pertaining to sigma^2 of data
id_lu = round([0.025 0.975]*N);     % 2.5 and 97.5 percentile
for j = 1:p
    a = sort(Bopt(j,1:N)); bs_95(j,1:2) = [a(id_lu(1)) a(id_lu(2))];
end
% -----------------------------
% Postprocess bootstrap results for plotting
Bopt = sortrows(Bopt',p+1);         % sort in order of ascending SSR
id_up = round(gam*N);               % percentiles
% Now find which ones satisfy the threshold
bs_50id = 1:id_up(4);
bs_90id = id_up(4)+1:id_up(3);      % no double plotting of points
bs_95id = id_up(3)+1:id_up(2);      % no double plotting of points
bs_99id = id_up(2)+1:id_up(1);      % no double plotting of points

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3: Confidence region using F-distribution: for p <= 3
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_gam = numel(gam);                 % number of gamma values
beta = c_ls;                        % set beta to least squares solution
b_low = c_ls - 2*dc;                % lower values
b_up = c_ls + 2*dc;                 % upper values
for j = 1:p-1
    b(j,1:M) = [ b_low(j):(b_up(j)-b_low(j))/(M-1):b_up(j) ];
end
% method = 1: F-distribution
method = 1; maxiter = 1000; tolfun = 1e-10; % Define maxiter and tolfun for secant method
% Solve for ellipsoid for given alfa's (critical points)
% P = chol(V); Z = P\y_meas; Q = P\D;       % Get nx1 vector Z and nxp matrix Q
% [semimajor_axisB,semiminor_axisB,X0B,Y0B,phiB] = ...
% get_ellipse(gam,beta,dc,p,n,M,b,s2_eps,Q,Z,method,maxiter,tolfun);
% --> from this output we can plot the ellipsoidal regions

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4: Confidence region using chi-square distribution: for p <= 3
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% method 2: chi-square distribution and SSmin
method = 2; maxiter = 1000; tolfun = 1e-10;     % Define maxiter and tolfun for secant method
% Solve for ellipsoid for given gamma values (critical points)
P = chol(V); Z = P\y_meas; Q = P\D;             % Get nx1 vector Z and nxp matrix Q
[semimajor_axisC,semiminor_axisC,X0C,Y0C,phiC] = ...
    get_ellipse(gam,beta,dc,p,n,M,b,s2_eps,inv(C),Z,method,maxiter,tolfun);
% We submit for Q the covariance matrix, C, of the parameters
% --> from this output we can plot the ellipsoidal regions

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5: Confidence region from gridding of GLS: for p = 2
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:2
    b(j,1:M) = beta(j)-5*dc(j) : 10*dc(j)/(M-1) : beta(j)+5*dc(j);
end
[b1,b2] = meshgrid(b(1,1:M),b(2,1:M)); SSR = nan(M,M);
counter = 1; G = nan(M^2,p+1);
for i = 1:M
    for j = 1:M
        res = y_meas - f(t,[b1(i,j) b2(i,j)]);
%       SSR(i,j) = (W*res)'*(W*res);                    % Generalized least squares
%       SSR(i,j) = res'*inv(S2_eps)*res;                % Generalized least squares
%                = 1/s2_eps * res'*Vi*res
        SSR(i,j) = res'/S2_eps * res;                   % Generalized least squares
        % now also sort in different way
        G(counter,1:3) = [ b1(i,j) b2(i,j) SSR(i,j) ]; counter = counter + 1;
    end
end
G = sortrows(G,p+1);             
delta_SSR = chi2inv(gam,p);        % Now determine tolerable increment
% Now find which ones satisfy the threshold
gr_50id = (G(1:M^2,p+1)<(SSR_opt+delta_SSR(4)));
gr_90id = (G(1:M^2,p+1)<(SSR_opt+delta_SSR(3))) - gr_50id;
gr_95id = (G(1:M^2,p+1)<(SSR_opt+delta_SSR(2))) - gr_90id - gr_50id;
gr_99id = (G(1:M^2,p+1)<(SSR_opt+delta_SSR(1))) - gr_95id - gr_90id - gr_50id;
gr_50id = find(gr_50id==1); gr_90id = find(gr_90id==1); gr_95id = find(gr_95id==1);
gr_99id = find(gr_99id==1);
% --> from this output we can plot the ellipsoidal regions

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6: JAV approach
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[semimajor_axisT,semiminor_axisT,X0T,Y0T,phiT] = ...
    get_ellipse(gam,beta,dc,p,n,M,b,s2_eps,Vi,D,4,maxiter,tolfun);
% Q = Vi; and Z = D; --> this way we do not have to transpose original regression model
% MAIN EQUATIONS LISTED BELOW THAT GO WITH APPROACH JAV
c_ls = (D'*Vi*D)\(D'*Vi*y_meas);    % Least squares solution - sigma2_eps_max multiplier cancels
res = y_meas - f(t,c_ls);           % residual
s2_eps = res'*Vi*res/(n-p);         % variance of residual (should be sigma2_eps for n --> infty)
S2_eps = s2_eps * V;                % approximation of covariance matrix measurement errors
% C = s2_eps * inv(D'*Vi*D);          % variance-covariance matrix of parameters
% % dc = tinv(0.975,n-p)*sqrt(diag(C)); % range at 95%
% % ls_95 = [ c_ls-dc , c_ls+dc ];      % 95% ranges of parameters
% % % SSR_opt = res'*inv(S2_eps)*res      % weigted SSR with GLS parameter values
% %                                     % now determine 
% % % vareps = W*res/sqrt(s2_eps);      % numerator removes correlation from residuals
% Still need to look at normalized decorrelated residuals; their
% computation as V only does provide uncorrelated residuals but not with unit variance
% now lets look at projections on axix --> 
jv_95 = [c_ls - sqrt(diag(C))*sqrt(chi2inv(0.95,p)) , ...
    c_ls + sqrt(diag(C))*sqrt(chi2inv(0.95,p)) ];
jv_95


% Print to screen for now
[ ls_95 , zeros(p,1) , bs_95 , zeros(p,1) , jv_95 ]
fprintf('                        Results of different methods\n');
fprintf('------------- GLS method ------ Bootstrap ----- Projection ----- \n');
fprintf('              2.5%%  97.5%%     2.5%%    97.5%%    2.5%%   97.5%% \n');
fprintf('parameter 1: %4.3f  %4.3f     %4.3f   %4.3f    %4.3f  %4.3f \n',...
    ls_95(1,1:2),bs_95(1,1:2),jv_95(1,1:2));
fprintf('parameter 2: %4.3f  %4.3f     %4.3f   %4.3f    %4.3f  %4.3f \n',...
    ls_95(2,1:2),bs_95(2,1:2),jv_95(2,1:2));
fprintf('---------------------------------------------------------------- \n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  PLOT FIGURES 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch fig1 
    case 1
        plot_figure_1
end
switch fig11
    case 1
        plot_figure_11
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXTRA: CHECK OF DISTRIBUTION
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LETS LOOK AT NUMERATOR: THE ELLIPSOID
U = D'*Vi*D; % pxp matrix
for i = 1:N
    e(i,1) = (beta' - Bopt(i,1:p))*U*(beta - Bopt(i,1:p)');
end
% Make a histogram
Nbins = 30; figure(3)
[Nedges,edges] = histcounts(e,Nbins,'Normalization','pdf'); 
cnt_bin = 1/2*(edges(1:Nbins) + edges(2:Nbins+1));
subplot(1,2,1),bar(cnt_bin,Nedges/trapz(cnt_bin,Nedges));
% Now plot chi2 with p degrees of freedom
u = sigma2_eps*chi2rnd(p,N,1); 
[Nedges,edges] = histcounts(u,Nbins,'Normalization','pdf'); 
cnt_bin = 1/2*(edges(1:Nbins) + edges(2:Nbins+1));
subplot(1,2,2),bar(cnt_bin,Nedges/trapz(cnt_bin,Nedges));
%% This 'proofs' that the numerator is ~ sigma2_eps*chi2_(p)

%%% end of script