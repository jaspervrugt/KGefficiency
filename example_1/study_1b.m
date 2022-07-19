%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% Vrugt, J.A. and de Oliveira, D.Y., 2022. Confidence intervals of the
%%  Kling-Gupta efficiency, Journal of Hydrology, 612, 
%%  10.1016/j.jhydrol.2022.127968
%% https://doi.org/10.1016/j.jhydrol.2022.127968
%%
%% GLS versus KG efficiency for linear regression function
%% 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Synthetic case study: Create measurement data using one of four methods
approach = 1;   % [1] homoscedastic error only
                % [2] heteroscedastic error only
                % [3] homoscedastic + correlated error
                % [4] heteroscedastic + correlated error
load_data = 2;  % [1] remake data 
                % [2] load data of trial (avoid variations in fig settings due to randomness)

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
                    R(i,idr) = phi.^([1:numel(idr)]);   % upper triangle
                end
                V = R;
                switch approach
                    case 3 % V matrix is done
                        sigma2_eps = 0.5;
                    case 4 % Must update diagonal of V to account for heteroscedasticity
                        a = 0.1; b = 0.01;
                        sigma2_eps_all = (a*abs(y) + b).^2;
                        sigma2_eps = min(sigma2_eps_all);
                        for i = 1:n, V(i,i) = sigma2_eps_all(i)/sigma2_eps; end
                        clear sigma2_eps_all; 
                end
        end
        Sigma2_eps = sigma2_eps * V;    % Create measurement error covariance matrix
        Vi = inv(V);                    % Get the inverse of matrix V 
       % P = chol(V);                    % Now get P'*P=V or P*P=V or P^2 = V
        meas_err = mvnrnd(zeros(n,1),...    
            Sigma2_eps,1)';             % Realization of measurement errors
        y_meas = y + meas_err;          % Measurements of dgp
       % evalstr = strcat('save',{' '},'data_',...
       % num2str(approach),{' '},'t n p meas_err D c_dgp Sigma2_eps ...
       % sigma2_eps V y_meas'); eval(char(evalstr));
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
s2_eps = res'*Vi*res/(n-p);         % variance of residual (should be sigma2_eps for n --> infty)
S2_eps = s2_eps * V;                % approximation of covariance matrix measurement errors
C = s2_eps * inv(D'*Vi*D);          % variance-covariance matrix of parameters
   % = inv(D'*inv(S2_eps)*D)
   % = s2_eps*inv(Q'*Q)
dc = tinv(0.975,n-p)*sqrt(diag(C)); % range at 95%
ls_95 = [ c_ls-dc , c_ls+dc ];      % 95% ranges of parameters
SSR_opt = res'*inv(S2_eps)*res;     % weigted SSR with GLS parameter values

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0: Sum of squared residuals + KG for contour plot
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define xmin etc. works here
xmin = 0.97; xmax = 1.025; ymin = 1.30; ymax = 2.60;
b(1,1:M) = xmin:(xmax-xmin)/(M-1):xmax;
b(2,1:M) = ymin:(ymax-ymin)/(M-1):ymax;
s_meas = std(y_meas); m_meas = mean(y_meas); 
beta = c_ls;
[b1,b2] = meshgrid(b(1,1:M),b(2,1:M)); SSR = nan(M,M); KGE_cont = nan(M,M);
for i = 1:M
    for j = 1:M
        y_sim = f(t,[b1(i,j) b2(i,j)]);
        res = y_meas - y_sim;
%        SSR(i,j) = (W*res)'*(W*res);                   % Generalized least squares
%        SSR(i,j) = res'*inv(S2_eps)*res;               % Generalized least squares
%                 = 1/s2_eps * res'*Vi*res
        SSR(i,j) = res'/S2_eps * res;                   % Generalized least squares
        s_sim = std(y_sim); m_sim = mean(y_sim); 
        alpha = s_sim/s_meas; gamma = m_sim/m_meas;
        r = corrcoef(y_sim,y_meas); r = r(1,2);
        % Compute KG efficiency
        KGE_cont(i,j) = 1 - sqrt((r-1)^2 + (alpha-1)^2 + (gamma-1)^2);
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1: JAV approach
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_gam = numel(gam);                             % number of alfa values
method = 4; maxiter = 1000; tolfun = 1e-10;     % Define maxiter and tolfun for secant method
[semimajor_axisT,semiminor_axisT,X0T,Y0T,phiT] = ...
    get_ellipse(gam,c_ls,dc,p,n,M,b,s2_eps,Vi,D,4,maxiter,tolfun);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2: KG efficiency: Bootstrap
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[c_kg,kgmax] = getoptKGE(f,t,y_meas);           % Now get optimum for KG efficiency
s_meas = std(y_meas); m_meas = mean(y_meas);    % Now get statistics of measured record for KGE 
KGEopt = nan(p+2,N); resb = nan(n,N);           % Bootstrap method
R = y_meas + mvnrnd(zeros(n,1),S2_eps,N)';      % Create ensemble of N measured records
for i = 1:N
    [ab_kg,kg] = getoptKGE(f,t,R(1:n,i)); KGEopt(1:p+1,i) = [ab_kg kg]';
    % Evaluate KGE for measured record
    y = f(t,KGEopt(1:p,i)); s_sim = std(y); m_sim = mean(y); 
    r = corrcoef(y,y_meas); 
    alpha = m_sim/m_meas; gamma = s_meas/s_sim; r = r(1,2);
    KGEopt(p+2,i) = 1 - sqrt((r-1)^2 + (alpha-1)^2 + (gamma-1)^2);
end
KG = -sortrows(-KGEopt',p+2);   % Get percentiles
ii_perc = round(N*gam);         % Now pick percentiles
kg_95 = nan(p,2);               % Now get 95% pertaining to sigma^2 of data
id_lu = round([0.025 0.975]*N); % 2.5 and 97.5 percentile
for j = 1:p
    a = sort(KGEopt(j,1:N)); kg_95(j,1:2) = [a(id_lu(1)) a(id_lu(2))];
end
% Get joint 95% confidence region of parameters
id_up = round(gam*N);           % gamma (percentiles)
% Now find which ones satisfy the threshold
kg_50id = 1:id_up(4);
kg_90id = id_up(4)+1:id_up(3);  % no double plotting of points
kg_95id = id_up(3)+1:id_up(2);  % no double plotting of points
kg_99id = id_up(2)+1:id_up(1);  % no double plotting of points

%%% end of script