function [semimajor_axis,semiminor_axis,X0,Y0,phi] = ...
    get_ellipse(gam,beta,dc,p,n,M,b,s2_eps,Q,Z,method,maxiter,tolfun)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% Solves for the extreme points of the 100gamma confidence region and fits
%% ellipsoids through the data points
%% 
%% Vrugt, J.A. and de Oliveira, D.Y., 2022. Confidence intervals of the
%%  Kling-Gupta efficiency, Journal of Hydrology, 612, 
%%  10.1016/j.jhydrol.2022.127968
%% https://doi.org/10.1016/j.jhydrol.2022.127968
%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_gam = numel(gam);         % how many alfa values
B_gam = cell(n_gam,1);      % Initialize cell array with points of ellipse 
for z = 1:n_gam
    B = nan(2*M,p+1);       % Initialize matrix B with b1, b2 and 
                            % exitflag for lower + upper perimeter
    for u = 1:2
        switch u
            case 1                  % lower end ellipse -> walk up secant
                bp_init = beta(p)-4*dc(p);
            case 2                  % upper end ellipse -> walk down secant
                bp_init = beta(p)+4*dc(p);
        end
        [bp,exitflag] = deal(nan(M,1)); % Initialize b2 and exitflag as 
                                        % we do loop 2 times
        for i = 1:M                     % Loop over the different b1 values 
                                        % to get b1,b2 pairs of perimeter
            [bp(i,1),exitflag(i,1)] = secant(bp_init,maxiter,tolfun,...
                b(1:p-1,i),beta,s2_eps,Q,Z,gam(z),p,n-p,method);
        end
        % pool points together for figure
        B((u-1)*M+1:u*M,1:p+1) = [ b(1:p-1,1:M)' bp exitflag ]; 
    end
    % keep only those pairs of b1,b2 with exitflag of secant equal to 1
    idx = B(:,p+1) == 1; B_gam{z} = B(idx,1:p+1);
end

% Now derive the axes of the ellipse
[semimajor_axis,semiminor_axis,X0,Y0,phi] = deal(nan(1,n_gam)); 
% Get angle, major + minor axis by fitting ellipse to B_gam{z} data pairs
for zz = 1:n_gam
    [semimajor_axis(zz), semiminor_axis(zz), X0(zz), Y0(zz), phi(zz)] = ...
        ellipse_fit( B_gam{zz}(:,1) , B_gam{zz}(:,2) );   
end