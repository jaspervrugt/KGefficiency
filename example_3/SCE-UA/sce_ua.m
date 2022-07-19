function [D,x,ParSet] = sce_ua(SCEPar,ModelName,ParRange,Measurement,Extra,option,plugin);
% The SHUFFLED COMPLEX EVOLUTION GLOBAL OPTIMIZATION METHOD (SCE-UA)
%
% Written by Jasper A. Vrugt, Tucson, AZ, July 2001
% ------------------------------------------------------------------------------------

% Define number of points in each complex
SCEPar.m = (2*SCEPar.n+1);

% Define number of points to draw Simplex inferences
SCEPar.q = SCEPar.n+1;

% How many offspring should be generated
SCEPar.Beta = (2*SCEPar.n+1);

% Bayesian inference scheme
[SCEPar.Wb,SCEPar.Cb] = calccbwb(SCEPar.Gamma);

% Computes the sample size (Step 1 in WRR Duan et al);
SCEPar.s = SCEPar.p * SCEPar.m; Iter = SCEPar.s;

% Do latin hypercube sampling of the initial space
[x] = lhsu(ParRange.minn,ParRange.maxn,SCEPar.s);

% Step 2b: Compute the objective function for each point x(i)
[OF] = compOF(x,SCEPar,ModelName,Measurement,Extra,option,plugin);

% Step 3: Sort the points in order of increasing function value
D = sortje(OF);

% Save the memory
ParSet = [x OF(:,1)];

% Initialize waitbar
% hwait = waitbar(0,'Shuffled Complex Evolution - Optimization Algorithm...');

% Now iterate
while Iter < SCEPar.ndraw,

    % Step 4: Partition D into p complexes
    [A] = part(D,SCEPar.p,x,SCEPar.n);
    % Loop over individual complexes
    for pi = 1:SCEPar.p
        % Loop in individual complexes
        for w = 1:SCEPar.Beta,
            % Step 5: Competitive complex evolution algorithm
            [rho] = asswght(SCEPar.m);
            % Step 5(3): Select parents 
            [B,L] = choose(rho,A(:,:,pi),SCEPar.n,SCEPar.q); 
            % Create offspring
            [B,L,ParSet,Iter] = offspring(A(:,:,pi),B,L,SCEPar,ModelName,ParRange,Iter,Measurement,Extra,ParSet,option,plugin);               
            % Step 5(5): Replace parents by offspring
            [A(:,:,pi)] = sce_replace(A(:,:,pi),B,L,SCEPar.n);                                               
        end;
    end;
    
    if 0 
%     if size(ParSet,1)>250
        visualize(ParSet,ParRange)
    end
    
    % Step(6): Reshuffle Complexes
    [D,x] = reshuffle(A,SCEPar.p,SCEPar.n,SCEPar.m);                                                     
    % Update waitbar
%     waitbar(Iter/SCEPar.ndraw); 
    disp([num2str(100.*Iter/SCEPar.ndraw,'%.0f') '%']) 
end;

% Remove the waitbar
% close(hwait);       


function visualize(ParSet,ParRange)

figure(111)
set(gcf,'Position',[50 50 1000 600])

% number of parameters
np = length(ParRange.minn);

[nRows,nCols] = size(ParSet);

subplot(1,2,1)
hPlotObj = plot(1:nRows,ParSet(1:nRows,nCols),'bs');
set(hPlotObj,'markersize',3,'markerfacecolor','b','markeredgecolor','b')
xlabel('iteration')
ylabel('objective score')

hPlotPar = repmat(NaN,[np,1]);
m=1;
for k=2:2:2*np

    subplot(np,2,k)
    hPlotPar(m) = plot(1:nRows,ParSet(1:nRows,m),'ks');
    set(hPlotPar(m),'markerfacecolor','k',...
        'markeredgecolor','k',...
        'markersize',3)
    set(gca,'ylim',[ParRange.minn(m),ParRange.maxn(m)])
            
    m=m+1;

end
drawnow