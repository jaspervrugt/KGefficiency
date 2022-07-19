function [B,L,ParSet,Iter] = offspring(A,B,L,SCEPar,ModelName,ParRange,Iter,Measurement,Extra,ParSet,option,plugin);
% Generates offspring using complex evolution algorithm

% Loop over the number of offspring to be generated in each complex
for w = 1:SCEPar.alpha,
    % Sort B and L so that the q points are arranged in order of increasing function value
    B = -sortrows(-B,[SCEPar.n+1]); L = -sortrows(-L,[2]);
    % Compute centroid of original vector B,
    g = mean(B(1:SCEPar.q-1,1:SCEPar.n));
    % Define the worst point in the simplex
    uq = B(SCEPar.q,1:SCEPar.n); fq = L(SCEPar.q,2); OF = fq;
    %Compute reflection step,
    r = (2*g - uq);
    % Test whether r is within the feasible parameter space
    if sum((r-ParRange.minn) > 0) == SCEPar.n & sum((ParRange.maxn-r) > 0) == SCEPar.n,
        % Compute the objective function value
        [fr] = compOF(r,SCEPar,ModelName,Measurement,Extra,option,plugin); Iter = Iter + 1;
        % Save ParSet
        ParSet = [ParSet; r fr(1)];
    else
        % Mutation step
        [r] = mutation(A,SCEPar.n);
        % Compute the objective function value
        [fr] = compOF(r,SCEPar,ModelName,Measurement,Extra,option,plugin); Iter = Iter + 1;
        % Save ParSet
        ParSet = [ParSet; r fr(1)];
    end;
    % Test whether new Simplex is improvement
    if (fr(1) >= fq),
        uq = r; OF = fr(1);
    else
        % Compute contraction step,
        c = (g + uq)/2;
        % Compute the objective function value
        [fc] = compOF(c,SCEPar,ModelName,Measurement,Extra,option,plugin); Iter = Iter + 1;
        % Save ParSet
        ParSet = [ParSet; c fc(1)];
        % Now do evaluations
        if (fc(1) >= fq),
            uq = c; OF = fc(1);
        else
            % Generate mutation step
            [z] = mutation(A,SCEPar.n);
            % Compute the objective function value
            [fz] = compOF(z,SCEPar,ModelName,Measurement,Extra,option,plugin); Iter = Iter + 1;
            % Save ParSet
            ParSet = [ParSet; z fz(1)];
            % And accept
            uq = z; OF = fz(1);
        end;
    end;
    B(SCEPar.q,1:SCEPar.n+1) = [uq(1,1:SCEPar.n) OF];
end;