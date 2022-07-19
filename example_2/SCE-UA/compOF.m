function [pset] = compOF(x,SCEPar,ModelName,Measurement,Extra,option,plugin);
% Compute the objective function value

% Unpack the measurement data 

MeasData = Measurement;
% MeasData = Measurement.MeasData; 
% Sigma = Measurement.Sigma;

% Unpack the parameters of the exponential power density function
Cb = SCEPar.Cb; 
Wb = SCEPar.Wb;
% Compute the number of measurements
N = length(MeasData(:));
% Compute the size of x
[NrParSets,r] = size(x);

% For each samplecompute function value and calculate posterior density ...
for ii=1:NrParSets,
    % Call model to generate simulated data
    evalstr = [ModelName,'(x(ii,:),plugin);'];
    SimData = eval(evalstr);

    if option == 1,
        % Model directly computes posterior density
        p = SimData;
        pset(ii,1:2) = [p ii];
    end;

    if option == 2,
        % Compute log likelihood value: Taking the logarithm is a good 
        % way of handling small values for large N
        lnp = N.*log(Wb./Sigma) - Cb.*(sum((abs(Err./Sigma)).^(2/(1+SCEPar.Gamma))));
        pset(ii,1:2) = [lnp ii];
    end;

    if option == 3, % RMSE

        a = [Extra.calPeriod(1):Extra.calPeriod(2)]';
        b = MeasData(a);
        c = SimData;

        if Extra.plotYN
            figure(222)
            subplot(3,1,3)
            plot(Extra.numTime(a),[b,c])
            ylabel('discharge')
            set(gca,'yscale','log')
            datetick('x',1)
            drawnow
        end
        
        % Sum of squared error likelihood function
        Err = b-c;
        % Part of equation 8 from the manual
        M = sum(abs(Err).^(2/(1+SCEPar.Gamma))); 
        % Write RMSE to screen to follow progress
%        RMSE = sqrt(M/N)
        % Minus is added so minimum M corresponds to highest likelihood
        pset(ii,1:2) = [1./M ii];
        
    end;
    
    if option == 4,
        % Model directly computes log posterior density
        p = SimData;
        pset(ii,1:2) = [p ii];
    end;
 
    if option == 5, % KGE

        a = [Extra.calPeriod(1):Extra.calPeriod(2)]';
        b = MeasData(a); 
        c = SimData; 
        M = parKGE(c,b);

        if Extra.plotYN
            figure(222)
            subplot(3,1,3)
            plot(Extra.numTime(a),[b,c])
            ylabel('discharge')
            set(gca,'yscale','log')
            datetick('x',1)
            drawnow
        end
       
        pset(ii,1:2) = [M ii];
        
    end;
    
    if option == 6, % heteroscedastic likelihood

        a = [Extra.calPeriod(1):Extra.calPeriod(2)]';
        b = MeasData(a); 
        c = SimData; 
        res = b - c;
        invC = Extra.invC;
        M = - 1/2 * sum(res.*(invC*res));
        
        if Extra.plotYN
            figure(222)
            subplot(3,1,3)
            plot(Extra.numTime(a),[b,c])
            ylabel('discharge')
            set(gca,'yscale','log')
            datetick('x',1)
            drawnow
        end
        
        pset(ii,1:2) = [M ii];
        
    end;
    
end;