function [lowerBound, upperBound] = ConfidenceIntervalBoundsLogEpsV2(list,N,logEpsilon)
    arguments
        list (:,:,:) double {mustBeNonnegative, mustBeLessThanOrEqual(list,1)}
        N (1,1) double {mustBePositive}
        logEpsilon (1,1) double {mustBeNegative}
    end
    % list is a list of observed propobailities of events 
    % N is total number of trials
    % logAlpha is the confidence interval
    
    %Store shape of input list
    shape = size(list);

    if logEpsilon >= -12
        %Covert log of alpha back into alpha 
        alpha = 10^(logEpsilon);

        %observations, we round to nearest integer
        kobs = round(N*list(:));
        
        %Bounds given by Clopper pearson intervals from Matlab
        [~,bounds] = binofit(kobs,N,alpha); 
        
        lowerBounds = bounds(:,1);
        upperBounds = bounds(:,2);
        
        %Reshape to original shape of input
        lowerBound = reshape(lowerBounds, shape);
        upperBound = reshape(upperBounds, shape);

    else
        %Implement bounds using a Gaussian bound on binomial distribution
        %from: A Complete Proof of Universal Inequalities for the Distribution
        % Function of the Binomial Law 
        % https://doi.org/10.1137/S0040585X97986138

        %Store symbolic version of confidence level
        epsilonSym = 10^(sym(logEpsilon));

        %recast list into a flattened vector
        Fflat = list(:);
        
        %predefine vectors
        lowerBound = zeros(size(Fflat));
        upperBound = ones(size(Fflat));
        
        for index=1:numel(Fflat)        
            %Calculate initial mu
            [lowerBound(index), upperBound(index)] = binBoundsGauss(Fflat(index),N,epsilonSym);
        end

        %Reshape to original shape of input
        lowerBound = reshape(lowerBound, shape);
        upperBound = reshape(upperBound, shape);
    end

end

function [pLower,pUpper] = binBoundsGauss(F,N,alphaSym)
        %Calculate inverse of normal distribution 
        phiInvAlpha = -sqrt(2)*erfcinv(alphaSym);
        phiInv1MinusAlpha = -sqrt(2)*erfcinv(2*(1-alphaSym/2));

        %Convert N into double
        Ndbl = double(N);

        %% Upper Bound 
        %Upper range of the accepted parameters
        kUpper = double(ceil(N*F));
        cUpper = double(phiInvAlpha^2/2);

        funUpper = @(x) ((kUpper+1)*(log((kUpper+1)/Ndbl) - log(x)) + ...
                (Ndbl - kUpper - 1)*(log(1-(kUpper+1)/Ndbl) - log(1-x)) - cUpper);

        if kUpper >= N
            pUpper = 1;
        elseif kUpper <= 0
            pUpper = min(max(double(1 - nthroot(alphaSym/2,sym(N))),0),1);
        else
            try
                pUpper = min(max(fzero(funUpper,[(kUpper+1)/Ndbl,1-10^(-14)]),0),1);
            catch
                warning('Problem using fzero. Assigning a value of 1 to pUpper.');
                pUpper = 1;
            end          
        end

        %% Lower Bound 
        kLower = double(floor(N*F));
        cLower = double(phiInv1MinusAlpha^2/2);

        funLower = @(x) ((kLower-1)*(log((kLower-1)/Ndbl) - log(x)) + ...
                (Ndbl - kLower + 1)*(log(1-(kLower-1)/Ndbl) - log(1-x)) - cLower);
        
        if kLower >= N
            pLower = min(max(double(nthroot(alphaSym/2,sym(N))),0),1);
        elseif kLower <= 1
            pLower = 0;
        else
            try
                pLower = min(max(fzero(funLower,[eps(0),(kLower-1)/Ndbl]),0),1);
            catch
                warning('Problem using fzero. Assigning a value of 0 to pLower.');
                pLower = 0;
            end
        end

end
