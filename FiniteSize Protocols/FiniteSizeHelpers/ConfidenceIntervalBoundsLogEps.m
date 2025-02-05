function [lowerBound, upperBound] = ConfidenceIntervalBoundsLogEps(list,N,logAlpha)
    arguments
        list (:,:,:) double {mustBeNonnegative, mustBeLessThanOrEqual(list,1)}
        N (1,1) double {mustBePositive}
        logAlpha (1,1) double {mustBeNegative}
    end
    % list is a list of observed propobailities of events 
    % N is total number of trials
    % logAlpha is the confidence interval
    
    %Store shape of input list
    shape = size(list);

    if logAlpha >= -12
        %Covert log of alpha back into alpha 
        alpha = 10^(logAlpha);

        %observations, we round to nearest integer
        kobs = round(N*list(:));
        
        %Bounds given by Clopper pearson intervals from Matlab
        [~,bounds] = binofit(kobs,N,alpha); 
        
        lowerBounds = bounds(:,1);
        upperBounds = bounds(:,2);
        
        %Respahe to original shape of input
        lowerBound = reshape(lowerBounds, shape);
        upperBound = reshape(upperBounds, shape);

    else
        %Implement bounds using a Gaussian bound on binomial distribution
        %from On binomial quantile and proportion bounds: With applications in engineering and informatics
        %https://doi.org/10.1080/03610926.2021.1986540

        %Store symbolic version of confidence level
        alphaSym = 10^(sym(logAlpha));
        R = (1 - alphaSym);
        
        %Calculate inverse of normal distribution of R
        phiInvR = -sqrt(2)*erfcinv(2*R);

        %Observed number of outcomes
        kobs = sym(round(N*list(:)));

        %% Upper Bound see Corollary 3 following same notation
        %Constants
        cU1 = sym(5/6)*phiInvR^2 + 1;
        cU2 = sym(7/12)*phiInvR^2 + 1;
        cU3 = sym(2/3)*phiInvR^2 - 2;
        cU4 = sym(1/9)*phiInvR^4 + 2/3*phiInvR^2 + 1;
        cU5 = phiInvR^2;
    
        %Calculate upper bound and return it as double
        %If kobs = N, upperbound = 1, because formula will result in
        %complex values
        upperBound = zeros(size(kobs));

        for index=1:numel(kobs)
            kobsIterate = kobs(index);
            if kobsIterate < N
                upperBound(index) = (kobsIterate + cU1)./(N + cU5) +...
                    sqrt(cU5).*sqrt((N*kobsIterate - kobsIterate.^2 + cU2*N - cU3*kobsIterate - cU4)./(N*(N+cU5)^2));
            else
                upperBound(index) = 1;
            end
        end
        
        %Reshape to original shape of input
        upperBound = reshape(double(upperBound),shape);

        %% Lower Bound see Corollary 4 following same notation
        %Constants
        cL1 = sym(1/6)*phiInvR^2 + 1;
        cL2 = sym(1/12)*phiInvR^2 - 4;
        cL3 = sym(2/3)*phiInvR^2 - 4;
        cL4 = sym(1/9)*phiInvR^4 + 2;
        cL5 = phiInvR^2;
       
        %Calculate lower bound and return it as double
        lowerBound = zeros(size(kobs));

        %If kobs = 0, lowerbound = 0, because formula will result in
        %complex values
        for index=1:numel(kobs)
            kobsIterate = kobs(index);
            if kobsIterate > 0
                lowerBound(index) = (kobsIterate + cL1)./(N + cL5) - ...
                    sqrt(cL5).*sqrt((N*kobsIterate - kobsIterate.^2 - cL2*N + cL3*kobsIterate - cL4)./(N*(N+cL5)^2));
            else 
                lowerBound(index) = 0;
            end
        end
        
        %Reshape to original shape of input
        lowerBound = reshape(double(lowerBound),shape);
    end

end
