function [lowerbnd, upperbnd, indbnd] = ProbAcceptCalc(N,expec,siftProb,t,tsift)
    upperbnd = ProbUpperAccept(N,expec,siftProb,t,tsift);
    lowerbnd = ProbLowerAccept(N,expec,siftProb,t,tsift);
    indbnd = ProbIndAccept(N,expec,siftProb,t,tsift);
end


function upperbnd = ProbUpperAccept(N,expec,siftProb,t,tsift)
    %Calculates upper bound on accept probability for entry-wise constraints
    
    % N:        Total number of rounds
    % expec:    Full matrix of expectations in test rounds, i.e. p(a,b,test)
    % siftProb: Probability of detecting a signal and surviving sifting Sum_a p(a,det+sift,gen)
    % t:        tolerated fluctuations in test rounds
    % tsift:    tolerated fluctuation from sifted signals in generation rounds

    %Iterate over single events

    %% Dev edit : The exep data is 3 dimensions (intensity choices).
    % so we unroll 
    mat_single_events = zeros(numel(expec),1);

    for index = 1 : numel(expec)
         mat_single_events(index) = P_single_succ(N,expec(index),t);
    end
    
    %Accept probability for sifting calculated separately
    probAcceptSift = P_sift_accept(N,siftProb,tsift);

    vecSingleSuccess = [mat_single_events; probAcceptSift];

    %Calculate lower bound
    upperbnd = min(vecSingleSuccess,[],"all");
end

function lowerbnd = ProbLowerAccept(N,expec,siftProb,t,tsift)
    %Calculates lower bound on accept probability for entry-wise constraints

    % N:        Total number of rounds
    % expec:    Full matrix of expectations in test rounds, i.e. p(a,b,test)
    % siftProb: Probability of detecting a signal and surviving sifting Sum_a p(a,det+sift,gen)
    % t:        tolerated fluctuations in test rounds
    % tsift:    tolerated fluctuation from sifted signals in generation rounds

    % %Iterate over single events
    mat_single_events = zeros(numel(expec),1);
    for index = 1 :numel(expec)
        mat_single_events(index) = P_single_succ(N,expec(index),t);
    end

    vecSingleSuccess = [mat_single_events; P_single_succ(N,siftProb,tsift)];

    %Calculate lower bound
    lowerbnd = max(0,sum(vecSingleSuccess,"all") - (numel(vecSingleSuccess)-1));
end

function indbnd = ProbIndAccept(N,expec,siftProb,t,tsift)
    %Calculates accept probability for entry-wise constraints assuming
    %independent events

    % N:        Total number of rounds
    % expec:    Full matrix of expectations in test rounds, i.e. p(a,b,test)
    % siftProb: Probability of detecting a signal and surviving sifting Sum_a p(a,det+sift,gen)
    % t:        tolerated fluctuations in test rounds
    % tsift:    tolerated fluctuation from sifted signals in generation rounds

    %Iterate over single events
    mat_single_events = zeros(numel(expec),1);
    for index = 1 :numel(expec)
        mat_single_events(index) = P_single_succ(N,expec(index),t);
    end

    vecSingleSuccess = [mat_single_events; P_single_succ(N,siftProb,tsift)];


    %Calculate lower bound
   indbnd = prod(vecSingleSuccess,"all");
end

function value = helper_minus(x)
    %helper function to make array only have positive entries
    value = zeros(size(x));
    for index = 1:length(x)
        if x(index) < 0
            value(index) = 0;
        else
            value(index) = x(index);
        end
    end
end

function prob = P_single_succ(N,p,t)
    %Calculates accept probability Pr(|X/N - p| <= t)
    % given a probability p for input t and N

    %N: Total number of rounds
    %p: probability
    %t: Acceptance parameter

    %tolerance to 0
    tol = 1e-14;

    %Upper range of the accepted parameters
    m1 = floor(N*(p+t));

    %Lower range of the accepted parameters
    m2 = helper_minus(ceil(N*(p-t)));

    if p <= tol
        prob = 1;
    elseif abs(p-1) <= tol
        prob = 1;
    else
        prob = betainc(1-p, helper_minus(N-m1), m1+1) - betainc(1-p, helper_minus(N-m2+1), m2);
    end
end

function prob = P_sift_accept(N,p,t)
    %Calculates accept probability Pr( X/N >= p - t)
    % given a probability p for input t and N

    %N: Total number of rounds
    %p: probability
    %t: Acceptance parameter

    %tolerance to 0
    tol = 1e-14;

    %Lower range of the accepted parameters
    m2 = helper_minus(ceil(N*(p-t)));

    if p <= tol
        prob = 1;
    elseif abs(p-1) <= tol
        prob = 1;
    else
        prob = 1 - betainc(1-p, helper_minus(N-m2+1), m2);
    end
end