function [eps_PA, eps_EC, eps_AT, epsBar] = optimalEpsilon(epsilonSecurity,gains,params)
%Calculates the optimal values for epsilon PA eps_EC, eps_AT, epsBar
%by finding the root of the derivative of l wrt epsilon PA
% epsilonSecurity:  total security parameter
% nsift:            number of detected signals usable for key generation,
%                   e.g. detection and basis match
% N:                total number of signals sent
% t:                acceptance parameter
% dimX:             dimension of the alphabet X
    
    %testing probability
    ptest = params.ptest;
    
    %key generation probability
    pgen = 1 - ptest;

    %total signals sent
    N = params.N;
    
    %t parameter for allowed fluctuations
    t = params.tsift;
    
    %nsift
    totalGain = sum(gains,1);
    nsift = N*pgen*totalGain;

    %Dimension of alphabet
    dimX = params.alphabet;
    
    %define drivative as single valued function
    der = @(x) derivative_l_normalised(x,epsilonSecurity,nsift,N,t,dimX);

    %Find root
    optz = fzero(der,[eps(0),1]);
    eps_PA = epsilonSecurity*(1-optz) - 8*epsilonSecurity^2*(1-optz)^2;

    %Calculate remaining epsilon values
    % Smooting parameter
    epsBar = 1/2*(epsilonSecurity - 8*eps_PA^2-eps_PA);
    %Error correction parameter
    eps_EC = 8*eps_PA^2;
    %Acceptance test parameter (or parameter estimation)
    eps_AT = eps_PA + 2*epsBar;
end

function dl = derivative_l_normalised(z,epsilonSecurity,nsift,N,t,dimX)
%Calculates derivative of key length with respect to epsilon_PA
% eps_PA:           security parameter from privacy amplification
% epsilonSecurity:  total security parameter
% nsift:            number of detected signals usable for key generation,
%                   e.g. detection and basis match
% N:                total number of signals sent
% t:                acceptance parameter
% dimX:             dimension of the alphabet X

    c = sqrt(N)*sqrt(max(nsift/N-t,eps(0))).*2.*log(1+dimX);
    dl = z.*sqrt(1-2*log2(z)-2*log2(1/2*epsilonSecurity)) + z.*c - c;
end