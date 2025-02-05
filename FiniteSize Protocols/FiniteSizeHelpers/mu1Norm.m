function fval = mu1Norm(epsilon,Sigma,m)
    fval = sqrt(2/m)*sqrt(log(1/epsilon.PE) + Sigma*log(m+1));  
end