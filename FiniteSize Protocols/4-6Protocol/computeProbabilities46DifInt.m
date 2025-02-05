function [pA, pTestCondSingle, pGenCondSingle, pTestCondSingleandA, pGenCondSingleandA, pACondTestandSingle, pACondSingle, pSinglePhoton,pSingleCondGen,pMuandACondTest] = computeProbabilities46DifInt(pGen, pTest, probSignalsGeneration, probSignalsTesting, decoyMatrix, decoyProbsTest, decoyProbsGen)
    %Takes P( A = 1,2,3 | G = g/t ), P(g/t) and decoy signal probabilities and intensities as argument
    %Returns P( G = g/t | A = 1,2,3, n = 1)
    %The one time R would be useful :(
    
    %P( A = 1,2,3,4 | G = gen )
    p1CondGen = probSignalsGeneration(1); 
    p2CondGen = probSignalsGeneration(2); 
    p3CondGen = probSignalsGeneration(3);
    p4CondGen = probSignalsGeneration(4); 

    %P( A = 1,2,3,4 | G = test )
    p1CondTest = probSignalsTesting(1); 
    p2CondTest = probSignalsTesting(2); 
    p3CondTest = probSignalsTesting(3); 
    p4CondTest = probSignalsTesting(4); 
    
    %p(a)
    pA = [p1CondTest*pTest + p1CondGen*pGen; p2CondTest*pTest + p2CondGen*pGen; p3CondTest*pTest + p3CondGen*pGen;p4CondTest*pTest + p4CondGen*pGen];

    %photonProbDensity = @Thermal;
    photonProbDensity = @Poisson;

    %P(n=1 | a,test)
    pSingleCondAandTest = zeros(max(numel(probSignalsTesting),numel(probSignalsGeneration)),1);
    for indexA=1:numel(probSignalsTesting)
        IntensitiesA = decoyMatrix(:,indexA);
        pSingleCondIntensity = photonProbDensity(IntensitiesA,1);
        pSingleCondAandTest(indexA) = cell2mat(decoyProbsTest)*pSingleCondIntensity;
    end
    
    %P(n=1 | a,gen)
    pSingleCondAandGen = zeros(max(numel(probSignalsTesting),numel(probSignalsGeneration)),1);
    for indexA=1:numel(probSignalsGeneration)
        IntensitiesA = decoyMatrix(:,indexA);
        pSingleCondIntensity = photonProbDensity(IntensitiesA,1);
        pSingleCondAandGen(indexA) = cell2mat(decoyProbsGen)*pSingleCondIntensity;
    end

    %P(n=1 | gen)
    pSingleCondGen = probSignalsGeneration*pSingleCondAandGen;

    %P(n=1 | test)
    pSingleCondTest = probSignalsTesting*pSingleCondAandTest;

    %P(n=1)
    pSinglePhoton = pTest*pSingleCondTest + pGen*pSingleCondGen;

    %P(test|n=1) = P(1|test)*p(test)/p(1)
    pTestCondSingle = pSingleCondTest*pTest/pSinglePhoton;

    %P(gen|n=1) = P(1|gen)*p(gen)/p(1)
    pGenCondSingle = pSingleCondGen*pGen/pSinglePhoton;

    %P(a,n=1|test)
    pAandSingleCondTest = pSingleCondAandTest.*(probSignalsTesting.');

    %P(a,n=1|gen)
    pAandSingleCondGen = pSingleCondAandGen.*(probSignalsGeneration.');

    %P(a|test,n=1)
    pACondTestandSingle = pAandSingleCondTest/pSingleCondTest;

    %P(a|gen,n=1)
    pACondGenandSingle = pAandSingleCondGen/pSingleCondGen;

    %P(a|n=1)
    pACondSingle = pACondGenandSingle*pGenCondSingle + pACondTestandSingle*pTestCondSingle;

    %P(a,n=1)
    pAandSingle = pTest*pAandSingleCondTest + pGen*pAandSingleCondGen;

    %P(test|a,1)
    pTestCondSingleandA = pAandSingleCondTest.*pTest./pAandSingle;

    %P(gen|a,1)
    pGenCondSingleandA = pAandSingleCondGen.*pGen./pAandSingle;
    
    %P(mu,a|test)
    pMuandACondTest = zeros([length(probSignalsTesting), length(cell2mat(decoyProbsTest))]);
    for k = 1:length(cell2mat(decoyProbsTest))
        pMuandACondTest(:, k) = probSignalsTesting*decoyProbsTest{k};
    end
end

function p = Poisson(intensity, n)
    p = exp(-intensity).*intensity.^n./factorial(n);
end
function p = Thermal(intensity, n)
    p = intensity.^n./(intensity + 1).^(n+1);
end