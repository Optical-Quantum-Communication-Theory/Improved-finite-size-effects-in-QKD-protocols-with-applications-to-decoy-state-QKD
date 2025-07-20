function [keyRate, modParser] = FourSixWCPDifIntAdaptiveKeyRateFunc(params,options,mathSolverFunc,debugInfo)
% WLOG, this function works with the assumption that the z basis is used
% exclusively for generation rounds and the x basis is used exclusively for
% parameter estimation in the finite domain. Currently hard coded for two
% basis choices z and x.
% Aodhan Corrigan
%
% We also compute Adaptive key rates.
%% Input parameters:
% * dimA: dimension of Alice's system.
% * dimB: dimension of Bob's system.
% * announcementsA: Array of announcements made for each measurement Alice
%   made ordered in the same way as the columns of expectationsJoint.
% * announcementsB: Array of announcements made for each measurement Bob
%   made ordered in the same way as the rows of expectationsJoint.
% * keyMap: An array of KeyMapElement objects that contain pairs of accepted
%   announcements and an array dictating the mapping of Alice's measurement
%   outcome to key bits (May be written with Strings).
% * krausOps: A cell array of matrices. The Kraus operators that form the G
%   map on Alice and Bob's joint system. These should form a completely
%   postive trace non-increasing linear map. Each Kraus operator must be
%   the same size.
% * keyProj:  A cell array of projection operators that perform the pinching map 
%   key on  G(\rho). These projection operators should sum to identity.
% * f: error correction effiency. Set to 1 means for Shannon limit. 
% * observablesJoint: The joint observables of Alice and Bob's
%   measurments. The observables must be hermitian and each must be the size 
%   dimA*dimB by dimA*dimB. The observables assume the spaces are ordered A \otimes B.
%   They also should be positive semi-definite and should sum to identity. 
% * expectationsJoint: The joint expectations (as an array) from Alice and
%   Bob's measurements that line up with it's corresponding observable in
%   observablesJoint. These values should be betwen 0 and 1.
% * rhoA (nan): The fixed known density matrix on Alice's side for
%   prepare-and-measure protocols.
%% Following parameters are specific to Finite protocol %%
% * epsilons: Struct containing the values of different epsilons
% * 
%% Outputs:
% * keyrate: Key rate of the QKD protocol measured in bits per block
%   processed.
% Options:
% * verboseLevel: (global option) See makeGlobalOptionsParser for details.
% DebugInfo:
% * 
% 
% Reviewed by Devashish Tupkary 2023/09/18
% See also QKDKeyRateModule, FiniteBB84_4DAliceDescriptionFunc, makeGlobalOptionsParser
arguments
    params (1,1) struct
    options (1,1) struct
    mathSolverFunc (1,1) function_handle
    debugInfo (1,1) DebugInfo
end

%% options parser
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.addOptionalParam("decoyTolerance",1e-14,@(x) x>=0);
optionsParser.addOptionalParam("decoyPrecision","high");
optionsParser.addOptionalParam("decoyForceSep",false, @islogical);
optionsParser.addOptionalParam("decoySolver","MOSEK");;
optionsParser.addOptionalParam("photonCutOff",3, @(x) mustBeInteger(x));


optionsParser.parse(options);
options = optionsParser.Results;


%% modParser
modParser = moduleParser(mfilename);

%modParser.addAdditionalConstraint(@isEqualSize,["observablesJoint","expectationsJoint"]);

modParser.addRequiredParam("krausOps", @isCPTNIKrausOps);
modParser.addRequiredParam("keyProj", @(x) mustBeAKeyProj(x));

modParser.addRequiredParam("dimA",@mustBeInteger);
modParser.addRequiredParam("dimB", @mustBeInteger);
modParser.addAdditionalConstraint(@mustBePositive,"dimA")
modParser.addAdditionalConstraint(@mustBePositive,"dimB")
modParser.addRequiredParam("dimR",@mustBeInteger);

modParser.addRequiredParam("observablesJoint",@(x) allCells(x,@ishermitian));

modParser.addAdditionalConstraint(@observablesAndDimensionsMustBeTheSame,["observablesJoint","dimA","dimB"])

modParser.addRequiredParam("announcementsA")
modParser.addRequiredParam("announcementsB")
modParser.addRequiredParam("keyMap",@(x)mustBeA(x,"KeyMapElement"))
%modParser.addAdditionalConstraint(@mustBeSizedLikeAnnouncements,["expectationsJoint","announcementsA","announcementsB"])

%Epsilon(s) and postselection lift dimensions
modParser.addRequiredParam("log2targetepsSec");
modParser.addRequiredParam("log2targetepsCor");
modParser.addRequiredParam("epsBiasFactor", @(x) mustBeInRange(x, 0, 1));
modParser.addRequiredParam("liftDimension");
modParser.addOptionalParam("liftType","IID");

%Error correction
modParser.addRequiredParam("f", @(x) mustBeGreaterThanOrEqual(x,1));

% Finite-size
modParser.addRequiredParam("N", @(N) mustBeInteger(N)); % SignalCount

modParser.addRequiredParam("probsB", @(p) mustBeProbDist(p));
modParser.addRequiredParam("pTest", @(p) mustBeInRange(p, 0, 1)); 
modParser.addOptionalParam("rhoA", nan, @(x) all(isnan(x),"all") ||isDensityOperator(x));


modParser.addOptionalParam("blockDimsA", nan);
modParser.addOptionalParam("blockDimsB", nan);
modParser.addAdditionalConstraint(@(x,y) blockDimsMustMatch(x,y),["blockDimsA","dimA"]);
modParser.addAdditionalConstraint(@(x,y) blockDimsMustMatch(x,y),["blockDimsB","dimB"]);
modParser.addAdditionalConstraint(@(blockDimsA,blockDimsB) ~xor(all(isnan(blockDimsA),"all"),all(isnan(blockDimsB),"all")),["blockDimsA","blockDimsB"]);


modParser.addRequiredParam("expectationsConditional",@(x) eachRowMustBeAProbDist(x)); %%P( b | a, mu_i, test *or* gen, coherent state sent into channel)
modParser.addRequiredParam("detectorMat", @(detectorMat) all(size(detectorMat) == [6,2])); %Matrix describing Bob's detection setup

%Signal intensity
modParser.addRequiredParam("decoysSignal",@(x) allCells(x,@(y) y>=0));
modParser.addAdditionalConstraint(@(x) allCells(x,@isscalar),"decoysSignal");
modParser.addAdditionalConstraint(@(x) allCells(x,@(y) y>=0),"decoysSignal");
%Decoy intensity 1
modParser.addRequiredParam("decoys1",@(x) allCells(x,@(y) y>=0));
modParser.addAdditionalConstraint(@(x) allCells(x,@isscalar),"decoys1");
modParser.addAdditionalConstraint(@(x) allCells(x,@(y) y>=0),"decoys1");
% %Decoy intensity 2
% modParser.addRequiredParam("decoys2",@(x) allCells(x,@(y) y>=0));
% modParser.addAdditionalConstraint(@(x) allCells(x,@isscalar),"decoys2");
% modParser.addAdditionalConstraint(@(x) allCells(x,@(y) y>=0),"decoys2");

%Select decoy methods
modParser.addOptionalParam("decoyMethods","LP");
modParser.addOptionalParam("FlagSquash",true, @(x) islogical(x));

modParser.addRequiredParam("decoyProbsTest",@mustBeProbDistCell); %,must sum to 1.
modParser.addRequiredParam("decoyProbsGen",@mustBeProbDistCell); %,must sum to 1. 
modParser.addRequiredParam("probSignalsAgen",@mustBeProbDist); %Add something for same number of announcementsA or something
modParser.addRequiredParam("probSignalsAtest",@mustBeProbDist); %Add something for same number of announcementsA or something
modParser.addRequiredParam("signalsAlice", @(x) length(x) == 4); 
modParser.addRequiredParam("POVMB", @(x) length(x) == 8); 
modParser.addRequiredParam("pTestandACondSingle");
modParser.addRequiredParam("pMuandACondTest");
modParser.addRequiredParam("pSinglePhoton", @(p) p >= 0);
modParser.addRequiredParam("pGenCondSingle", @(p) p >= 0); 
modParser.addRequiredParam("pACondSingle", @mustBeProbDist); 
modParser.addRequiredParam("patternOrder"); 
modParser.addRequiredParam("pTestCondSingleandA");
modParser.addRequiredParam("pSingleCondGen");

modParser.parse(params);
params = modParser.Results;




%% simple setup
debugMathSolver = debugInfo.addLeaves("mathSolver");
mathSolverInput = struct();



%% Squashing

squashingMap = fourSixPostProcessingSquashingMap();
squashedConExp = pagemtimes(params.expectationsConditional,squashingMap.');

%Add another postprocessing required for WFSS
matWFSS = WFSSPosprocessing(squashedConExp);
squashedConExp = pagemtimes(squashedConExp,matWFSS);


pTest = params.pTest;
pGen = 1- pTest; 

%decoy probs p(mu|test)
decoyProbsTest = cell2mat(params.decoyProbsTest); 

%p(a|test)
probSignalsTesting = params.probSignalsAtest;

pMuACondTest = params.pMuandACondTest;
squashedExpCondionedOnTest = zeros(size(squashedConExp(:, :, 1))); %initialize for storage
%start P(b | a, mu, test) NO PHOTON CONDITIONING
for i=1:numel(decoyProbsTest)
    squashedExpCondionedOnTest(:,:,i) = diag(pMuACondTest(:, i))*squashedConExp(:,:,i);% P(a,b,mu | test)
end
squashedJointExp = (1-pGen)*squashedExpCondionedOnTest; %P(a,b,mu,test)

%% Error correction

%p(a|gen)
probSignalsGeneration = params.probSignalsAgen; 

%decoy probs p(mu|gen)
decoyProbsGen = cell2mat(params.decoyProbsGen);

%P(b|a,gen)
squashedExpCondAandGen = sum(bsxfun(@times,squashedConExp,reshape(decoyProbsGen,1,1,[])),3);

%P(a,b|gen)
squashedExpCondGen = diag(probSignalsGeneration)*squashedExpCondAandGen;

genAnnouncementsA = params.announcementsA;
squashedExpCondGen = squashedExpCondGen(:,:); 

[deltaLeak, gains] = errorCorrectionCost(genAnnouncementsA,params.announcementsB,...
    squashedExpCondGen,params.keyMap,params.f); %%squashedExpCondGen is conditioned on pGen 

debugInfo.storeInfo("deltaLeak",deltaLeak);
debugInfo.storeInfo("gains", gains);

%% translate for the math solver
%now we have all the parts we need to get a key rate from the a math
%solver, but we need to put it into a form it can understand.
%first we give it the kraus operators for the G map and the projection
%operators for the key map (Z).
mathSolverInput.krausOps = params.krausOps;
mathSolverInput.keyProj = params.keyProj;
% also include rhoA from the description if it was given

if ~isnan(params.rhoA)
    mathSolverInput.rhoA = params.rhoA;
end

%% PostSelection Technique usage

log2targetepsSec = params.log2targetepsSec;
log2targetepsCor = params.log2targetepsCor;

liftDimension = params.liftDimension;

switch params.liftType
    case 'IID'
        %IID version
        logCost = 0;
        params.logCost = logCost;
        logepsilon.EC = log2targetepsCor; % failure probability for error-correction
        logepsilon.AT = log2targetepsSec +log2(1-params.epsBiasFactor); % failure probability for acceptance test
        logepsilon.PA = log2targetepsSec+ log2(params.epsBiasFactor); % failure probability for privacy amplification
        keyCost = 0;
    case 'PS'
        %Full postselection technique
        logCostOld = (liftDimension-1)*log2(params.N+1);
        logCostNew = (liftDimension-1)*(log2(exp(1))+log2(params.N)+(liftDimension-1)/(params.N*log(2)) - log2(liftDimension-1));
        if (logCostNew > logCostOld)
            warning("The new log cost is HIGHER \n");
        end
        logCost = min(logCostNew,logCostOld);
        params.logCost = logCost;
        
        logepsilon.EC = log2targetepsCor; %failure probability for error-correction
        
        %Assuming \sqrt{8 \epsSec} = \epsilontilde / 2 = targetEpsSec /
        %2g_{n,x} as stated in the paper.
                                    
        logepsilon.PA = 2*(log2targetepsSec - logCost) + log2(params.epsBiasFactor) -5;
        logepsilon.AT = 2*(log2targetepsSec - logCost) + log2(1-params.epsBiasFactor) -5;
        
        keyCost = 4*logCost - 2*log2targetepsSec; %Assuming \sqrt{8 \epsSec} = \epsilontilde / 2 = targetEpsSec / 2g_{n,x}
    otherwise
        disp("Lift type does not match!")
end

%store key cost and resulting epsilons
params.keyCost = keyCost;
params.logepsilon = logepsilon;

%% Calculating Kappas / Confidence intervals (Size of set to optimize over).

squashedJointExpTestOnly = squashedJointExp; % We use entire data since each round can be tested with some probability

SigmaSize = numel(squashedJointExpTestOnly) + 1; % The 
% the 1 is needed for bounds on key and generation. 
params.SigmaSize = SigmaSize;

% we construct muball for entrywise using clopper peason confidence
% intervals.
debugInfo.storeInfo("SigmaSize",SigmaSize);
[TestingLowerBound, TestingUpperBound] = ConfidenceIntervalBoundsLogEpsV2(squashedJointExpTestOnly ,params.N, 1/log2(10)*logepsilon.AT - log10(SigmaSize)); %The last entry has to be log10(epsAT / SigmaSize);

muLower = max(squashedJointExpTestOnly - TestingLowerBound , 0);
muUpper = max(TestingUpperBound - squashedJointExpTestOnly , 0); %sometimes mu is a tiny negative value, which we overwrite to zero.

debugInfo.storeInfo("muLower",muLower);
debugInfo.storeInfo("muUpper",muUpper);


fprintf("\n Finite Size concentration inequality calculations are done \n ");

%% Decoy Analysis %%
%Perform detector decomposition to find bound on subspace
[~,cn,~,~] = detectorDecompFlag(params.detectorMat,1,true);

%2 decoy
reshapedIntensities = {cell2mat(params.decoysSignal);cell2mat(params.decoys1)};
decoyIntensityFlag = reshapedIntensities{2};

% %3 decoy
% reshapedIntensities = {cell2mat(params.decoysSignal);cell2mat(params.decoys1);cell2mat(params.decoys2)};
% decoyIntensityFlag = reshapedIntensities{2};


%Generate photon bounds for decoy analysis
%observations for flag-state squasher used in Decoy SDP
obsFlagDecoy = 1 - sum(  squashedConExp(:,1:6,:)  - muLower(:,1:6,:),2) - squashedConExp(:,end,:)  + muLower(:,end,:); %%make sure you understand this. 
for indexA = 1:numel(params.signalsAlice)  
    photonBoundDecoy(:,2) = 1 - obsFlagDecoy(:,:,2)./(decoyIntensityFlag(indexA)*exp(-decoyIntensityFlag(indexA))*cn); % 1-photon
    photonBoundDecoy(:,1) = 1 - obsFlagDecoy(:,:,2)./(exp(-decoyIntensityFlag(indexA))*cn); % 0-photon
end


tau = zeros(size(muLower)); % These don't play a role, but we have set it to zero because we are using 
% code written for fixed-length implementations, where tau / t was the size of
% acceptance test. 

if params.decoyMethods == "SDP"
    % SDP
    [CondExpectationsL, CondExpectationsU] = decoySDPFiniteDifInt(squashedJointExpTestOnly,...
        params.signalsAlice,params.POVMB,reshapedIntensities',probSignalsTesting,...
        decoyProbsTest,1-pGen,tau,muLower,muUpper, photonBoundDecoy, params.FlagSquash,...
       "decoySolver",options.decoySolver,"decoyTolerance",options.decoyTolerance,...
        "decoyPrecision",options.decoyPrecision,"verboseLevel",options.verboseLevel,"photonCutOff",options.photonCutOff);
elseif params.decoyMethods == "LP"
    % LP
    [CondExpectationsL, CondExpectationsU] = decoyLPFiniteDifInt(squashedJointExpTestOnly,...
        reshapedIntensities',probSignalsTesting,decoyProbsTest,1-pGen,tau,muLower,muUpper,...
       "decoySolver",options.decoySolver,"decoyTolerance",options.decoyTolerance,...
        "decoyPrecision",options.decoyPrecision,"verboseLevel",options.verboseLevel,"photonCutOff",options.photonCutOff);
else
    disp("Decoy methods don't match!")
end

%P(a,b,t| n) = P(a,t| n) P(b|a,t,n) , cond
JointExpectationsL = diag(params.pTestandACondSingle)*CondExpectationsL;  %P(a,t|n)P(b|a,t,n)=P(a,b,t|n)
JointExpectationsU = diag(params.pTestandACondSingle)*CondExpectationsU;


%Now, we add constraints. Recall that the POVMs correspond to Bob and Alice
%given Test . (and trivially intensity, and single photon)
pTestCondSingleandA = params.pTestCondSingleandA; 

observablesJointTesting = params.observablesJoint;

%View testing as post-processing
for a = 1:size(squashedConExp,1)
    for b = 1:size(squashedConExp,2)
        observablesJointTesting{a, b} = pTestCondSingleandA(a)*observablesJointTesting{a,b};
    end
end


observablesTest = observablesJointTesting(:,:);

numObs = numel(observablesTest);

% Add the constraint
decoyConstraints = arrayfun(@(index)InequalityConstraint(...
    observablesTest{index},JointExpectationsL(index),...
    JointExpectationsU(index)), 1:numObs);

%% Generate photon number constraints for a = 1, 2, 3,4

pU = params.pACondSingle(:); 
pL = zeros([4,1]); 


for i = 1:4
    %CondExpectationsL  S_(H V D A R L) D_( HV DA RL ) CC_(ANY) NON
    singleClicksM = sum(CondExpectationsL(i, 1:6)); 
    vacM = CondExpectationsL(i, end); 
    multiClicksM = sum(CondExpectationsU(i, 7)); 
    mObs = min(1 - singleClicksM - vacM, multiClicksM); 
    
    pL(i) = params.pACondSingle(i) * (1 - mObs/cn); 
end
%Constraint ops to project onto n <= NB subspace
projectionOps = cell([1, 4]);
Pi = blkdiag(eye(3), zeros(params.dimB-3)); 
for i = 1:4
    aa = zket(4, i)*zket(4, i)';
    projectionOps{i} = kron(aa, Pi); 
end
photonNumberConstraints = arrayfun(@(index) InequalityConstraint(projectionOps{index}, pL(index), pU(index)), 1:4); 

%Add all constraints to mathsolver
ineqConstraintsTotal = [decoyConstraints, photonNumberConstraints]; 
mathSolverInput.inequalityConstraints = ineqConstraintsTotal; 


% if block diag information was give, then pass it to the solver.
if ~any(isnan(params.blockDimsA),"all")
    mathSolverInput.blockDimsA = params.blockDimsA;
    mathSolverInput.blockDimsB = params.blockDimsB;
end


% now we call the math solver function on the formulated inputs, with
% it's options.
[relEnt,~] = mathSolverFunc(mathSolverInput,debugMathSolver);

%store the key rate (even if negative)
keyRate = computeFiniteAdaptiveKeyRate(relEnt, deltaLeak, gains, params); 

if options.verboseLevel>=1
    %ensure that we cut off at 0 when we display this for the user.
    fprintf("Key rate: %e\n",max(keyRate,0));
end

%set the linearization estimate key rate as well for debuging
if isfield(debugMathSolver.info,"relEntStep2Linearization")
    keyRateStep2Linearization = debugMathSolver.info.relEntStep2Linearization - deltaLeak ; %change our 
    debugInfo.storeInfo("keyRaterelEntStep2Linearization",keyRateStep2Linearization)

    if options.verboseLevel>=2
        fprintf("Key rate using step 2 linearization intial value: %e\n",max(keyRateStep2Linearization,0))
    end
end
end



function [finiteKeyRate] = computeFiniteAdaptiveKeyRate(relEnt, deltaLeak , gains, params)
    N = params.N; 
    pTest  = params.pTest; 
    pGen = 1-pTest; 

    
    dimR = params.dimR; 
    logepsilon = params.logepsilon;

    % pSinglePhotonCondGen = params.decoys{1}*exp(-params.decoys{1}); 
   % pSinglePhotonCondGen = params.pSingleCondGen;
   % pGenCondSingle = params.pGenCondSingle;
    ObservedNsift = pGen*sum(gains)*N; 
    
    ObservedSiftProb = pGen*sum(gains); 

    [~, siftUpperBound] = ConfidenceIntervalBoundsLogEpsV2(ObservedSiftProb,N, 1/log2(10)*logepsilon.AT - log10(params.SigmaSize)); %% DO SOMETHING HERE>.
   % debugInfo.storeInfo("siftLowerBound",siftLowerBound);
   % debugInfo.storeInfo("siftUpperBound",siftUpperBound);

    alpha = 1 + (1/(log2(dimR+1)) )*sqrt( -logepsilon.PA / ObservedNsift   ); %This is optimal alpha.
    %Technically this should be chosen BEFORE the protocol is run. For now
    %the above code assumes that we chose a "good" one in advance. 
    keyCost = params.keyCost / N;
    
    prefactor = params.pSinglePhoton*(ObservedSiftProb)/(siftUpperBound); % Note this is prob of Single photon PERIOD without any conditioning on test or gen ...
    ECCost = pGen * deltaLeak + (1-logepsilon.EC)/N; 
    CorrectionTermOne = ObservedNsift*(alpha-1)*(log2(dimR+1))^2 / N; %basically AEP Correction Term
    CorrectionTermTwo = (alpha / (alpha-1) * (-2 - logepsilon.PA+2/alpha)) / N; %PA correction term..   
    finiteKeyRate = prefactor*relEnt - ECCost - CorrectionTermOne - CorrectionTermTwo - keyCost; %check if this needs to be key cost vs params.keycost!!!!

end

%Checker that each row is a probability distribution
function eachRowMustBeAProbDist(expectationsConditional)

% get the dimensions of the conditional expectations. Then based on that
% pick a strategy to handle it
dimExpCon = size(expectationsConditional);

errorID ="FourSixWCPKeyRateFunc:InvalidRowsAreNotProbDists";
errorTXT = "A row in the conditional distribution is not a valid probability distribution.";

if numel(dimExpCon) == 2 % Matlab's minimum number of dimensions is 2.
    % The array is 2d and the slicing is easy
    for index = 1:dimExpCon(1)
        if~isProbDist(expectationsConditional(index,:))
           throwAsCaller(MException(errorID,errorTXT));
        end
    end
else
    % We have some tricky slicing to do for 3 plus dimensions.
    % We need to index the first dimension and the combination of
    % dimensions 3 and up. The second dimension will just use :.
    maskedDims = [dimExpCon(1),prod(dimExpCon(3:end))];

    for index = 1:prod(maskedDims)
        vecIndex = ind2subPlus(maskedDims,index);
        if ~isProbDist(expectationsConditional(vecIndex(1),:,vecIndex(2)))
            throwAsCaller(MException(errorID,errorTXT));
        end
    end
end
end

function observablesAndDimensionsMustBeTheSame(observables,dimA,dimB)
if ~allCells(observables,@(x) size(x,1) == dimA*dimB)
    throwAsCaller(MException("FiniteKeyRateFunc:ObservablesAndDimensionsMustBeTheSame","The Observables must have the same dimensions as Alice and Bob multiplied together."));
end
end
function mustBeProbDistCell(input)
mustBeProbDist([input{:}])
end
function mustBeSizedLikeAnnouncements(jointExpectations,announcementsA,announcementsB)
if ~isequal(size(jointExpectations),[numel(announcementsA),numel(announcementsB)])
    throwAsCaller(MException("FiniteKeyRateFunc:jointKeyDoesNotHaveSameSizeAsAnnouncements",...
        "The joint key distribution must have size numel(announcementsA) by numel(announcementsB)."))
end
end

function mat = WFSSPosprocessing(squashedConExp)
    [~,numB,~] = size(squashedConExp);
    mat = zeros(numB,8);
    %Single clicks stay the same
    for index = 1:6
        mat(:,index) = zket(numB,index);
    end
    %No clicks stay the same
    mat(:,8) = zket(numB,11);
    %Double and cross clicksget combined to one multiclick
    mat(7:10,7) = 1;
end