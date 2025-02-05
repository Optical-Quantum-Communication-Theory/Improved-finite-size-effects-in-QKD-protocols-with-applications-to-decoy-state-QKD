function [keyRate, modParser, debugInfo] = FiniteDecoyBB84KeyRateFunc(params,options,mathSolverFunc,debugInfo)
% FiniteDecoyBB84KeyRateFunc A finite size key rate function for a 
% lossy qubit BB84 protocol withloss. Uses a 3 dimensional Bob which has a
% vacuum dimension perpendicular to qubit dimensions. 
%
% Input parameters:
% * dimA: dimension of Alice's system.
% * dimB: dimension of Bob's system.
% * f: error correction effiency. A 1 means we are correcting at the
%   Shannon limit. A more practical value would be around 1.16.
% * observablesJoint: The joint observables from Alice and Bob's
%   measurments. These are organized in a 4x6 table. The
%   observables must be hermitian and each must be the size dimA*dimB by
%   dimA*dimB. The observables assume the spaces are ordered A \otimes B.
%   They also should be positive semi-definite and should sum to identity,
%   but this is hard to check because of machine precision issues.
% * N: (finite size) the total number of signals Alice sent
% * eps: (finite size) a struct containing four epsilons used in finite key
%   rate analysis: PE (parameter estimation), bar (smoothing min entropy),
%   EC (error correction), PA (privacy amplification)
% * ptest: the fraction of signals sent that are used for testing.
%   Typically should be fairly low. This cuts into key rate, but also
%   entirely destroys key rate if it is too low.
% * t: a finite size parameter that controls the looseness of the allowed
%   frequencies compared to expected ideal probabilities. This should be
%   optimized over in the presence of experimental data.
% * physDimAB: physical dimensions of Alice and Bob's measurements, used
%   for postselection
% * expectationsConditional : Conditional expectations of Bob's outcome
%   Alice sent a given signal at a given intensity. (this is also
%   conditioned on it being a test round).
% * ObservablesJoint : POVMs corresponding to Alice sending a given signal
%   and Bob obtaining a given outcome, conditioned on test and
%   single-photon being sent.
% Outputs:
% * keyrate: Key rate of the QKD protocol.
% Options:
% * verboseLevel: (global option) See makeGlobalOptionsParser for details.
% * errorHandling: (global option) See makeGlobalOptionsParser for details.
% DebugInfo:
% * krausSum: sum_i K^\dagger_i*K_i. For a completely positive trace
%   non-increasing map, this sum should be <=I. 
%
% See also QKDKeyRateModule, PM46DescriptionFunc, makeGlobalOptionsParser
arguments
    params (1,1) struct
    options (1,1) struct
    mathSolverFunc (1,1) function_handle
    debugInfo (1,1) DebugInfo
end

%% opitions parser
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.addOptionalParam("decoyTolerance",1e-14,@(x) x>=0);
optionsParser.addOptionalParam("decoySolver","SDPT3");
optionsParser.addOptionalParam("decoyPrecision","high");
optionsParser.parse(options);
options = optionsParser.Results;

%% module parser
modParser = moduleParser(mfilename);
modParser.addRequiredParam("observablesJoint",@(x) allCells(x,@ishermitian));
modParser.addRequiredParam("expectationsConditional",@(x) eachRowMustBeAProbDist(x));
modParser.addRequiredParam("decoys",@(x) allCells(x,@(y) y>=0));
modParser.addRequiredParam("decoyProbs",@mustBeProbDistCell); %,must sum to 1. 
% modParser.addAdditionalConstraint(@isEqualSize,["observablesJoint","expectationsConditional"]);

modParser.addRequiredParam("signalsAlice", @(x) length(x) == 4); 
modParser.addRequiredParam("POVMB", @(x) length(x) == 5); 
modParser.addRequiredParam("krausOps", @isCPTNIKrausOps);
modParser.addRequiredParam("keyProj", @(x) mustBeAKeyProj(x));
modParser.addRequiredParam("probSignalsA",@mustBeProbDist); %Add something for same number of announcementsA or something

%Select decoy methods
modParser.addOptionalParam("decoyMethods","LP",@(x)ismember(x,["SDP","LP"]));
modParser.addOptionalParam("FlagSquash",false, @(x) islogical(x));

modParser.addRequiredParam("dimA",@mustBeInteger);
modParser.addRequiredParam("dimB", @mustBeInteger);
modParser.addAdditionalConstraint(@mustBePositive,"dimA")
modParser.addAdditionalConstraint(@mustBePositive,"dimB")
modParser.addAdditionalConstraint(@observablesAndDimensionsMustBeTheSame,["observablesJoint","dimA","dimB"])

modParser.addRequiredParam("announcementsA")
modParser.addRequiredParam("announcementsB")
modParser.addRequiredParam("keyMap",@(x)mustBeA(x,"KeyMapElement"))
% modParser.addAdditionalConstraint(@mustBeSizedLikeAnnouncements,["expectationsConditional","announcementsA","announcementsB"])

modParser.addRequiredParam("f", @(x) mustBeGreaterThanOrEqual(x,1));
modParser.addRequiredParam("rhoA",@isDensityOperator)
%modParser.addOptionalParam("rhoA", nan, @isDensityOperator);
modParser.addRequiredParam("alphabet", @(x) mustBeInteger(x));


%% finite key analysis parameters
modParser.addRequiredParam("N", @(x) mustBeGreaterThan(x, 0));
modParser.addRequiredParam("ptest", @(x) mustBeInRange(x, 0, 1));
modParser.addRequiredParam("epsilonTotal", @(x) mustBeInRange(x, 0, 1));
modParser.addOptionalParam("epsilon", @(x) mustBeInRange(x, 0, 1));
modParser.addRequiredParam("tExp", @(x) mustBeLessThanOrEqual(x, 0));
modParser.addRequiredParam("tSiftExp", @(x) mustBeLessThanOrEqual(x,0));

% modParser.addOptionalParam("blockDimsA", nan);
% modParser.addOptionalParam("blockDimsB", nan);
% modParser.addAdditionalConstraint(@(x,y) blockDimsMustMatch(x,y),["blockDimsA","dimA"]);
% modParser.addAdditionalConstraint(@(x,y) blockDimsMustMatch(x,y),["blockDimsB","dimB"]);
% modParser.addAdditionalConstraint(@(blockDimsA,blockDimsB) ~xor(all(isnan(blockDimsA),"all"),all(isnan(blockDimsB),"all")),["blockDimsA","blockDimsB"]);

modParser.parse(params);

params = modParser.Results;

%% extra parameters
%modParser.addOptionalParam("acceptanceSetChoice", "none", @(x) mustBeMember(x, ["upper", "lower", "independent", "none"]));


%% parse parameters
modParser.parse(params);

params = modParser.Results;


debugMathSolver = debugInfo.addLeaves("mathSolver");
mathSolverInput = struct();


%%% we start computations %%%

%% Squashing map %%
% We use a simple squashed to squash Bob to Qubit + Vac. 
squashingMap = BB84StandardSquashingPostProccessingMap();
squashedConExp = pagemtimes(params.expectationsConditional,squashingMap.');

% The pagemtimes applies the squashing post-processing map on each "slice"
% of the conditional expectations (for each intensity). Note that since we
% condition on Alice's signal and intensity before squashing, we
% conditional expectations after squashing. 


%% Error correction
% We compute error-correction cost on the squashed probability distribution 
% for the signal intensity, conditioned on the generation rounds

squashedJointExpSignalOnly = diag(params.probSignalsA)*squashedConExp(:,:,1); %these are joint between Alice and Bob, conditioned on signal and generation
[deltaLeak, gains] = errorCorrectionCost(params.announcementsA,params.announcementsB, squashedJointExpSignalOnly,params.keyMap,params.f); 
debugInfo.storeInfo("deltaLeak",deltaLeak);
debugInfo.storeInfo("gains",gains);


%% Optimal epsilons 
%Retrieve t values
params.t = 10.^(params.tExp)*ones(size(squashedConExp));
params.tsift = 10.^(params.tSiftExp);

%Compute Optimal epsilons
epsilonSecurity = params.epsilonTotal;
[eps_PA, eps_EC, eps_PE, epsBar] = optimalEpsilon(epsilonSecurity,gains,params);

%Store optimal epsilon values
params.epsilon.PE = eps_PE; %failure probability for parameter estimation 
params.epsilon.bar = epsBar; %uncertainty from "smoothing" of min-entropy among states similar to feasible state
params.epsilon.EC = eps_EC; %failure probability for error-correction
params.epsilon.PA = eps_PA; %failure probability for privacy amplification


%% Compute mu .
%in order to compute mu, we need to first convert out conditional
%expectations into joint expectations on Alice and Bob and Intensity AND
%test.

decoyProbs = cell2mat(params.decoyProbs);
for i=numel(decoyProbs):-1:1
    temp = diag(params.probSignalsA)*squashedConExp(:,:,i);
    squashedJointExpCondionedOnTest(:,:,i)= decoyProbs(i)*temp;
end

squashedJointExp = params.ptest*squashedJointExpCondionedOnTest;

%then we compute mu as follows
tau = params.t.*ones(size(squashedJointExp));
[muLower,muUpper] = mubetaUpperLower(params.N, squashedJointExp, tau, log10(params.epsilon.PE));  



%% Perform the decoy analysis. 
%
decoyIntensities = cell2mat(params.decoys); %these expectations are bounds on prob (Bob | Alice, 1, key/test, intensity); 
reshapedIntensities = mat2cell(transpose(ones(numel(params.probSignalsA),numel(decoyIntensities)).*decoyIntensities),...
    ones(numel(decoyIntensities),1));
%no Flagstate squasher used => no photon number bounds
photonBoundDecoy = zeros(numel(params.probSignalsA),2);

switch params.decoyMethods
    case "SDP"
        % SDP
        [CondExpectationsL, CondExpectationsU] = decoySDPFiniteDifInt(squashedJointExp,...
        params.signalsAlice,params.POVMB,reshapedIntensities',params.probSignalsA,...
        decoyProbs,params.ptest,tau,muLower,muUpper,photonBoundDecoy, params.FlagSquash, ...
       "decoySolver",options.decoySolver,"decoyTolerance",options.decoyTolerance,...
        "decoyPrecision",options.decoyPrecision,"verboseLevel",options.verboseLevel);
    case "LP"
        % LP
        [CondExpectationsL, CondExpectationsU] = decoyLPFiniteDifInt(squashedJointExp, ...
        reshapedIntensities',params.probSignalsA,decoyProbs,params.ptest,tau,muLower,muUpper, ...
       "decoySolver",options.decoySolver,"decoyTolerance",options.decoyTolerance,...
        "decoyPrecision",options.decoyPrecision,"verboseLevel",options.verboseLevel);
    otherwise
        throw(MException("funcName:NotDecoyMethod","Decoy methods don't match!"))
end

% We convert them to Prob (Bob and Alice | one-photon, key/test, intensity);
JointExpectationsL = diag(params.probSignalsA)*CondExpectationsL;
JointExpectationsU = diag(params.probSignalsA)*CondExpectationsU;

%Now, we add constraints. Recall that the POVMs correspond to Bob and Alice
%given Test . (and trivially intensity, and single photon).

numObs = numel(params.observablesJoint);

%Add constraints to math solver
mathSolverInput.inequalityConstraints = arrayfun(@(index)InequalityConstraint(...
    params.observablesJoint{index},JointExpectationsL(index),...
    JointExpectationsU(index)), 1:numObs);


%% Translate for math solver
mathSolverInput.krausOps = params.krausOps;
mathSolverInput.keyProj = params.keyProj;
mathSolverInput.rhoA = params.rhoA;


%Calculate relative entropy from math solver
[relEnt,~] = mathSolverFunc(mathSolverInput, debugMathSolver);


%Calculate relative entropy from math solver
keyRate = finiteKeyRate(relEnt, deltaLeak, gains, params, options);


%Store debug info
if isfield(debugMathSolver.info,"relEntStep2Linearization")
    relEntStep2Linearization = debugMathSolver.info.relEntStep2Linearization; 
    
    keyRateStep2Linearization = finiteKeyRate( relEntStep2Linearization, deltaLeak, gains, params, options);
    if options.verboseLevel>=2
        fprintf("Key rate using step 2 linearization intial value: %e\n",max(keyRateStep2Linearization,0))
    end
    debugInfo.storeInfo("keyRateRelEntStep2Linearization",keyRateStep2Linearization)
    
end


%Print final key rate
if options.verboseLevel>=1
    fprintf("Key rate: %e\n",keyRate);
end
    
end



%%%%%%%%%%%  FINITE Key Rate Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Finite-size key rate function from Lars Kamin

function [keyRate, debugInfo] = finiteKeyRate(relEnt, deltaLeak, gains, params, options, debugInfo)
%computes the finite size keyrate.

%Alice's key alphabet size
dA = params.alphabet;

%testing probability
ptest = params.ptest;

%key generation probability
pgen = 1 - ptest;

%security parameters
epsilon = params.epsilon;

%total signals sent
N = params.N;

%t parameter for allowed fluctuations
t = params.t;
tsift = params.tsift;

%AEP correction in purified distance
AEPCorrectionPD = 2*log2(2*dA+1)*sqrt(1-2*log2(epsilon.bar));

%Privacy Amplification
privacyAmplification = 2*log2(1/2/epsilon.PA)/N;

%Error correction including error verification
% deltaleak calculated conditioned on having a generation round, hence the
% prefactor of pgen
ECLeakage = pgen*deltaLeak + log2(2/epsilon.EC)/N; %actually epsilon.EV this is!

%Total gain
totalGain = sum(gains,1);

%Sifted number of signals
n_sift = pgen*totalGain*N; 

%Mu sift
[~,musift] = mubetaUpperLower(N,totalGain*pgen,tsift,log10(epsilon.PE));

%Probabilty of Alice sending a single photon 
pSinglePhoton = params.decoys{1}*exp(-params.decoys{1});

%Calculate prefactor of relative entropy
prefactor = pgen*pSinglePhoton/N*floor(n_sift-tsift*N)/(totalGain*pgen + tsift + musift); 


%Final resulting finite-size key rate
keyRate = prefactor*relEnt - sqrt((totalGain*pgen - tsift)/N)*AEPCorrectionPD - ECLeakage - privacyAmplification; 

end


function observablesAndDimensionsMustBeTheSame(observables,dimA,dimB)
if ~allCells(observables,@(x) size(x,1) == dimA*dimB)
    throwAsCaller(MException("BasicKeyRateFunc:ObservablesAndDimensionsMustBeTheSame","The Observables must have the same dimensions as Alice and Bob multiplied together."));
end
end

function mustBeSizedLikeAnnouncements(jointExpectations,announcementsA,announcementsB)
if ~isequal(size(jointExpectations),[numel(announcementsA),numel(announcementsB)])
    throwAsCaller(MException("BasicKeyRateFunc:jointKeyDoesNotHaveSameSizeAsAnnouncements",...
        "The joint key distribution must have size numel(announcementsA) by numel(announcementsB)."))
end
end

function mustBeProbDistCell(input)
mustBeProbDist([input{:}])
end


%%%%%%%%%%%  Squashing map %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mapping = BB84StandardSquashingPostProccessingMap()
% squashes detector patterns from H,V,D,A (as bit string patterns 0000,
% 1000, ..., 1111) to qubit+vac values H,V,D,A,vac.

    function mapping = quickMap(mapping,pattern,remaping)
        mapping(:,sub2indPlus(2*ones(1,numel(pattern)),pattern+1)) = remaping;
    end

mapping = zeros(5,16);

% The vast majority of the squashed bits are cross clicks that are mapped
% to vac for discarding. We will replace patterns that don't represent
% cross clicks in later steps.
mapping(5,:) = 1;

% vacume to vacume
mapping = quickMap(mapping,[0,0,0,0],[0,0,0,0,1]);

% single clicks to single clicks
mapping = quickMap(mapping,[1,0,0,0],[1,0,0,0,0]); % H
mapping = quickMap(mapping,[0,1,0,0],[0,1,0,0,0]); % V
mapping = quickMap(mapping,[0,0,1,0],[0,0,1,0,0]); % D
mapping = quickMap(mapping,[0,0,0,1],[0,0,0,1,0]); % A

% double clicks
mapping = quickMap(mapping,[1,1,0,0],[0.5,0.5,0,0,0]); % Z (HV)
mapping = quickMap(mapping,[0,0,1,1],[0,0,0.5,0.5,0]); % X (DA)
end


%Checker that each row is a probability distribution
function eachRowMustBeAProbDist(expectationsConditional)

% get the dimensions of the conditional expectations. Then based on that
% pick a strategy to handle it
dimExpCon = size(expectationsConditional);

errorID ="BasicBB84_WCPKeyRateFunc:InvalidRowsAreNotProbDists";
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
