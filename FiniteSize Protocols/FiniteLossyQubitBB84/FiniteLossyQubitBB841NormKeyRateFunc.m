function [keyRate, modParser, debugInfo] = FiniteLossyQubitBB841NormKeyRateFunc(params,options,mathSolverFunc,debugInfo)
% FiniteLossyQubitBB84KeyRateFunc A finite size key rate function for a 
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
% * expectationsJoint: The joint expectations from Alice and Bob's 
%   measurements. These should line up with the corresponding observables 
%   at each entry.
% 
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

optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;


modParser = moduleParser(mfilename);

modParser.addRequiredParam("observablesJoint",@(x) allCells(x,@ishermitian)); %these are observations conditioned on test.
modParser.addRequiredParam("expectationsJoint", @(x) all(or(x>=0,x<=1),"all"));  %these are POVMs conditioned on test. 

% modParser.addRequiredParam("pDetection", @(x) x>=0);
modParser.addRequiredParam("f", @(x) x>=1);
modParser.addRequiredParam("announcementsA")
modParser.addRequiredParam("announcementsB")
modParser.addRequiredParam("keyMap",@(x)mustBeA(x,"KeyMapElement"))
modParser.addRequiredParam("keyProj", @(x) mustBeAKeyProj(x));

modParser.addAdditionalConstraint(@mustBeSizedLikeAnnouncements,["expectationsJoint","announcementsA","announcementsB"])
modParser.addOptionalParam("fullstat", 1, @(x) mustBeMember(x, [-1,0,1]));
modParser.addOptionalParam("rhoA", nan, @(x) abs(trace(x)-1)<=eps);


modParser.addOptionalParam("physDimAB", 2*3, @(x) mustBeInteger(x));
modParser.addRequiredParam("alphabet", @(x) mustBeInteger(x));

modParser.addRequiredParam("dimA",@(x) x==2);
modParser.addRequiredParam("dimB", @(x) x ==3);

modParser.addRequiredParam("krausOps",@(x) isCPTNIKrausOps(x));

%% finite key analysis parameters
modParser.addRequiredParam("N", @(x) mustBeGreaterThan(x, 0));
modParser.addRequiredParam("m", @(x) mustBeGreaterThan(x, 0));
modParser.addRequiredParam("Sigma", @(x) mustBeGreaterThan(x, 0));
modParser.addOptionalParam("epsilon", @(x) mustBeInRange(x, 0, 1));
modParser.addRequiredParam("tExp", @(x) mustBeLessThanOrEqual(x, 0));

%% extra parameters
%modParser.addOptionalParam("acceptanceSetChoice", "none", @(x) mustBeMember(x, ["upper", "lower", "independent", "none"]));


%% parse parameters
modParser.parse(params);

params = modParser.Results;


debugMathSolver = debugInfo.addLeaves("mathSolver");
mathSolverInput = struct();


%%% we start computations %%%

%% compute EC cost conditioned on the generation rounds
[deltaLeak, gains] = errorCorrectionCost(params.announcementsA,params.announcementsB, params.expectationsJoint,params.keyMap,params.f);
%note that this deltaLeak corresponds to us correcting errors only in Key
%rounds. 
debugInfo.storeInfo("deltaLeak",deltaLeak);
debugInfo.storeInfo("gains",gains);

%Retrieve t value
params.t = 10.^(params.tExp);

%% compute mu
muBall = mu1Norm(params.epsilon,params.Sigma,params.m);  
%multiplying by ptest gives us Prob(Alice event and Bob event AND test).
params.muBall = muBall;



%%% translate for the math solver
mathSolverInput.krausOps = params.krausOps;
mathSolverInput.keyProj = params.keyProj;
mathSolverInput.rhoA = params.rhoA;

% Create cell array of POVM elements and vecor for 1-norm constraint
numObs = numel(params.observablesJoint);

cellObs = cell(numObs,1);
vecObs = zeros(numObs,1);

for index=1:numObs
    cellObs{index} = params.observablesJoint{index};
    vecObs(index) = params.expectationsJoint(index);
end


%Add constraints to math solver
mathSolverInput.vectorOneNormConstraints = VectorOneNormConstraint(cellObs,vecObs,params.muBall + params.t);

%Calculate relative entropy from math solver
[relEnt,~] = mathSolverFunc(mathSolverInput, debugMathSolver);

%Calculate finite-size key rate
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
%Finite-size key rate function from Ian George et al.

function [keyRate, debugInfo] = finiteKeyRate(relEnt, deltaLeak, gains, params, options, debugInfo)
%computes the finite size keyrate.

%Alice's key alphabet size
dA = params.alphabet;

%total signals sent
N = params.N;

%testing fraction of signals
m = params.m;

%key generation rounds
n = N - m;

%security parameters
epsilon = params.epsilon;

%AEP correction in purified distance
AEPCorrectionTD = 2*log2(dA+3)*sqrt(log2(2/epsilon.bar)/n);

%Privacy Amplification
privacyAmplification = 2*log2(2/epsilon.PA);

%Error correction including error verification
% deltaleak calculated conditioned on having a generation round, hence the
% prefactor of pgen
ECLeakage = n*deltaLeak + log2(2/epsilon.EC); %actually epsilon.EV this is!

%Final resulting finite-size key length
keyLength = n*(max(relEnt,0) - AEPCorrectionTD) - ECLeakage - privacyAmplification; 

keyRate = keyLength/N;
end


%Checker that all observables have the same dimensions
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
