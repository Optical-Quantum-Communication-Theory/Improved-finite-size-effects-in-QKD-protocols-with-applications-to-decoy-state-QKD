function [newParams,modParser] = FourSixWCPDifInt_Description(params, options, debugInfo)
% BasicBB84_LossyDescriptionFunc A simple description function for a qubit BB84
% protocol with loss, using the Schmidt decomposition to turn Alice's
% 4d space of signals sent to a 2d space. 
%
% Input parameters:
% * pzA/B: The probability that Alice and Bob measure in the Z-basis. 
% Output parameters:
% * observablesJoint: The joint observables for Alice and Bob's measurement
%   of the signals.
% * dimA: dimension of Alice's system.
% * dimB: dimension of Bob's system.
% * probSignalsA: Probability of Alice selecting a signal to send to Bob.
%   In this protocol, it is half the probability of the basis choice.
% * rhoA: Alice's reduced density matrix for prepare-and-measure based
%   protocols.
% * POVMA: Alice's set of POVM operators which she measures her state with
%   in the source replacement scheme.
% * POVMB: Bob's set of POVM operators which he measures his state with.
% * announcementsA: Alice's announcements for each of her POVM operators.
%   Can be integers or strings.
% * announcementsB: Bob's announcements for each of his POVM operators.
%   Can be integers or strings.
% * keyMap: An array of KeyMap objects that contain pairs of accepted
%   announcements and an array dictating the mapping of Alice's measurement
%   outcome to key bits (May be written with Strings).
% * krausOps: A cell array of matrices. The Kraus operators that form the G
%   map on Alice and Bob's joint system. These should form a completely
%   postive trace non-increasing linear map. Each Kraus operator must be
%   the same size.
% * keyProj:  A cell array of projection operators that extract the key
%   from G(\rho). These projection operators should sum to identity. This
%   map is often called Z.
% Options:
% * none
% DebugInfo:
% * krausSum: sum_i K^\dagger_i*K_i which should be <= identity for
%   a CPTNI map.
%
% See also QKDDescriptionModule, makeGlobalOptionsParser
arguments
    params (1,1) struct
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end

%% options parser
%Parsing technical options for the module
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;

%% module parser
%Parsing parameters for the module
modParser = moduleParser(mfilename);
% modParser.addRequiredParam("pGen",@(x) mustBeInRange(x,0,1));
modParser.addRequiredParam("pTest",@(x) mustBeInRange(x,0,1));
modParser.addRequiredParam("probsB",@(x) mustBeProbDist(x));
modParser.addRequiredParam("N", @(N) mustBeInteger(N));
modParser.addRequiredParam("epsSound", @(e) mustBeNonnegative(e));
modParser.addRequiredParam("epsATfrac", @(f) mustBeInRange(f, 0, 1));
modParser.addRequiredParam("t", @(t) mustBeNonnegative(t)); 

%Signal intensity
modParser.addRequiredParam("decoysSignal",@(x) allCells(x,@(y) y>=0));
modParser.addAdditionalConstraint(@(x) allCells(x,@isscalar),"decoysSignal");
modParser.addAdditionalConstraint(@(x) allCells(x,@(y) y>=0),"decoysSignal");
%Decoy intensity 1
modParser.addRequiredParam("decoys1",@(x) allCells(x,@(y) y>=0));
modParser.addAdditionalConstraint(@(x) allCells(x,@isscalar),"decoys1");
modParser.addAdditionalConstraint(@(x) allCells(x,@(y) y>=0),"decoys1");
%Decoy intensity 2
modParser.addRequiredParam("decoys2",@(x) allCells(x,@(y) y>=0));
modParser.addAdditionalConstraint(@(x) allCells(x,@isscalar),"decoys2");
modParser.addAdditionalConstraint(@(x) allCells(x,@(y) y>=0),"decoys2");

modParser.addRequiredParam("decoyProbsTest",@mustBeProbDistCell); %,must sum to 1.
modParser.addRequiredParam("decoyProbsGen",@mustBeProbDistCell); %,must sum to 1. 
% modParser.addRequiredParam("detectorMat", @(detectorMat) size(detectorMat) == [2,6]); %Matrix describing Bob's detection setup

modParser.parse(params)
params = modParser.Results;

%% simple setup
%Setup parameters for later use
newParams = struct();

pTest = params.pTest; 
pGen = 1 - pTest; 

%Define + and - states for later use
ketH = [1;0];
ketV = [0;1];
ketD = [1;1]/sqrt(2);
ketA = [1;-1]/sqrt(2);
ketR = [1;1i]/sqrt(2);
ketL = [1;-1i]/sqrt(2);

dimA = 4; % |1,2,3,4>
dimB = 10; % Qubit+Vac + 6 single-click flags + Mult = 3 + 6 + 1 = 10 

newParams.dimA = dimA;
newParams.dimB = dimB;
probsB = params.probsB;
pzB = probsB(1);
pxB = probsB(2);
pyB = probsB(3); 

%% generate rhoA
% probSignalsAgen = [1/2, 1/2, 0, 0]; %p(a|gen)
% probSignalsAtest = [0, 0, 1/2, 1/2]; %p(a|test)

probSignalsAgen = [1/4, 1/4, 1/4, 1/4]; %p(a|gen)
probSignalsAtest = [1/4, 1/4, 1/4, 1/4]; %p(a|test)

newParams.probSignalsAgen = probSignalsAgen; 
newParams.probSignalsAtest = probSignalsAtest; 

%Store intensities in a matrix
% %2 decoy
% reshapedIntensities = [cell2mat(params.decoysSignal);cell2mat(params.decoys1)];

%3 decoy
reshapedIntensities = [cell2mat(params.decoysSignal);cell2mat(params.decoys1);cell2mat(params.decoys2)];

%Calculate probabilities
[pA, pTestCondSingle, pGenCondSingle, pTestCondSingleandA, pGenCondSingleandA, pACondTestandSingle, pACondSingle, pSinglePhoton ,pSingleCondGen, pMuandACondTest] = ...
    computeProbabilities46DifInt(pGen, pTest, probSignalsAgen, probSignalsAtest, reshapedIntensities, params.decoyProbsTest,params.decoyProbsGen);


newParams.pSinglePhoton = pSinglePhoton; 
debugInfo.storeInfo("pSinglePhoton", pSinglePhoton); 
newParams.pTestandACondSingle = pTestCondSingle*pACondTestandSingle; %p(a,test|n=1)
newParams.pGenCondSingle = pGenCondSingle; %p(gen|n=1)
newParams.pACondSingle = pACondSingle; %p(a|n=1)
newParams.pTestCondSingleandA = pTestCondSingleandA; %p(test|n=1,a) for a =1,2,3
newParams.pMuandACondTest = pMuandACondTest; %p(mu,a|test)
debugInfo.storeInfo("pMuandACondTest", pMuandACondTest);
newParams.pSingleCondGen = pSingleCondGen;

%Alice's states conditioned on single photons
psi1 = sqrt(pACondSingle(1))*kron(zket(4,1), ketH);
psi2 = sqrt(pACondSingle(2))*kron(zket(4,2), ketV);
psi3 = sqrt(pACondSingle(3))*kron(zket(4,3), ketD);
psi4 = sqrt(pACondSingle(4))*kron(zket(4,4), ketA);


psiAAp = psi1 + psi2 + psi3 + psi4;  %% |psi> conditioned on n = 1. 
isUnitVector(psiAAp); 

rhoAAp = psiAAp * psiAAp'; 
newParams.rhoAAp = rhoAAp; 

rhoA = PartialTrace(rhoAAp, 2, [4, 2]); 
newParams.rhoA = rhoA; 

%store Alice's sent states
newParams.signalsAlice = {ketH,ketV,ketD,ketA};

%% joint obserables
POVMsA = {diag([1,0,0,0]), diag([0,1,0,0]), diag([0,0,1,0]), diag([0,0,0,1])}; 

%Flag state POVMsB
%order S_( H V D A R L ) Multi-click NON
ketD = [ketD; 0]; 
ketA = [ketA; 0];
ketR = [ketR; 0]; 
ketL = [ketL; 0];
GammaH = blkdiag(pzB*diag([1,0,0]), diag(zket(7,1)));
GammaV = blkdiag(pzB*diag([0,1,0]), diag(zket(7,2)));
GammaD = blkdiag(pxB*(ketD*(ketD')), diag(zket(7,3)));
GammaA = blkdiag(pxB*(ketA*(ketA')), diag(zket(7,4)));
GammaR = blkdiag(pyB*(ketR*(ketR')), diag(zket(7,5)));
GammaL = blkdiag(pyB*(ketL*(ketL')), diag(zket(7,6)));
GammaMult = blkdiag(zeros(3), diag(zket(7,7)));
GammaVAC = blkdiag(diag([0,0,1]), zeros(7));

POVMsB = {GammaH, GammaV, GammaD, GammaA, GammaR, GammaL, GammaMult, GammaVAC};

mustBePOVM(POVMsB);

newParams.POVMA = POVMsA;
newParams.POVMB = POVMsB;


% each POVM element is assigned an announcement made by their respective
% party

newParams.announcementsA = ["Z", "Z", "X", "X"];
newParams.announcementsB = ["Z","Z","X","X","Y","Y", "Mult", "VAC"];

% newParams.keyMap = [ KeyMapElement("Z", "Z", [1, 2]) ]; %% Z basis used for key generation 
newParams.keyMap = [KeyMapElement("Z","Z",[1,2,1,2]), KeyMapElement("X","X",[1,2,1,2])];


%Define joint observables
observablesJoint = cell(numel(POVMsA), numel(POVMsB)); 
for a = 1:numel(POVMsA)
    for b = 1:numel(POVMsB)
        observablesJoint{a,b} = kron(POVMsA{a}, POVMsB{b}); 
    end
end     

newParams.observablesJoint = observablesJoint; 

%% Kraus Ops (for G map)

%K_X shrunk via iso = |0X2|_RA + |1X3|_RA. 
sqrtHV = diag([sqrt(pzB),sqrt(pzB),0,1,1,0,0,0,0,0]);
krausOpZ = sqrt(pGenCondSingleandA(1))*zket(2,1)*zket(4,1)' + sqrt(pGenCondSingleandA(2))*zket(2,2)*zket(4,2)'; 
krausOpZ = kron(kron(krausOpZ, sqrtHV),zket(2,1)); %only comp basis used for key gen

sqrtDA = diag([sqrt(pxB),sqrt(pxB),0,0,0,1,1,0,0,0]);
krausOpX = sqrt(pGenCondSingleandA(3))*zket(2,1)*zket(4,3)' + sqrt(pGenCondSingleandA(4))*zket(2,2)*zket(4,4)';
krausOpX = kron(kron(krausOpX, sqrtDA),zket(2,2));

krausOps = {krausOpZ,krausOpX}; 

% krausOps = {krausOpZ}; 
newParams.krausOps = krausOps; 

krausSum = 0;
for index = 1:numel(krausOps)
    krausSum = krausOps{index}'*krausOps{index};
end
debugInfo.storeInfo("krausSum", krausSum);

%% Z map 
proj0 = kron(diag([1,0]), eye(dimB*2) ); %ident on B' system (B - flag states)
proj1 = kron(diag([0,1]), eye(dimB*2) );
keyProj = {proj0,proj1};

%Dimension of Alice's key register
newParams.dimR = 2;

%% set key map, kraus ops, and block diagonal structure in new parameters
newParams.krausOps = krausOps;
newParams.keyProj = keyProj;
newParams.blockDimsA = 4;
newParams.blockDimsB = [3, ones([1, dimB-3])]; 

%%
newParams.n = floor(pGen*params.N); %gen
newParams.m = floor((1-pGen)*params.N); %test

newParams.coarseGrain = eye([dimA*dimB,dimA*dimB]);
newParams.Sigma = size(newParams.coarseGrain, 1);
newParams.Lambda = size(newParams.coarseGrain, 2); 

% epsilons = optimalEpsVals(params.N, n_sift, params.t, dimA, params.epsSound, params.epsATfrac);
epsilons = [params.epsSound/4,params.epsSound/4,params.epsSound/4,params.epsSound/4]; 
newParams.epsilons = epsilons; 
end

function mustBeProbDistCell(input)
mustBeProbDist([input{:}])
end

function isUnitVector(vec)
if ~ismembertol(vec'*vec, 1)
    throw(MException('FourSixWCP_Description:NotUnitVector', ...
        'Computed vector is not a unit vector.'));
end
end

function mustBePOVM(S)
    arguments
        S (1, :) cell
    end
    ident = zeros(size(S{1})); 
    for i = 1:numel(S)
        ident = ident + S{i}; 
    end
    if ~all(ismembertol(real(ident), eye(size(ident))), 'all') || ~all(ismembertol(imag(ident), zeros(size(ident)),"Datascale",1),'all')
        throw(MException('mustBePOVM:NotPOVM', ...
            'The given cell aray does not sum to identity.'));
    end
end
