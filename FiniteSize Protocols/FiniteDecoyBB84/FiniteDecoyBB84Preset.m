function qkdInput = FiniteDecoyBB84Preset()
% FiniteDecoyBB84Preset A preset for BB84 with WCP states with decoy
% analysis. We use the simple squashed on Bob's side, and key generation is
% only from single photon pulses
qkdInput = QKDSolverInput();

%% Parameters
qkdInput.addScanParameter("transmittance", num2cell(10.^(-linspace(0,4,41))));
% qkdInput.addFixedParameter("transmittance",1);

% qkdInput.addFixedParameter("ptest", 0.0288);
ptestOpt.lowerBound = 0.001;
ptestOpt.upperBound = 0.95;
ptestOpt.initVal = 0.0288;
qkdInput.addOptimizeParameter("ptest", ptestOpt); % proportion of signals used for testing


% decoys should be in one group, which can be created with these lines:
% qkdInput.addFixedParameter("GROUP_decoys_1", 0.9); %signal intensity %0.9 for EAT
qkdInput.addOptimizeParameter("GROUP_decoys_1", struct("lowerBound",0,"initVal",0.85,"upperBound",1)); %signal intensity

% qkdInput.addFixedParameter("GROUP_decoys_1", -log(0.9)); % log of signal intensity
% qkdInput.addOptimizeParameter("GROUP_decoys_1", struct("lowerBound",0,"initVal",-log(0.9),"upperBound",10)); %log of signal intensity

qkdInput.addFixedParameter("GROUP_decoys_2", 0.02); % decoy intensity 1 %0.02 for EAT
% qkdInput.addOptimizeParameter("GROUP_decoys_2", struct("lowerBound",0,"initVal",0.1,"upperBound",1)); % decoy intensity 1

qkdInput.addFixedParameter("GROUP_decoys_3", 0.001); % decoy intensity 2 (something slightly above 0 is usually optimal.) %0.001 for EAT

% qkdInput.addFixedParameter("GROUP_decoyProbs_1",1/2);
% qkdInput.addFixedParameter("GROUP_decoyProbs_2",1/2);

qkdInput.addFixedParameter("GROUP_decoyProbs_1",1/3);
qkdInput.addFixedParameter("GROUP_decoyProbs_2",1/3);
qkdInput.addFixedParameter("GROUP_decoyProbs_3",1/3);
%these are decoyProb conditioned on test. 


qkdInput.addFixedParameter("pz",1/2);
qkdInput.addFixedParameter("f",1.16);
qkdInput.addFixedParameter("misalignmentAngle", 0);


qkdInput.addFixedParameter("darkCountRate", 0);


%%Finite-size parameters
% qkdInput.addFixedParameter("N",1e7);
% qkdInput.addScanParameter("N",num2cell([1e7,1e8,1e9,1e10,1e11,1e12,1e13]));

qkdInput.addFixedParameter("alphabet", 2); % encoding alphabet size; for qubits, this is 2

% Store all epsilons in one struct
epsilon.PE = (13/16)*1e-8; %failure probability for parameter estimation 
epsilon.bar = (1/16)*1e-8; %uncertainty from "smoothing" of min-entropy among states similar to feasible state
epsilon.EC = (1/16)*1e-8; %failure probability for error-correction
epsilon.PA = (1/16)*1e-8; %failure probability for privacy amplification
epsilonTotal = 1e-8;

qkdInput.addFixedParameter("epsilon", epsilon);
qkdInput.addFixedParameter("epsilonTotal", epsilonTotal);


qkdInput.addFixedParameter("tExp", -Inf);
qkdInput.addFixedParameter("tSiftExp",-Inf);

qkdInput.addFixedParameter("physDimAB", 2*3); % physical dimensions of Alice's and Bob's outgoing/incoming signals, used for post-selection

%Select decoy methods and define squasher
qkdInput.addFixedParameter("decoyMethods","LP");
qkdInput.addFixedParameter("FlagSquash",false);

% description is the same as the lossy qubit description since we squash
% Bob's detector data down to a lossy qubit equivalent
descriptionModule = QKDDescriptionModule(@BasicWCPBB84LossyDescriptionFunc);
qkdInput.setDescriptionModule(descriptionModule);

% channel model is very different from normal qubit
channelModule = QKDChannelModule(@BasicBB84WCPDecoyChannelFunc);
qkdInput.setChannelModule(channelModule);

% Key rate module performs squashing and decoy analysis
keyRateOptions = struct();
keyRateOptions.decoyTolerance = 1e-14;
keyRateOptions.decoySolver = "Mosek";
keyRateOptions.decoyForceSep = true;
keyMod = QKDKeyRateModule(@FiniteDecoyBB84KeyRateFunc, keyRateOptions);
qkdInput.setKeyRateModule(keyMod);

optimizerMod = QKDOptimizerModule(@coordinateDescentFunc,struct("verboseLevel",0),struct("verboseLevel",0));
qkdInput.setOptimizerModule(optimizerMod);

% math solver options
mathSolverOptions = struct();
mathSolverOptions.initMethod = 1;
mathSolverOptions.maxIter = 20;
mathSolverOptions.maxGap = 1e-6;
mathSolverOptions.blockDiagonal = false;
% mathSolverOptions.infiniteDecoyIntensities = true; % Good for testing
% with an infinite number of decoy intensities (Upper and lower decoy
% bounds converge.)
mathSolverMod = QKDMathSolverModule(@FW2StepSolver,mathSolverOptions,mathSolverOptions);
qkdInput.setMathSolverModule(mathSolverMod);

% global options
qkdInput.setGlobalOptions(struct("errorHandling",3,"verboseLevel",1,"cvxSolver","mosek"));
