function qkdInput = FiniteLossyQubitBB84ExpecRatePreset()
%  A preset for a Lossy Qubit BB84 protocol with finite size 
%  
qkdInput = QKDSolverInput();

%% Parameters
% qkdInput.addScanParameter("transmittance", num2cell(10.^(-linspace(0,4,41))));
qkdInput.addFixedParameter("transmittance", 10^(-1))
qkdInput.addFixedParameter("misalignmentAngle",0); %check physical vs bloch sphere rotations. 
qkdInput.addFixedParameter("depolarization", 0);

qkdInput.addFixedParameter("pz",1/2);
qkdInput.addFixedParameter("f",1.16);

%choose variant.


%% Finite size parameters
% qkdInput.addFixedParameter("N", 1e9); % number of signals sent
% qkdInput.addScanParameter("N",num2cell([1e5,1e6,1e7,1e8,1e9,1e10,1e11]));

ptestOpt.lowerBound = 0.001;
ptestOpt.upperBound = 0.95;
ptestOpt.initVal = 0.0288;
qkdInput.addOptimizeParameter("ptest", ptestOpt); % proportion of signals used for testing
% qkdInput.addFixedParameter("ptest", 0.0288);

qkdInput.addFixedParameter("alphabet", 2); % encoding alphabet size; for qubits, this is 2

% Store all epsilons in one struct
epsilon.PE = (13/16)*1e-8; %failure probability for parameter estimation 
epsilon.bar = (1/16)*1e-8; %uncertainty from "smoothing" of min-entropy among states similar to feasible state
epsilon.EC = (1/16)*1e-8; %failure probability for error-correction
epsilon.PA = (1/16)*1e-8; %failure probability for privacy amplification
epsilonTotal = 1e-8;

qkdInput.addFixedParameter("epsilon", epsilon);
qkdInput.addFixedParameter("epsilonTotal", epsilonTotal);


% qkdInput.addFixedParameter("tExp", -Inf)
qkdInput.addScanParameter("tExp",num2cell(linspace(-9,-4,100)));
% qkdInput.addFixedParameter("tSiftExp", -Inf);

qkdInput.addFixedParameter("physDimAB", 2*3); % physical dimensions of Alice's and Bob's outgoing/incoming signals, used for post-selection

% end finite size parameters


% description and channel models are independent of finite size analysis,
% so we use the same description and channel as the lossy qubit
descriptionModule = QKDDescriptionModule(@BasicBB84LossyDescriptionFunc);
qkdInput.setDescriptionModule(descriptionModule);

channelModule = QKDChannelModule(@BasicBB84LossyChannelFunc);
qkdInput.setChannelModule(channelModule);

% keyrate function and solver choice change with finite size because we
% need to change certain observables to uncertain observables
keyMod = QKDKeyRateModule(@FiniteLossyQubitBB84ExpectedKeyRateFunc);
qkdInput.setKeyRateModule(keyMod);

optimizerMod = QKDOptimizerModule(@coordinateDescentFunc,struct("verboseLevel",0),struct("verboseLevel",0));
qkdInput.setOptimizerModule(optimizerMod);

% options for math solver
mathSolverOptions = struct();
mathSolverOptions.initMethod = 1;
mathSolverOptions.maxIter = 10;
mathSolverOptions.maxGap = 1e-6;
mathSolverOptions.cvxSolver = "mosek";
mathSolverMod = QKDMathSolverModule(@FW2StepSolver,mathSolverOptions);
qkdInput.setMathSolverModule(mathSolverMod);

% global options
qkdInput.setGlobalOptions(struct("errorHandling",3,"verboseLevel",1));