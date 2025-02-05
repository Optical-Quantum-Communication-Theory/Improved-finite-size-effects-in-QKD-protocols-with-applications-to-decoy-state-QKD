function qkdInput = FourSixWCPDifInt_AdaptivePreset()
%Four Six Protocol with WCP, Decoys
% We compute Adaptive Keyrate for the expected channel behaviour. 
qkdInput = QKDSolverInput();

%% Parameters

%Add misalignment
qkdInput.addFixedParameter("misalignmentAngle", 0.01);
%qkdInput.addScanParameter("misalignmentAngle", {0, pi/16, pi/8, 3*pi/16, pi/4})
%qkdInput.addFixedParameter("misalignmentAngle", pi/18);



%Add dark counts
% qkdInput.addFixedParameter("darkCountRate", 0);

%Add depolarization
% qkdInput.addFixedParameter("depolarization", 0);

%Add detector efficiencies as vector
qkdInput.addFixedParameter("detectorEfficiency", ones(6,1));

%Add loss
% lossdB = linspace(0, 40, 21);
% % lossdB = 0;
% lossEta = 10.^(-lossdB/10);
% qkdInput.addScanParameter("eta", num2cell(lossEta));
%qkdInput.addFixedParameter('eta', 1e-1);


%Bob's probabilities of choosing each basis
qkdInput.addFixedParameter("probsB", [1/3, 1/3, 1/3]);

%Define testing probability of Alice 
%fixed testing
% qkdInput.addFixedParameter("pTest", 0.05);

%optimze testing fraction
ptestOpt.lowerBound = 0.001;
ptestOpt.upperBound = 0.95;
ptestOpt.initVal = 0.0288;
qkdInput.addOptimizeParameter("pTest", ptestOpt);

%Error correction efficiency f >= 1, f=1 is at Shanon limit 
qkdInput.addFixedParameter("f",1.16);

%% Finite Correction parameters


numDecoy = 2;  %make sure this is total number of deocy intensities used
decoyPhotonCutoff = 2; %photon number cutoff for decoy analysis.
dA = 4; % Alice has 4 signal states.

%Bob has 1 block of dimension 2, 1 block of dimension 1 (vacuum),
% and total 7 classical flags (6 for single, 1 for multiclicks);
% We assume that no-click flag is in the vacuum component.


liftDimension = numDecoy^2*(decoyPhotonCutoff+2)*(dA^2*5  + dA^2*9); 
% %according to current formula
qkdInput.addFixedParameter("liftDimension",liftDimension); %


log2targetepsSec = -12 / log10(2); %10^{-12}
log2targetepsCor = -12 / log10(2);  % heuristic choice. 
qkdInput.addFixedParameter("log2targetepsSec",log2targetepsSec);
qkdInput.addFixedParameter("log2targetepsCor",log2targetepsCor);
qkdInput.addFixedParameter("liftType","IID"); %IID
% qkdInput.addFixedParameter("liftType","PS"); %Postselection


epsBiasFactor = 1/2;
qkdInput.addFixedParameter("epsBiasFactor",epsBiasFactor); %adjusts bias between (epsPA), and (epsAT)

%qkdInput.addFixedParameter("alpha", 1+1/sqrt(N)); %Must be chosen BEFORE
%protocol is run. For now, we can choose the best value of alpha (which
%depends on number of sifted signals observed) in the key rate function.


% qkdInput.addFixedParameter("N", 1e15);



%% Decoy Parameters
musig = 0.9;
qkdInput.addFixedParameter("GROUP_decoysSignal_1", 0.9); %signal intensity H
qkdInput.addFixedParameter("GROUP_decoysSignal_2", 0.5); %signal intensity V
qkdInput.addFixedParameter("GROUP_decoysSignal_3", 0.6); %signal intensity D
qkdInput.addFixedParameter("GROUP_decoysSignal_4", 1.1); %signal intensity A

qkdInput.addFixedParameter("GROUP_decoys1_1", 0.1);%decoy intensity 1 H
qkdInput.addFixedParameter("GROUP_decoys1_2", 0.2); %V
qkdInput.addFixedParameter("GROUP_decoys1_3", 0.01); %D
qkdInput.addFixedParameter("GROUP_decoys1_4", 0.09); %A

% qkdInput.addFixedParameter("GROUP_decoys2_1", 0.01);% decoy intensity 2 H
% qkdInput.addFixedParameter("GROUP_decoys2_2", 0.01); %V
% qkdInput.addFixedParameter("GROUP_decoys2_3", 0.01); %D
% qkdInput.addFixedParameter("GROUP_decoys2_4", 0.01); %A

%Decoy probabilities in test rounds
qkdInput.addFixedParameter("GROUP_decoyProbsTest_1", 1/2); %p(mu1|test)
qkdInput.addFixedParameter("GROUP_decoyProbsTest_2", 1/2);
% qkdInput.addFixedParameter("GROUP_decoyProbsTest_3", 1/3);

%Decoy probabilities in generation rounds
qkdInput.addFixedParameter("GROUP_decoyProbsGen_1", 1); %p(mu1|gen)
qkdInput.addFixedParameter("GROUP_decoyProbsGen_2", 0);
% qkdInput.addFixedParameter("GROUP_decoyProbsGen_3", 0);

%Select decoy methods
qkdInput.addFixedParameter("decoyMethods","SDP");
qkdInput.addFixedParameter("FlagSquash",true);


%% Modules
% channel model
defaultChannel = @FourSixWCPDecoyAdaptiveDifIntChannelFunc;
channelModule = QKDChannelModule(defaultChannel);

qkdInput.setChannelModule(channelModule);

% description 
descriptionModule = QKDDescriptionModule(@FourSixWCPDifInt_AdaptiveDescription);
qkdInput.setDescriptionModule(descriptionModule);


% key rate function
keyModuleOptions = struct(); 
keyModuleOptions.decoyTolerance = 1e-12;
keyModuleOptions.photonCutOff = decoyPhotonCutoff; 
keyModuleOptions.ChoiTolerance = 1e-10;
keyModuleOptions.decoySolver = "mosek";
keyModuleOptions.decoyForceSep = true;
keyModuleOptions.decoyPrecision = "high";
keyModule = QKDKeyRateModule(@FourSixWCPDifIntAdaptiveKeyRateFunc, keyModuleOptions);
qkdInput.setKeyRateModule(keyModule);

% optimization
optimizerMod = QKDOptimizerModule(@coordinateDescentFunc,struct("verboseLevel",0),struct("verboseLevel",0));
qkdInput.setOptimizerModule(optimizerMod);

% math solver options
mathSolverOptions = struct();
mathSolverOptions.initMethod = 1;
mathSolverOptions.maxIter = 40;
mathSolverOptions.maxGap = 1e-6;
mathSolverOptions.blockDiagonal = true;
mathSolverMod = QKDMathSolverModule(@FW2StepSolver,mathSolverOptions);
qkdInput.setMathSolverModule(mathSolverMod);

%% global options
%ErrorHandling.CatchWarn
qkdInput.setGlobalOptions(struct("errorHandling",ErrorHandling.DontCatch,"verboseLevel",1,"cvxSolver","mosek", "cvxPrecision", "default"));
% qkdInput.setGlobalOptions(struct("errorHandling",3,"verboseLevel",1,"cvxSolver", "sedumi", "cvxPrecision", "high"));
