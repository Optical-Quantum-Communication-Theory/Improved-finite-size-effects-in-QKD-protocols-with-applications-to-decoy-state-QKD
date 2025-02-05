function qkdInput = FourSixWCPDifInt_Preset()
%Four Six Protocol with WCP, Decoys
qkdInput = QKDSolverInput();

%% Parameters

%Add misalignment
qkdInput.addFixedParameter("misalignmentAngle", 0);
%qkdInput.addScanParameter("misalignmentAngle", {0, pi/16, pi/8, 3*pi/16, pi/4})
%qkdInput.addFixedParameter("misalignmentAngle", pi/18);

%Add birefringence
% qkdInput.addFixedParameter("birefringenceAngle", 0);
% qkdInput.addScanParameter("birefringenceAngle", num2cell(linspace(0, pi/4, 11)));

%Add birefringence with frequency
qkdInput.addFixedParameter("birefringenceFreq", 0);
% qkdInput.addScanParameter("birefringenceFreq", num2cell(linspace(0,0.04,10)));

%Add dark counts
% qkdInput.addFixedParameter("darkCountRate", 0);

%Add depolarization
% qkdInput.addFixedParameter("depolarization", 0);

%Add detector efficiencies as vector
qkdInput.addFixedParameter("detectorEfficiency", ones(6,1));

%Add loss
lossdB = linspace(0, 30, 3);
lossdB = 20;
lossEta = 10.^(-lossdB/10);
qkdInput.addScanParameter("eta", num2cell(lossEta));
% qkdInput.addScanParameter('eta', num2cell(1e-2));

%Bob's probabilities of choosing each basis
qkdInput.addFixedParameter("probsB", [1/3, 1/3, 1/3]);

%Define basis choice probability of Alice (X basis)
qkdInput.addFixedParameter("pTest", 0.5);
% ptestOpt.lowerBound = 0.001;
% ptestOpt.upperBound = 0.5;
% ptestOpt.initVal = 0.0288;
% qkdInput.addOptimizeParameter("pTest", ptestOpt); % proportion of signals used for testing

%Error correction efficiency f >= 1, f=1 is at Shanon limit 
qkdInput.addFixedParameter("f",1.16);

%% Finite Correction parameters 
qkdInput.addFixedParameter("epsSound", 1e-12);
qkdInput.addFixedParameter("epsATfrac", 0.95); 

%qkdInput.addFixedParameter("N", 1e15);
%qkdInput.addScanParameter("N", num2ell(logspace(8, 16, 21)));

% qkdInput.addFixedParameter("t", 1e-9);
% qkdInput.addFixedParameter("tsift", 1e-9);
% qkdInput.addFixedParameter("tau", 1e-9*ones(4,11,3));

qkdInput.addFixedParameter("t", 0);
qkdInput.addFixedParameter("tsift", 0);
qkdInput.addFixedParameter("tau", zeros(4,8,3));

%% Decoy Parameters
musig = 1;
qkdInput.addFixedParameter("GROUP_decoysSignal_1", musig); %signal intensity H
qkdInput.addFixedParameter("GROUP_decoysSignal_2", musig); %signal intensity V
qkdInput.addFixedParameter("GROUP_decoysSignal_3", musig); %signal intensity D
qkdInput.addFixedParameter("GROUP_decoysSignal_4", musig); %signal intensity A

qkdInput.addFixedParameter("GROUP_decoys1_1", 0.1);%decoy intensity 1 H
qkdInput.addFixedParameter("GROUP_decoys1_2", 0.1); %V
qkdInput.addFixedParameter("GROUP_decoys1_3", 0.1); %D
qkdInput.addFixedParameter("GROUP_decoys1_4", 0.1); %A

qkdInput.addFixedParameter("GROUP_decoys2_1", 0);% decoy intensity 2 H
qkdInput.addFixedParameter("GROUP_decoys2_2", 0); %V
qkdInput.addFixedParameter("GROUP_decoys2_3", 0); %D
qkdInput.addFixedParameter("GROUP_decoys2_4", 0); %A

%Decoy probabilities in test rounds
qkdInput.addFixedParameter("GROUP_decoyProbsTest_1", 1/3); %p(mu1|test)
qkdInput.addFixedParameter("GROUP_decoyProbsTest_2", 1/3);
qkdInput.addFixedParameter("GROUP_decoyProbsTest_3", 1/3);

%Decoy probabilities in generation rounds
qkdInput.addFixedParameter("GROUP_decoyProbsGen_1", 1); %p(mu1|test)
qkdInput.addFixedParameter("GROUP_decoyProbsGen_2", 0);
qkdInput.addFixedParameter("GROUP_decoyProbsGen_3", 0);

%Select decoy methods
qkdInput.addFixedParameter("decoyMethods","LP");
qkdInput.addFixedParameter("FlagSquash",true);


%% Modules
% channel model
defaultChannel = @FourSixWCPDecoyDifIntChannelFunc;
channelModule = QKDChannelModule(defaultChannel); %Channel for decoys

qkdInput.setChannelModule(channelModule);

% description 
descriptionModule = QKDDescriptionModule(@FourSixWCPDifInt_Description);
qkdInput.setDescriptionModule(descriptionModule);


% key rate function
keyModuleOptions = struct(); 
keyModuleOptions.decoyTolerance = 1e-12;
keyModuleOptions.ChoiTolerance = 1e-10;
keyModuleOptions.decoySolver = "mosek";
keyModuleOptions.decoyForceSep = true;
keyModuleOptions.decoyPrecision = "high";
keyModule = QKDKeyRateModule(@FourSixWCPDifIntKeyRateFunc, keyModuleOptions);
qkdInput.setKeyRateModule(keyModule);

% optimization
optimizerMod = QKDOptimizerModule(@coordinateDescentFunc,struct("verboseLevel",0),struct("verboseLevel",0));
qkdInput.setOptimizerModule(optimizerMod);

% math solver options
mathSolverOptions = struct();
mathSolverOptions.initMethod = 1;
mathSolverOptions.maxIter = 20;
mathSolverOptions.maxGap = 1e-6;
mathSolverOptions.blockDiagonal = true;
mathSolverMod = QKDMathSolverModule(@FW2StepSolver,mathSolverOptions);
qkdInput.setMathSolverModule(mathSolverMod);

%% global options
%ErrorHandling.CatchWarn
qkdInput.setGlobalOptions(struct("errorHandling",ErrorHandling.DontCatch,"verboseLevel",1,"cvxSolver","mosek", "cvxPrecision", "high"));
% qkdInput.setGlobalOptions(struct("errorHandling",3,"verboseLevel",1,"cvxSolver", "sedumi", "cvxPrecision", "high"));