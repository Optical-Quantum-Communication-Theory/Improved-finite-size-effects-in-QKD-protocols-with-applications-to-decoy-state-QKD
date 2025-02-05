function [newParams, modParser]= FourSixWCPDecoyOverpassChannel_DifInt_Func(params,options,debugInfo)
% FourSixWCPDecoyChannelFunc A channel function for the 4-6 protocol using WCP
% states, supporting decoy intensities. Given a collection of decoy
% intensities, this channel produces a group of 4x64 tables of
% expectations, one for each decoy intensity, which are the conditional
% probability for each of Bob's 64 detector patterns given Alice's signal
% sent (and the decoy intenisty).
%
% Input parameters:
% * probsB (1/3*ones(3,1)): probabilities of Bob choosing each basis in
%   order Z,X,Y
% * decoys: a cell of the intensities used in decoy analysis. These are the
%   mean photon numbers that Alice can choose from when performing the
%   decoy protocol. The first element in the cell array is treated as the
%   intensity used for key generation.
% * lossfile: model of the loss of a overpass
% * timestep (0): time step in of overpass
% * detectorEfficiency (1): the efficiency of Bob's detectors. Must be
%   between 0 and 1 inclusive.
% * misalignmentAngle (0):  Physical angle of misalignment between Alice
%   and Bob's measurements around Y axix. This angle is measured as the
%   physical rotation of the device (period 2pi). Although calculations are
%   done on the Bloch sphere, angles should not be given in that form
%   (period 4pi).
% * birefringenceFreq (0):  Rotation frequency of birefringence between Alice
%   and Bob's measurements around Z axix. The resulting angle is measured as the
%   physical rotation of the device (period 2pi). Although calculations are
%   done on the Bloch sphere, angles should not be given in that form
%   (period 4pi).
% * detectorEfficiency (ones(6,1)): detector efficiency of each of Bob's
%   detectors
% * darkCountRate (0): The probability that a detector that recieves no
%   photons will still randomly click anyway. Must be between 0 and 1.
% Output parameters:
% * expectationsConditional: The conditional expectations (as a 3D array)
%   from Alice and Bob's measurements. This should be organized as a 4 x 64
%   x n array, where 4 = number of signals Alice sent, 64 = Bob's detector
%   click patterns, and n = the number of intensities used in the decoy
%   protocol. The Table is conditioned on the signal Alice sent, and the
%   intensity she chose. Therefore, each row should sum to 1.
% Options:
% * None.
% DebugInfo:
% * transMat: The linear operator that trasforms the mode operators from
%   what Alice sent, to what Bob recieves. This includes Bob's detector
%   setup, except for the non-linear dark counts. See Coherent for more
%   details.
% * probDetectorClickCon: The probability of each individual detector
%   clicking given the signal choice Alice sent (and intensity in dim 3).
%
% 
% See also QKDChannelModule, Coherent
arguments
    params (1,1) struct
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end

%% options parser
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;

%% module parser
modParser = moduleParser(mfilename);

%Decoy intensities
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

%Bob's basis choices
modParser.addRequiredParam("probsB", @(x) mustBeProbDist(x));

%Time step in overpass
modParser.addRequiredParam("timeStart", @(x) mustBeReal(x));
modParser.addRequiredParam("timeEnd", @(x) mustBeReal(x));
modParser.addRequiredParam("timeInterpol", @(x) mustBeReal(x));

%Detector efficiencies as vector
modParser.addOptionalParam("detectorEfficiency", ones(6,1), @(x) mustBeInRange(x, 0, 1));
modParser.addAdditionalConstraint( @(x) mustBeEqualSize(x, ones(6,1)) ,"detectorEfficiency");

%Misalingment angle
modParser.addOptionalParam("misalignmentAngle",0,@mustBeReal);
modParser.addAdditionalConstraint(@isscalar,"misalignmentAngle");

%Birefringence angle
modParser.addOptionalParam("birefringenceFreq",0,@mustBeReal);
modParser.addAdditionalConstraint(@isscalar,"birefringenceFreq");

%Darkcount rate
modParser.addOptionalParam("darkCountRate", 0, @(x) mustBeInRange(x, 0, 1));
modParser.addAdditionalConstraint(@isscalar,"darkCountRate");

%Channel loss
modParser.addOptionalParam("eta", 1, @(x) mustBeInRange(x, 0, 1));
modParser.addAdditionalConstraint(@isscalar,"eta");

%Load loss file
modParser.addOptionalParam("LoadLossFile", @(x) mustBeMember(x, [true,false]));


modParser.parse(params);

params = modParser.Results;

if params.LoadLossFile == true
    disp("!!! Loading Loss File !!!")

    % %load loss file from Jennewein group
    % matData = load("Jennewein Data\LinkData.csv");
    % lossData = matData(:,1);
    % timeintervals = matData(:,3);
    % 
    % %construct model from loss data using spline interpolation and evaluate at
    % %time = timestep
    % params.eta = 10.^(-1/10*interp1(timeintervals,lossData,params.timestep,'spline'));
    
    %load loss file from Jas
    matData = csvread("ReFQ_code\data.csv",1,0);
    lossData = matData(:,4);
    timeintervals = matData(:,1);
else
    timeintervals = params.timeStart:params.timeInterpol:params.timeEnd;
    lossData = params.eta*ones(size(timeintervals));
end


% construct the signals
    
% we will use an intensity of 1 for now as we can scale that up after the
% fact.
signals = {Coherent.pauliCoherentState(1,1,1);... %H
        Coherent.pauliCoherentState(1,1,2);... %V
        Coherent.pauliCoherentState(1,2,1);... %D
        Coherent.pauliCoherentState(1,2,2)}; %A

%total number of intensities
numtotdecoy = numel([params.decoysSignal{1},params.decoys1{1},params.decoys2{1}]);

probEachDetectorClickConFinal = zeros(numel(signals),size(params.detectorEfficiency,1),numtotdecoy);
expectationsConFinal = zeros(numel(signals),2^size(params.detectorEfficiency,1),numtotdecoy);

for timestep = params.timeStart:params.timeInterpol:params.timeEnd
    %Total number of intervals
    nIntervals = floor((params.timeEnd-params.timeStart)/params.timeInterpol)+1;

    %construct model from loss data using spline interpolation and evaluate at
    %time = timestep (works also for constant loss)
    eta = interp1(timeintervals,lossData,timestep,'spline');
    
    %Determine birefringence angle from timestep
    params.birefringenceAngle = params.birefringenceFreq * timestep;    
    
    % build the (sub) isometry transition matrix that represents the channel
    % and Bob's measurement except for the dark counts which must be handled
    % later.
    [transMat,detectorMat] = sixstateLinearOpticsSetup(eta,params.misalignmentAngle,params.birefringenceAngle,params.detectorEfficiency,params.probsB);
    debugInfo.storeInfo("transMat",transMat);
    debugInfo.storeInfo("detectorMat",detectorMat);
    
    newParams.detectorMat = detectorMat;
    
    
    %Calculate the conditional probabilities for the click patterns
    reshapedIntensities = {params.decoysSignal,params.decoys1,params.decoys2};
    
    probEachDetectorClickCon = zeros(numel(signals),size(params.detectorEfficiency,1),numtotdecoy);
    expectationsCon = zeros(numel(signals),2^size(params.detectorEfficiency,1),numtotdecoy);
    
    for indexDecoySet = 1:numtotdecoy
        %scale the signal states for the intensity
        sqrtInt = cellfun(@sqrt, reshapedIntensities{indexDecoySet}, 'UniformOutput', false);
        signalsDecoy = cellfun(@times, signals, sqrtInt.', 'UniformOutput', false); 
        [expectationsCon(:,:,indexDecoySet), probEachDetectorClickCon(:,:,indexDecoySet),patternOrder] = simulateChannel(signalsDecoy,transMat,params.darkCountRate);
    end
    
    %Add to final expectations
    expectationsConFinal = 1/nIntervals*expectationsCon + expectationsConFinal;
    probEachDetectorClickConFinal = 1/nIntervals*probEachDetectorClickCon + probEachDetectorClickConFinal;

end

debugInfo.storeInfo("probEachDetectorClickCon",probEachDetectorClickConFinal);

debugInfo.storeInfo("expectationsConditional",expectationsConFinal);

newParams.expectationsConditional = expectationsConFinal;

newParams.patternOrder = patternOrder;


end

function [transMat,detectorMat] = sixstateLinearOpticsSetup(eta,misalignmentAngle,birefringenceAngle,detectorEfficiency,probsB)

% probDist = (pzB, pxB, pyB)
% detectorEfficiency = (deteffH, ..., deteffL), 6x1 vector

%% construct channel transition marix
%loss/transmittance
channelMat = Coherent.copyChannel(Coherent.transmittanceChannel(eta),2);

%misalignment rotation (rotation around Y)
channelMat = Coherent.rotateStateZXY(misalignmentAngle,[0,0,1],"angleOnBlochSphere",false)*channelMat;

%birefringence rotation (rotation around Z)
channelMat = Coherent.rotateStateZXY(birefringenceAngle,[1,0,0],"angleOnBlochSphere",false)*channelMat;

%% Build up Bob's detector transition matrix

% Bob applies a beam splitter to send signals to each detector basis setup
detectorMat = Coherent.copyChannel(Coherent.singleInputMultiBeamSpliter(probsB),...
    2,"weaveCopies",true);

% Bob applies a rotation to convert A and D back to H and V for easier
% measurement
detectorMat = blkdiag(pauliBasis(1,false)',pauliBasis(2,false)', pauliBasis(3,false)')*detectorMat;

% We have to handle dark counts after we get the click probabilities for
% each detector, so no changes here.

% Each detector has a different efficiency so we can't pull it right to the
% start.
detectorMat = Coherent.transmittanceChannel(detectorEfficiency)*detectorMat;

transMat = detectorMat*channelMat;
end

function probDetectorClickCon = applyDarkCounts(probDetectorClickCon,darkCountRate)
probDetectorClickCon = 1-(1-probDetectorClickCon)*(1-darkCountRate);
end

function [probDetectorClickPatternCon, probEachDetectorClicksCon, patternOrder]  = simulateChannel(signals,transMat,darkCountRate)
    %Construct the independent detector click probabilities for each signal
    probEachDetectorClicksCon = detectorClickProbabilities(signals,transMat);
    %simulate the effects of dark counts
    probEachDetectorClicksCon  = applyDarkCounts(probEachDetectorClicksCon,darkCountRate);
    %Construct all combinations of detector firing patterns from the
    %independent detectors.
    [probDetectorClickPatternCon, patternOrder]  = detectorClickPatterns(probEachDetectorClicksCon);
end



function probDetectorClickCon = detectorClickProbabilities(signals,transMat)

probDetectorClickCon = zeros(numel(signals),size(transMat,1));

for index = 1:numel(signals)
    bobsSignal = transMat*signals{index};
    probDetectorClickCon(index,:) = 1-Coherent.fockCoherentProb(zeros(size(transMat,1),1),bobsSignal,"combineModes",false);
end
end


function [probDetectorClickPatternCon,patternOrder] = detectorClickPatterns(probClickCon)
%Because we have coherent states, each detector acts independently. This
%function takes the independent results from each detector and computes the
%probabilities of each click pattern outcome.
numSignals = size(probClickCon,1);
numDetectors = size(probClickCon,2);

probDetectorClickPatternCon = zeros(numSignals,2^numDetectors);

sizeDetectorPatterns = 2*ones(1,numDetectors);
patternOrder = cell(1, 2^numDetectors);

clickProbSwitch = @(click, clickProb) (click==0).*(1-clickProb) + (click~=0).*clickProb;
for signalIndex = 1:numSignals
    for indexPat = 1:2^numDetectors
        patternVec = ind2subPlus(sizeDetectorPatterns,indexPat)-1;
        patternOrder{indexPat} = patternVec; 
        probDetectorClickPatternCon(signalIndex,indexPat) = prod(clickProbSwitch(patternVec,probClickCon(signalIndex,:)));
    end
end
end