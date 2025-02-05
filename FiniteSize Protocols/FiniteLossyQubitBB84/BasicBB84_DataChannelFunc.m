function [newParams, modParser]= BasicBB84_DataChannelFunc(params,options, debugInfo)
% BB84LossyChannelFunc A channel function for qubibt BB84 with loss added
% as a perpendicular dimension. Alice's dimensions have been reduced from 4
% to 2 via the Schmidt decomposition, similar to BB84ChannelFunc
%
% Input parameters:
% * dimA: dimension of Alice's system.
% * dimB: dimension of Bob's system.
% * observablesJoint: The joint observables from Alice and Bob's
%   measurments which they perform on the (idealy) max entangled state. The
%   observables must be hermitian and each must be the size dimA*dimB by
%   dimA*dimB. The observables assume the spaces are ordered A \otimes B.
%   They also should be positive semi-definite and should sum to identity,
%   but this is hard to check because of machine precision issues.
% * eta (1): the transmissivity of the quantum channel. What fraction of
%   the signal reaches Bob's detectors. Equivalent to 1 - loss. Must be
%   between 0 and 1 inclusive.
% * depolarization (0): The amount of depolarization applied to the signal
%   Alice sends to Bob. At maximum depolarization (depolariztion =1) a pure
%   qubit state is converted to a maximally mixed state. Depolarization
%   should be between 0 and 1.
% * misalignmentAngle (0): Angle Alice and Bob's bases are misaligned by
%   around the Y-axis. For example, Bob's detectors could be slightly
%   rotated away from the incoming signals. Although calculations are done
%   on the Bloch sphere, angles should not be given in that form (period
%   4pi). This angle is measured as the physical rotation of the device
%   (period 2pi).
% Output parameters:
% * expectationsJoint: The joint epxectations for Alice and Bob's
%   measurement of the signals. Simply formed by taking the
%   observablesJoint and applying them to a simulated rhoAB.
% Options:
% * none
% DebugInfo:
% * rhoAB: Alice and Bob's shared density matrix after the channel has
%   acted on it. Usefull for checking the channel has been applied
%   correctly.
%
% See also QKDChannelModule, BB84ChannelFunc, makeGlobalOptionsParser
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
modParser.addRequiredParam("observablesJoint",@(x) allCells(x,@ishermitian))
modParser.addOptionalParam("eta", 1, @(x) mustBeInRange(x, 0, 1));
modParser.addRequiredParam("dimA",@(x) x==2);
modParser.addRequiredParam("dimB", @(x) x ==3); %sanity check
modParser.addOptionalParam("depolarization",0,@(x) mustBeInRange(x,0,1));
modParser.addOptionalParam("misalignmentAngle",0,@(x) mustBeReal(x));

modParser.addAdditionalConstraint(@observablesAndDimensionsMustBeTheSame,["observablesJoint","dimA","dimB"]);

modParser.parse(params);

params = modParser.Results;

%% simple setup

newParams = struct();

observablesJoint = params.observablesJoint;
dimA = params.dimA;


%% load data

expectationsLower = load("filenameLower");
expectationsUpper = load("filenameUpper");

%Store data
newParams.expectationsLowerBound = expectationsLower;
newParams.expectationsUpperBound = expectationsUpper;

end

function observablesAndDimensionsMustBeTheSame(observables,dimA,dimB)
if ~allCells(observables,@(x) size(x,1) == dimA*dimB)
    throwAsCaller(MException("BasicBB84_LossyChannelFunc:ObservablesAndDimensionsMustBeTheSame",...
        "The Observables must have the same dimensions as Alice and Bob multiplied together."));
end
end