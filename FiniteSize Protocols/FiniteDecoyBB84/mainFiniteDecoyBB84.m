%pick the preset file to use In this case we start with a basic qubit BB84
%protocol with no loss. Feel free to open up the module and look at what it
%sets!
qkdInput = FiniteDecoyBB84Preset();

%iterate over total number of signals sent
N_list = [1e6,1e7,1e8,1e9,1e10,1e11];

%example for single number of signals
% N_list = [1e12];

for index = 1:numel(N_list)
    %Add total signals sent from list above
    qkdInput.addFixedParameter("N",N_list(index))

    % run the QKDSolver with this input
    results = MainIteration(qkdInput);
    filestr = sprintf("output3Decoy_%.2e",N_list(index)) + "_t=0.mat";

    % save the results and preset to a file.;
    save(filestr,"results","qkdInput");
end


% plot the result
QKDPlot.simple1DPlot(qkdInput ,results, xScaleStyle="dB",yScaleStyle="log")
