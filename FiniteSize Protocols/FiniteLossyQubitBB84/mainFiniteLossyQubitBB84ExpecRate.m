%pick the preset file to use In this case we start with a basic qubit BB84
%protocol with no loss. Feel free to open up the module and look at what it
%sets!
qkdInput = FiniteLossyQubitBB84ExpecRatePreset();

%iterate over total number of signals sent
% N_list = [1e6,1e7,1e8,1e9,1e10,1e11,1e12];

%Only one total number of signals sent used
N_list = [1e9];

for index = 1:numel(N_list)
    %Add total signals sent from list above
    qkdInput.addFixedParameter("N",N_list(index))

    % run the QKDSolver with this input
    results = MainIteration(qkdInput);    

    % save the results and preset to a file
    filestr = sprintf("Entrywise_output_%.2e",N_list(index)) + "ExpecRate.mat";
    save(filestr,"results","qkdInput");
end

% plot the result
QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log")
