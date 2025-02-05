% Decoy 4-6 protocol prepare and measure. We use Z-basis for key generation
% and Z-basis for testing only
qkdInput = FourSixWCPDifInt_Preset();

% etaFull = 10.^(-linspace(0,4,21));
% etaFull(1) = 0.95;

% N_list = [1e7,1e8,1e9,1e10,1e11,1e12];
% etaIndex = [7,12,21,21,21,21];
% N_list = [1e9,1e10,1e11,1e12];
% etaIndex = [21,21,21,21];

N_list = [1e13];
% etaIndex = [8];

for index = 1:numel(N_list)
    %Add total signals sent from list above
    
    qkdInput.addFixedParameter("N",N_list(index))
    
    % qkdInput.addFixedParameter("eat",etaFull(1))
    % qkdInput.addScanParameter("eta", num2cell(etaFull(1:etaIndex(index))));

    % run the QKDSolver with this input
    results = MainIteration(qkdInput);
    filestr = sprintf("output2Decoy46prot_%.2e",N_list(index)) + ".mat";

    %save the results and preset to a file.;
    % save(filestr,"results","qkdInput");
end

%% plot the result
QKDPlot.simple1DPlot(qkdInput ,results, xScaleStyle="dB",yScaleStyle="log")