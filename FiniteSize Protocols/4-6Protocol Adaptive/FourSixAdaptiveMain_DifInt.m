% Decoy 4-6 protocol prepare and measure with different intensities for
% signal states on the source side.

qkdInput = FourSixWCPDifInt_AdaptivePreset();

%iterate over number of signals
N_list = [1e9,1e10,1e11,1e12];

%example for single number of signals
% N_list = [1e12];

%Loss
%total array of loss values to iterate over
lossdB = linspace(0, 40, 21);
lossEta = 10.^(-lossdB/10);

% We only use array lossEta up until e.g. 5th entry because key rates are 0 for higher loss
% lossList stores maximum entry of lossEta used for key rate calculation. This is different for iid and postselection (PS).

%When switching between iid and postselection, lifttype in preset needs to be changed!

%Postselection
% lossList = [5,10,13,17]; %PS

%IID
lossList = [15,18,21,21]; %iid

%Example for single value
% lossList = [21];

for index = 1:numel(N_list)
    disp(N_list(index))
    %Add total signals sent from list above
    qkdInput.addFixedParameter("N",N_list(index));
    
    %Add loss until element from list above
    qkdInput.addScanParameter("eta", num2cell(lossEta(1:lossList(index))));

    % run the QKDSolver with this input
    results = MainIteration(qkdInput);

    % save the results and preset to a file
    filestr = sprintf("AdaptiveDecoy46Protocol_iid_%.2e",N_list(index)) + ".mat";
    save(filestr,"results","qkdInput");
end


QKDPlot.plotParameters({results},"eta","_keyRate_","addLegend",false,"xScaleStyle","dB","yScaleStyle","log");
