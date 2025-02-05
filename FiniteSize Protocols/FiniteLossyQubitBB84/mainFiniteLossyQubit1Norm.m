%pick the preset file to use In this case we start with a basic qubit BB84
%protocol with no loss. Feel free to open up the module and look at what it
%sets!
qkdInput = FiniteLossyQubitBB841NormPreset();

%iterate over number of signals
N_list = [1e8,1e10];

%example for single number of signals
% N_list = [1e12];

for index = 1:numel(N_list)
    %Add total signals sent from list above
    qkdInput.addFixedParameter("N",N_list(index))

    %Optimize testing fraction based on total signals sent
    qkdInput.addOptimizeParameter("m" ,struct("lowerBound",1,"initVal",0.5*N_list(index),"upperBound",N_list(index)));

    %Fixed testing ratio
    % qkdInput.addFixedParameter("m" ,0.5*N_list(index));


    % run the QKDSolver with this input
    results = MainIteration(qkdInput);
    filestr = sprintf("1Norm_output_%.2e",N_list(index)) + "_t=0.mat";

    %save the results and preset to a file.;
    save(filestr,"results","qkdInput");
end

% plot the result
QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log")
