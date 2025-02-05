%helper function that performs decoy state analysis
function [Y1L,Y1U] = decoySDPFiniteDifInt(expecJointTest,statesAlice,POVMsBob,decoys,probListAlice,probsDecoy,ptest,t,muLower,muUpper,photonBound,FlagStateSquash,options)
    arguments
        expecJointTest (:,:,:) double {mustBeGreaterThanOrEqual(expecJointTest,0)}
        statesAlice (1,:) cell;
        POVMsBob (1,:) cell;
        decoys (:,:) cell;        
        probListAlice (:,1) double {mustBeProbDist}
        probsDecoy (:,1) double {mustBeProbDist}
        ptest (1,1) double {mustBeInRange(ptest,0,1),mustBeSameTotalProb(ptest,expecJointTest)}        
        t (:,:,:) double {mustBeGreaterThanOrEqual(t,0),SameNumberOfDecoysAndPages(decoys,t)}
        muLower (:,:,:) double {mustBeGreaterThanOrEqual(muLower,0),SameNumberOfDecoysAndPages(decoys,muLower)}
        muUpper (:,:,:) double {mustBeGreaterThanOrEqual(muUpper,0),SameNumberOfDecoysAndPages(decoys,muUpper)}
        photonBound (:,2) double {mustBeNonnegative,mustBeLessThanOrEqual(photonBound,1)}
        FlagStateSquash (1,1) logical = false;
        options.decoyTolerance (1,1) double {mustBeGreaterThanOrEqual(options.decoyTolerance,0)}= 1e-10;
        options.ChoiTolerance (1,1) double {mustBeGreaterThanOrEqual(options.ChoiTolerance,0)}= 1e-10;
        options.decoySolver (1,1) string = "Mosek";
        options.decoyPrecision (1,1) string = "medium";
        options.photonCutOff (1,1) double {mustBeInteger,mustBePositive} = 10;
        options.verboseLevel (1,1) double {mustBeInteger,mustBeNonnegative} = 1;
    end
   
    %number of intensities used
    n_decoy=numel(decoys);

    %Dimension of Bob
    dimB_SDP = size(POVMsBob{1},1);

    %Dimension of Alice
    dimA_SDP = size(statesAlice{1},1);

    %m = #number of rows of observations
    %n = #number of colums of observations
    m = size(expecJointTest,1);
    n = size(expecJointTest,2);
    
    %Empty vectors for upper and lower bounds
    Y1L = zeros(m,n);
    Y1U = ones(m,n);
    
    %cut-off for decoy
    n_photon= options.photonCutOff;
    
    %Number of intensities
    n_signal_int = size(decoys{1},2);

    % %Possonian distribution
    probDist = @(intensity,n) exp(-intensity).*intensity.^n./factorial(n);

    %Thermal distribution
    % probDist=@(meanPh,n) meanPh.^n/(1+meanPh).^(1+n);
    
    %Tolerance in Decoy SDP
    decoy_tolerance = options.decoyTolerance;
    choi_tolerance = options.ChoiTolerance;

    %Pick n=0,1,...-photon components
    M = {1,n_photon+1};
    for i = 1:n_photon+1
        Mi = zeros(m, m*(n_photon + 1));
        indr = [1:m]';
        indc = [i:n_photon+1:m*(n_photon + 1)]';
        indx = sub2ind(size(Mi),indr,indc);
        Mi(indx) = 1;
        M{i}=Mi;
    end

    %Signal density matrices
    rhoAlice = cell(1,numel(statesAlice));
    for index = 1:numel(statesAlice)
        rhoAlice{index} = statesAlice{index}*statesAlice{index}';
    end

    %Calculate upper and lower bounds on the expectations
    upperbndsExpec = ones(m,n,numel(decoys));
    lowerbndsExpec = zeros(m,n,numel(decoys));

    for kdecoy = 1:numel(decoys)
        for row = 1:numel(probListAlice)
            upperbndsExpec(row,:,kdecoy) = min(max((expecJointTest(row,:,kdecoy) + t(row,:,kdecoy) + muUpper(row,:,kdecoy))/(probListAlice(row)*probsDecoy(kdecoy)*ptest),0),1);
            lowerbndsExpec(row,:,kdecoy) = min(max((expecJointTest(row,:,kdecoy) - t(row,:,kdecoy) - muLower(row,:,kdecoy))/(probListAlice(row)*probsDecoy(kdecoy)*ptest),0),1);
        end
    end
    
    %Projection onto n<=N_B subspace
    Proj_sub = blkdiag(eye(3),zeros(dimB_SDP-3));
    
    %solve for upper bound 1-photon
    for i=0:m-1
        for j=1:n
            try
                %Select the correct rows and columns from Y for the objective
                %function
                Objcolumn = zeros(n,1);
                Objcolumn(j) = 1;
                
                Objrow = zeros(1,m*(n_photon+1));
                Objrow(i*(n_photon+1)+2) = 1;

                %Set up SDP for upper bounds
                % Y is the matrix representing all yields
                % J0 is the Choi matrix for the 0-photon component
                % J1 is the Choi matrix for the 1-photon component
                % In both Choi matrices Bob's system is first and Alice's
                % system is second to satisfy e.g. Y_0 = Tr[J_0 (F_y otimes (rho_x)^T)] 
                cvx_begin quiet
                    %CVX solver
                    cvx_solver(convertStringsToChars(options.decoySolver));
                    variable Y(m*(n_photon+1),n) nonnegative
                    variable J0(dimB_SDP,dimB_SDP) hermitian semidefinite
                    variable J1(dimA_SDP*dimB_SDP,dimA_SDP*dimB_SDP) hermitian semidefinite
                    maximize( real(Objrow*Y*Objcolumn))
                    subject to
                        % 0 <= Y <=1
                        vec(Y) >= vec(zeros(m*(n_photon+1),n));
                        vec(Y) <= vec(ones(m*(n_photon+1),n));
                        
                        %0- and 1-photon component treated seperately with add.
                        %constraints
                        Y0 = M{1}*Y;
                        Y1 = M{2}*Y;
                        
                        %Partial trace of Choi matrix = id
                        norm(PartialTrace(J1,[1],[dimB_SDP,dimA_SDP])-eye(dimA_SDP)) <= choi_tolerance;
                        norm(trace(J0) - 1) <= choi_tolerance;
                        
                        %Usual decoy bounds rewritten as matrix times vector
                        for k = 1:n_decoy
                            %total probability for given intensity to have
                            %<= n_photon photons
                            %create empty vector which is filled for each
                            %intensity again
                            pmu_tot = zeros(m,n);

                            %itereate over each subintensity
                            for kint = 1: n_signal_int
                                %select subintensity
                                intensity = decoys{k};
                                %calculate prob. distribution
                                pmu = probDist(intensity(kint),0:1:n_photon);
                                %create list with all subintensities
                                pmu_list{kint} = pmu;
                                pmu_tot(kint,:) = sum(pmu);
                            end
                            %create matrix with vector of pmu on diagonal
                            Pmu = blkdiag(pmu_list{:});
                            vec(Pmu*Y) >= vec(lowerbndsExpec(:,:,k) - (1-pmu_tot) - decoy_tolerance);
                            vec(Pmu*Y) <= vec(upperbndsExpec(:,:,k) + decoy_tolerance);                   
                        end

                        %Additional constraints for 0- and 1-photon components
                        %in terms of Choi matrices
                        for indexrow = 1:m
                            for indexcol = 1:n
                                norm(Y1(indexrow,indexcol) - trace(J1*kron(POVMsBob{indexcol},transpose(rhoAlice{indexrow})))) <= choi_tolerance;
                                norm(Y0(indexrow,indexcol) - trace(J0*kron(POVMsBob{indexcol},1))) <= choi_tolerance;
                            end
                        end
    
                        %Upper and lower bounds for n<=N_B subspace
                        if FlagStateSquash == true
                            for indexrow = 1:m
                                %1-photon
                                real(trace(J1*kron(Proj_sub,transpose(rhoAlice{indexrow})))) >= photonBound(indexrow,2) - choi_tolerance;                       
                                real(trace(J1*kron(Proj_sub,transpose(rhoAlice{indexrow})))) <= 1 + choi_tolerance;
                                %0-photon
                                real(trace(J0*kron(Proj_sub,1))) >= photonBound(indexrow,1) - choi_tolerance;
                                real(trace(J0*kron(Proj_sub,1))) <= 1 + choi_tolerance;
                            end
                        end
                        
                cvx_end

            if options.verboseLevel>=2
                disp(cvx_status)
            end

            if strcmp(cvx_status, 'Infeasible') | strcmp(cvx_status, 'Failed')
                fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
            end
            catch error 
              rethrow(error);
            end

            %Store upper bound
            Y1U(i+1,j) = Y(i*(n_photon+1)+2,j);            
        end
    end
            
    %solve for lower bounds
    for i=0:m-1
        for j=1:n
            try
                %Select the correct rows and columns from Y for the objective
                %function
                Objcolumn = zeros(n,1);
                Objcolumn(j) = 1;
                
                Objrow = zeros(1,m*(n_photon+1));
                Objrow(i*(n_photon+1)+2) = 1;

                %Set up SDP for upper bounds
                % Y is the matrix representing all yields
                % J0 is the Choi matrix for the 0-photon component
                % J1 is the Choi matrix for the 1-photon component
                % In both Choi matrices Bob's system is first and Alice's
                % system is second to satisfy e.g. Y_0 = Tr[J_0 (F_y otimes (rho_x)^T)] 
                cvx_begin quiet
                    %CVX solver
                    cvx_solver(convertStringsToChars(options.decoySolver));
                    variable Y(m*(n_photon+1),n) nonnegative
                    variable J0(dimB_SDP,dimB_SDP) hermitian semidefinite
                    variable J1(dimA_SDP*dimB_SDP,dimA_SDP*dimB_SDP) hermitian semidefinite
                    minimize( real(Objrow*Y*Objcolumn))
                    subject to
                        % 0 <= Y <=1
                        vec(Y) >= vec(zeros(m*(n_photon+1),n));
                        vec(Y) <= vec(ones(m*(n_photon+1),n));
                        
                        %0- and 1-photon component treated seperately with add.
                        %constraints
                        Y0 = M{1}*Y;
                        Y1 = M{2}*Y;
                        
                        %Partial trace of Choi matrix = id
                        norm(PartialTrace(J1,[1],[dimB_SDP,dimA_SDP]) - eye(dimA_SDP)) <= choi_tolerance;
                        norm(trace(J0) - 1) <= choi_tolerance;
                        
                        %Usual decoy bounds rewritten as matrix times vector
                        for k = 1:n_decoy
                            %total probability for given intensity to have
                            %<= n_photon photons
                            %create empty vector which is filled for each
                            %intensity again
                            pmu_tot = zeros(m,n);

                            %itereate over each subintensity
                            for kint = 1: n_signal_int
                                %select subintensity
                                intensity = decoys{k};
                                %calculate prob. distribution
                                pmu = probDist(intensity(kint),0:1:n_photon);
                                %create list with all subintensities
                                pmu_list{kint} = pmu;
                                pmu_tot(kint,:) = sum(pmu);
                            end
                            %create matrix with vector of pmu on diagonal
                            Pmu = blkdiag(pmu_list{:});
                            vec(Pmu*Y) >= vec(lowerbndsExpec(:,:,k) - (1-pmu_tot) - decoy_tolerance);
                            vec(Pmu*Y) <= vec(upperbndsExpec(:,:,k) + decoy_tolerance);                   
                        end

                        %Additional constraints for 0- and 1-photon components
                        %in terms of Choi matrices
                        for indexrow = 1:m
                            for indexcol = 1:n
                                norm(Y1(indexrow,indexcol) - trace(J1*kron(POVMsBob{indexcol},transpose(rhoAlice{indexrow})))) <= choi_tolerance;
                                norm(Y0(indexrow,indexcol) - trace(J0*kron(POVMsBob{indexcol},1))) <= choi_tolerance;
                            end
                        end
    
                        %Upper and lower bounds for n<=N_B subspace
                        if FlagStateSquash == true
                            for indexrow = 1:m
                                %1-photon
                                real(trace(J1*kron(Proj_sub,transpose(rhoAlice{indexrow})))) >= photonBound(indexrow,2) - choi_tolerance;                       
                                real(trace(J1*kron(Proj_sub,transpose(rhoAlice{indexrow})))) <= 1 + choi_tolerance;
                                %0-photon
                                real(trace(J0*kron(Proj_sub,1))) >= photonBound(indexrow,1) - choi_tolerance;
                                real(trace(J0*kron(Proj_sub,1))) <= 1 + choi_tolerance;
                            end
                        end
                        
                cvx_end

            if options.verboseLevel>=2
                disp(cvx_status)
            end

            if strcmp(cvx_status, 'Infeasible') | strcmp(cvx_status, 'Failed')
                fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
            end
            catch error 
              rethrow(error);%  fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
            end

            %Store upper bound
            Y1L(i+1,j) = Y(i*(n_photon+1)+2,j);            
        end
    end

    %make sure 0 <= Y1U <= 1 and replace nan with 1 (can be violated due to numerics) 
    Y1U(isnan(Y1U)) = 1;
    Y1U = min(max(Y1U,0),1);

    %make sure 0 <= Y1L <= 1 and replace nan with 1 (can be violated due to numerics) 
    Y1L(isnan(Y1L)) = 0;
    Y1L = min(max(Y1L,0),1);
end

%% validation functions
function SameNumberOfDecoysAndPages(decoys,conditionalExpectations)
if numel(decoys) ~= size(conditionalExpectations,3)
    throwAsCaller(MException("decoyAnalysisIndependentLP:DecoysAndPagesDontMatch",...
        "The number of decoy intensities does not match the number of pages for conditionalExpectations."))
end
end

function mustBeSameTotalProb(totProb,dist)
if ~ismembertol(sum(dist,"all"),totProb)
    throw(MException("decoyAnalysisConstrFinite:TotalTestProbViolated",...
        "The joint distribution sum to the probability of testing up to tolerance."))
end
end