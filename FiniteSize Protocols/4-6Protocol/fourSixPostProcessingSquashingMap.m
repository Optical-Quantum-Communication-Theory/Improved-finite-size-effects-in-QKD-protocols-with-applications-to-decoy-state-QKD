function mapping = fourSixPostProcessingSquashingMap()
%Passive Squashing Map for 4-6 protocol

function mapping = quickMap(mapping,pattern,remaping)
    mapping(:,sub2indPlus(2*ones(1,numel(pattern)),pattern+1)) = remaping;
end

mapping = zeros(11,64);

% A^T = Lambda * E^T => A = E*Lambda^T
% => Lambda 8 by 16 since E 4x16 and A 4x11

% The vast majority of the squashed bits are cross clicks that are mapped
% to vac for discarding. We will replace patterns that don't represent
% cross clicks in later steps.
%mapping(5,:) = 1;

%order S_( H V D A R L) D_( HV DA RL ) CC_(ANY) NON
% zero clicks to vac
mapping = quickMap(mapping,[0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,1]);

% single clicks to single clicks
mapping = quickMap(mapping,[1,0,0,0,0,0],[1,0,0,0,0,0,0,0,0,0,0]); % H
mapping = quickMap(mapping,[0,1,0,0,0,0],[0,1,0,0,0,0,0,0,0,0,0]); % V
mapping = quickMap(mapping,[0,0,1,0,0,0],[0,0,1,0,0,0,0,0,0,0,0]); % D
mapping = quickMap(mapping,[0,0,0,1,0,0],[0,0,0,1,0,0,0,0,0,0,0]); % A
mapping = quickMap(mapping,[0,0,0,0,1,0],[0,0,0,0,1,0,0,0,0,0,0]); % R
mapping = quickMap(mapping,[0,0,0,0,0,1],[0,0,0,0,0,1,0,0,0,0,0]); % L

% double clicks
mapping = quickMap(mapping,[1,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0,0]); % Z (HV)
mapping = quickMap(mapping,[0,0,1,1,0,0],[0,0,0,0,0,0,0,1,0,0,0]); % X (DA)
mapping = quickMap(mapping,[0,0,0,0,1,1],[0,0,0,0,0,0,0,0,1,0,0]); % Y (RL)

%cross clicks
sameBasisClicks = [[1,1,0,0,0,0]; [0,0,1,1,0,0]; [0,0,0,0,1,1]];
twoClicks = [1,1,0,0,0,0];
permuteTwo = perms(twoClicks);
previousPatterns = [0,0,0,0,0,0]; 
for j = 1:size(permuteTwo, 1)
    clickPattern = permuteTwo(j, :);
    if ~ismember(clickPattern, sameBasisClicks, 'rows') && ~ismember(clickPattern, previousPatterns, 'rows')%ensure not double click or previously chosen
        previousPatterns = [previousPatterns; clickPattern];
        mapping = quickMap(mapping,clickPattern,[0,0,0,0,0,0,0,0,0,1,0]);
    end
end
threeClicks = [1,1,1,0,0,0];
permuteThree = perms(threeClicks); 
previousPatterns = [0,0,0,0,0,0]; 
for j = 1:size(permuteTwo, 1)
    clickPattern = permuteThree(j, :);
    if ~ismember(clickPattern, previousPatterns, 'rows')%ensure not previously chosen
        previousPatterns = [previousPatterns; clickPattern];
        mapping = quickMap(mapping,clickPattern,[0,0,0,0,0,0,0,0,0,1,0]);
    end
end 
fourClicks = [1,1,1,1,0,0];
permuteFour = perms(fourClicks); 
previousPatterns = [0,0,0,0,0,0]; 
for j = 1:size(permuteTwo, 1)
    clickPattern = permuteFour(j, :);
    if ~ismember(clickPattern, previousPatterns, 'rows')%ensure not previously chosen
        previousPatterns = [previousPatterns; clickPattern];
        mapping = quickMap(mapping,clickPattern,[0,0,0,0,0,0,0,0,0,1,0]);
    end
end
fiveClicks = [1,1,1,1,1,0];
permuteFive = perms(fiveClicks); 
previousPatterns = [0,0,0,0,0,0]; 
for j = 1:size(permuteTwo, 1)
    clickPattern = permuteFive(j, :);
    if ~ismember(clickPattern, previousPatterns, 'rows')%ensure not double click or previously chosen
        previousPatterns = [previousPatterns; clickPattern];
        mapping = quickMap(mapping,clickPattern,[0,0,0,0,0,0,0,0,0,1,0]);
    end
end
mapping = quickMap(mapping,[1,1,1,1,1,1],[0,0,0,0,0,0,0,0,0,1,0]);

end