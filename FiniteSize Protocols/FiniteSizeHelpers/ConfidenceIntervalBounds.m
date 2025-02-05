function [lowerBound, upperBound] = ConfidenceIntervalBounds(list,N,alpha)
% List is a list of observed propobailities of events 
% N is total number of trials
% epsAT is the security parameter

shape = size(list);


[~,bounds] = binofit(round(N*list(:)),N,alpha); % we round to nearest integer

lowerBounds = bounds(:,1);
upperBounds = bounds(:,2);

lowerBound  = reshape(lowerBounds, shape);
upperBound = reshape(upperBounds, shape);
end
