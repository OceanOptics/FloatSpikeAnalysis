function [qc, alpha] = combine_features(c1_features, c2_p, THRESHOLD)
% FIND_MESO_ORGA look for matching features in two channels and keep them
%   THRESHOLD can be used to define the minimum number of events in c2.
%   return qc for features_c1

if nargin < 3; THRESHOLD = 1; end

alpha = zeros(height(c1_features),1);
qc=false(height(c1_features),1);

if isempty(c2_p); return; end
if ~any(c2_p); return; end

for i=1:height(c1_features)
%   alpha(i) = sum(abs(features_c1.p(i) - features_c2.p) < max(features_c1.p_sd(i), features_c2.p_sd));
  alpha(i) = sum(c1_features.p_shallow(i) < c2_p & c2_p < c1_features.p_deep(i));
  qc(i) = alpha(i) >= THRESHOLD;
end

end