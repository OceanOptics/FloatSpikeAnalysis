function features = compute_features(idc, p, x)
% Compute features for each cluster of spikes within a profile
%   cluster ids should be consecutive numbers starting from 1
%
% author: Nils Haentjens
% created: Sept 12, 2019



if size(idc) ~= size(p); error('Input vectors are different size.'); end
if nargin > 2
  if size(idc) ~= size(x); error('Input vectors are different size.'); end
end

n = max(idc); if isnan(n); n=0; end
if nargin > 2
  features = array2table([NaN(n,10)], 'VariableNames',...
                    {'p', 'p_sd', 'n', 'p_shallow', 'p_deep', 'density', 'intensity', 'intensity_sd', 'intensity_norm', 'intensity_norm_sd'});
else
  features = array2table([NaN(n,6)], 'VariableNames',...
                    {'p', 'p_sd', 'n', 'p_shallow', 'p_deep', 'density'});
end

for i=1:n
  sel = idc == i;
  if isempty(sel); continue; end
  
  features.p(i) = nanmedian(p(sel));
  features.p_sd(i) = nanstd(p(sel));
  features.n(i) = sum(sel);
  features.p_shallow(i) = min(p(sel));
  features.p_deep(i) = max(p(sel));
  
  features.density(i) = features.n(i) / (features.p_deep(i) - features.p_shallow(i));
  
  if nargin > 2
    % Compute absoulte intensity
    features.intensity(i) = nanmedian(x(sel));
    features.intensity_sd(i) = nanstd(x(sel));
    % Compute relative intensity
    x = x - min(x); % Adjust values to make sure positive signal
    si = x(sel) ./ nanmedian(x); % Normalize
    features.intensity_norm(i) = nanmedian(si);
    features.intensity_norm_sd(i) = nanstd(si);
  end
end


end