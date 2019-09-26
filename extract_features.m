function features = extract_features(idc, p, x)
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
  features = array2table([NaN(n,11)], 'VariableNames',...
                    {'p', 'p_sd', 'n', 'p_shallow', 'p_deep', 'n_x_neg', 'density', 'intensity', 'intensity_sd', 'intensity_norm', 'intensity_norm_sd'});
else
  features = array2table([NaN(n,7)], 'VariableNames',...
                    {'p', 'p_sd', 'n', 'p_shallow', 'p_deep', 'n_x_neg', 'density'});
end


x0 = x - min(x); % Adjust values to make sure positive signal
% x_snr_floor = snr(x(~isnan(x) & isnan(idc)));
for i=1:n
  sel = idc == i;
  if isempty(sel); continue; end
  
  features.p(i) = nanmedian(p(sel));
  features.p_sd(i) = nanstd(p(sel));
  features.n(i) = sum(sel);
  features.p_shallow(i) = min(p(sel));
  features.p_deep(i) = max(p(sel));
  
  features.density(i) = features.n(i) / (features.p_deep(i) - features.p_shallow(i));
  
  features.n_x_neg(i) = sum(x(features.p_shallow(i) <= p & p <= features.p_deep(i)) < 0);
  
%   % Signal to Noise Ratio
%   snr_sel = features.p_shallow(i) <= p & p <= features.p_deep(i) &...
%                             ~isnan(x) & isnan(idc);
%   if any(snr_sel)
%     features.snr_x(i) = snr(x(snr_sel));
%     features.snr_x_norm(i) = features.snr_x(i) - x_snr_floor;
%   else
%     features.snr_x(i) = NaN;
%     features.snr_x_norm(i) = NaN;
%   end
  
  if nargin > 2
    % Compute absoulte intensity
    features.intensity(i) = nanmedian(x(sel));
    features.intensity_sd(i) = nanstd(x(sel));
    % Compute relative intensity
    si = x0(sel) ./ nanmedian(x0); % Normalize
    features.intensity_norm(i) = nanmedian(si);
    features.intensity_norm_sd(i) = nanstd(si);    
  end
end


end