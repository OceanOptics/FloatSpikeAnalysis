function [mo, varargout] = find_mesopelagic_organism(p, beta, fchl, fdom, cfg, verbose)
% FIND_MESOPELAGIC_ORGANISM use spikes in bio-optical profiles to find the depth of mesopelagic
%     organism (MO) if they are present.
% The method is combining an outlier detection algorithm (Hampel Filter) and a clustering algorithm
% (Hierarchical clusterin) to identify mesopelagic organism which are though to be at the origin of
% spikes in bio-optical profiles. Information from multiple profiles are used
%     1. Extract spikes (already done)
%     2. Cluster spikes
%     3. Compute features on clusters
%     4. Quality Check spike clusters
%     5. Combine QC-ed clusters from beta and fchl
%
% INPUTS (units do not matter):
%     p <Nx1 double> pressure of profile (add sampling frequency check)
%     beta <Nx1 double> profile of backscattering
%     fzchl <Nx1 double> profile of backscattering
%
% OPTIONAL INPUTS:
%     fdom <Nx1 double> profile of backscattering
%
% OUTPUTS:
%     mo <MX10 table> wiht the average depth, density, and intensity of each spike layer
%         computed from beta and fchl channels, official BGC-Argo floats parameters
%     mo_fdom <MX10 table> wiht the average depth, density, and intensity of each spike layer
%         computed from fdom channel, more reliable, but not available on every BGC-Argo floats
%         returned only if argument fdom was passed in
%     beta_spike, fchl_spike, fdom_spike are also returned
%
% author: Nils Haentjens
% created: Sept 16, 2019   

% Check if got FDOM
if nargin > 3 && ~isempty(fdom); flag_fdom = true;
else flag_fdom = false; end

% Check input
if size(p,2) ~= 1; error('Input must be column vectors (Nx1).'); end
if size(p) ~= size(beta); error('Input vectors are different size.'); end
if size(p) ~= size(fchl); error('Input vectors are different size.'); end
if flag_fdom; if size(p) ~= size(fdom); error('Input vectors are different size.'); end; end

% Default input
if nargin < 5; cfg = []; end
if nargin < 6; verbose = false; end

% Return empty array if profiles are too short
if size(p,1) < 5
  mo = array2table([NaN(0,10)], 'VariableNames',...
                    {'p', 'p_sd', 'n', 'p_shallow', 'p_deep', 'density', 'intensity', 'intensity_sd', 'intensity_norm', 'intensity_norm_sd'});
  empty_spike = false(size(beta));
  if flag_fdom;
    varargout = {mo, empty_spike, empty_spike, empty_spike};
  else
    varargout = {empty_spike, empty_spike};
  end                
  return
end

% Detect mesopelagic organism using beta and fchl channels
%   Spikes in beta channel are produice by any large particle, thereafter they are thought to originate from
%   either marine snow or mesopelagic organism. Mesopelagic organism are thought to generate layer
%   of spikes however this alone is not enought to distinguish them in beta profile when there is
%   many spikes due to marine snow spread accross the profile. To help us in the detection of MO we take advantage of the
%   fchl channel that register spikes in the presence of mesopelagic organism either due to there
%   gut fluorescence or due out of band response of the ECO sensor used. 

beta_spike = get_spikes(beta, p);
beta_cluster = cluster_spikes(beta_spike, p);
beta_features = extract_features(beta_cluster, p, beta);
beta_features_qc = classify_spike_cluster(beta_features, cfg, verbose);
fchl_spike = get_spikes(fchl, p);
% fchl_cluster = cluster_spikes(p, fchl_spike);
% fchl_features = compute_features(fchl_cluster, p, fchl);
% fchl_features_qc = qc_features(fchl_features, cfg, verbose);


beta_features_mo = false(height(beta_features),1);
% f.beta.mo(f.beta.qc) = combine_features(f.beta(f.beta.qc,:), f.fchl.p(f.fchl.qc).p);
beta_features_mo(beta_features_qc) = combine_features(beta_features(beta_features_qc,:), p(fchl_spike), 3);
mo = beta_features(beta_features_mo,:);

if flag_fdom
  % Detect MO using FDOM profiles only
  %     FDOM profiles does not seem to be sensitive to marine snow like the beta channel which makes it
  %     an ideal candidate to detect MO by itself.
  fdom_spike = get_spikes(fdom, p);
  fdom_cluster = cluster_spikes(fdom_spike, p);
  fdom_features = extract_features(fdom_cluster, p, fdom);
  fdom_features_qc = classify_spike_cluster(fdom_features, cfg, verbose);

  mo_fdom = fdom_features(fdom_features_qc,:);

  varargout = {mo_fdom, beta_spike, fchl_spike, fdom_spike};
else
  varargout = {beta_spike, fchl_spike};
end

end