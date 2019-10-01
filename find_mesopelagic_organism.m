function [mo_features, x_spikes] = find_mesopelagic_organism(x, p, cfg, verbose)
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
%     x <Nx1 double> profile of backscattering, chlorophyll a fluorescence, or FDOM
%
% OPTIONAL INPUTS:
%     cfg <struct>
%     verbose <boolean> display classification information
%
% OUTPUTS:
%     mo <MX10 table> wiht the average depth, density, and intensity of each spike layer
%         computed from beta and fchl channels, official BGC-Argo floats parameters
%     x_spikes <Nx1 double> indexes of spikes in x profile
%
% author: Nils Haentjens
% created: Sept 16, 2019   

% Check input
if size(p,2) ~= 1; error('Input must be column vectors (Nx1).'); end
if size(p) ~= size(x); error('Input vectors are different size.'); end

% Default input
if nargin < 3
  cfg = struct('get_spikes', struct('xerr', [], 'pres_res', [], 'max_iter', []),...
               'cluster_spikes', struct('min_depth', [], 'distance_cut_off', []),...
               'classifier', [],...
               'min_obs', 30, 'min_depth', 50); end
if nargin < 4; verbose = false; end

% Check profile depth & observation number
if sum(~isnan(x) & ~isnan(p) & ~isinf(x) & ~isinf(p)) < cfg.min_obs ||...
   max(p) < cfg.min_depth
  mo_features = array2table(NaN(0,8), 'VariableNames',...
                    {'p', 'p_sd', 'n', 'p_shallow', 'p_deep', 'density', 'intensity', 'intensity_norm'});
  x_spikes = false(size(x));              
  return
end

% Find mesopelagic organisms
x_spikes = get_spikes(x, p, cfg.get_spikes.xerr, cfg.get_spikes.pres_res, cfg.get_spikes.max_iter);
x_clusters = cluster_spikes(x_spikes, p);
x_features = extract_features(x_clusters, p, x);
mo_present = classify_spike_clusters(x_features, cfg.classifier, verbose);
mo_features = x_features(mo_present,:);

end