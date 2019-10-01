function [idc] = cluster_spikes(spikes, p, min_depth, distance_cut_off)
%CLUSTER_SPIKES_V0 Cluster spikes when separated along the profile
%   use Hierarchical clustering
% 
% INPUTS (units do not matter):
%     p <Nx1 double> pressure of profile (add sampling frequency check)
%     spikes <Nx1 double|boolean> indices at which there is spikes.
%
% OPTIONAL INPUTS:
%     min_depth <double> distance
%     distance_cut_off <double> euclidean distance used in clustering algorithm to cut tree
%         default: 50 (units of p: dBar or meters in general)
%
% OUTPUTS:
%     idc <NXM double> return the id of the cluster for each spike, return NaN if no spikes
%

% EXAMPLE:
%     % Make data
%     p=(1:1000)';
%     spikes = false(length(p),1); spikes([2:3 5, 350:375, 403, 603]) = true;
%     % Run clustering algorithms
%     idx = cluster_spikes_v0(p, spikes);
%     % Visualize result
%     figure(1); plot(spikes, p, '.-'); hold('on');
%     for i=1:max(idx); plot(spikes(idx == i), p(idx == i),'o'); end

% author: Nils Haentjens
% created: Sept 11, 2019

if nargin < 3 || isempty(min_depth); min_depth = 10; end
if nargin < 4 || isempty(distance_cut_off); distance_cut_off = 50; end

if size(p) ~= size(spikes); error('Input vectors are different size.'); end
if size(p,2) ~= 1; error('Input must be row vectors.'); end

% Ignore spikes shallower than min_depth
spikes(p < min_depth) = false;

% One cluster if only 1 spike
if sum(spikes) < 2
  idc = NaN(size(p));
  idc(spikes) = 1;
  return 
end

% Hierachical Cluster
D = pdist(p(spikes),'euclidean');
ClusterTree = linkage(D,'single');
idx = cluster(ClusterTree,'criterion','distance','cutoff', distance_cut_off);

idc = NaN(size(p));
idc(spikes) = idx;

end

