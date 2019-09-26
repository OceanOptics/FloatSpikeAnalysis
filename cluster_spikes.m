function [idc] = cluster_spikes(spikes, p, DISTANCE_CUT_OFF)
%CLUSTER_SPIKES_V0 Cluster spikes when separated along the profile
%   use Hierarchical clustering
% 
% INPUTS (units do not matter):
%     p <Nx1 double> pressure of profile (add sampling frequency check)
%     spikes <Nx1 double|boolean> indices at which there is spikes.
%
% OPTIONAL INPUTS:
%     DISTANCE_CUT_OFF <Nx1> euclidean distance used in clustering algorithm to cut tree
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

if nargin < 2
  DISTANCE_CUT_OFF = 50; % dBar
end

if size(p) ~= size(spikes); error('Input vectors are different size.'); end
if size(p,2) ~= 1; error('Input must be row vectors.'); end

if sum(spikes) < 2
  idc = NaN(size(p));
  idc(spikes) = 1;
  return 
end

D = pdist(p(spikes),'euclidean');
ClusterTree = linkage(D,'single');
idx = cluster(ClusterTree,'criterion','distance','cutoff',50);

idc = NaN(size(p));
idc(spikes) = idx;

end

