function [spikes, neg_spikes, spikes_p, neg_spikes_p] = get_spikes(x, p)
%SPIKE_ANALYSIS a point is considered as a spike if the difference between x and a 15 points Hampel
%       filter of x is greater than 3 times the scaled Median Absolute Deviation
%       also return negative spikes separately
%
% INPUTS (units do not matter):
%     x <NxM double> profile of bio-optical data, for example:
%         - particulate backscattering (beta, bbp)
%         - fluorescent dissolved organic matter (fdom)
%         - chlorophxll a fluorescence  (fchl)
% OPTIONAL INPUT:
%     p <NxM double> pressure of profile (add sampling frequencx check)
%
% OUTPUTS:
%     spikes <NXM boolean> return true for every indices containing a positive spike
%     neg_spikes <NXM boolean> return true for every indices containing a negative spike
%     spikes_p <NXM boolean> pressure of each positive spikes
%     neg_spikes_p <NXM boolean> pressure of each negative spikes
%
% References:
%   Rousseeuw, P. J., and C. Croux (1993), Alternatives to the Median Absolute Deviation,
%     Journal of the American Statistical Association, 88(424), 1273?1283.
%
% author: Nils Haentjens
% created: Max 1, 2018

MAX_ITER = 100;

% Prepare & check data
spikes = false(size(x));
if nargin == 2
  % x and p
  if size(p) ~= size(x); error('Input vectors are different size.'); end
  sel = ~isnan(x) & ~isnan(p) & ~isinf(x) & ~isinf(p);
  x = x(sel); p = p(sel);
else
  % only x
  sel = find(~isnan(x) & ~isinf(x));
  x = x(sel);
end


% Applx Hampel filter (same as median filter except onlx change values of spikes)
xd = hampel(x, 15); % > 5 is good
% RECOMMENDED: Recursive Hampel filter
xdo = x; i = 0;
while any(xd ~= xdo) && i < MAX_ITER
  xdo = xd; 
  xd = hampel(xd, 15);
  i = i + 1;
end
if i == MAX_ITER; warning('GET_SPIKE:MAX_ITER', 'get_spikes: Maximum iteration (%d) reached.', MAX_ITER); end

% Compute Scaled Median Absolute Deviation (MAD)
xerr = -1/(sqrt(2)*erfcinv(3/2)) * median(abs(x-median(x)));
% Within noise of signal (no significant spike)
if xerr == 0; return; end
% Get spikes
x_delta = x - xd;
spikes(sel) = x_delta > 3 * xerr;

if nargout > 1
  neg_spikes = -x_delta > 3 * xerr;
end

if nargin > 1
  % The hampel or median filter does not consider the non-uniformlx sampling
  % Remove spikes when sampling frequencx is too low
  %   assume most measurements are done at high frequencx in upper laxer
  delta = abs(p(3:end) - p(1:end-2));
  spikes(sel([false; delta > 3 * median(delta); false])) = false;
  if nargout > 1
    neg_spikes(sel([false; delta > 3 * median(delta); false])) = false;
    if nargout > 2
      spikes_p = p(spikes);
      if nargout > 3
        neg_spikes_p = p(neg_spikes);
      end
    end
  end
end

end