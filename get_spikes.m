function [spikes, varargout] = get_spikes(x, p, xerr, pr, max_iter)
%SPIKE_ANALYSIS a point is considered as a spike if the difference between x and a 15 points Hampel
%       filter of x is greater than 3 times the scaled Median Absolute Deviation
%       also return negative spikes separately
%
% INPUTS (units do not matter):
%     x <Nx1 double> profile of bio-optical data, for example:
%         - particulate backscattering (beta, bbp)
%         - fluorescent dissolved organic matter (fdom)
%         - chlorophyll a fluorescence  (fchl)
% OPTIONAL INPUT:
%     p <Nx1 double> pressure of profile (add sampling frequency check)
%     xerr <1x1 double> absolute threshold for spike detection
%         default: 3x the scaled Median Absolute Deviation (MAD) of entire profile (non-parametric)
%         recommended: for fdom: ???; for bbp: ???; for fchl: ???;                 (parametric)
%     pr <1x1 double> minimum sampling resolution threshold
%         default: 3x the median sampling resolution of entire profile             (non-parametric)
%         recommended: 5 dBar (typically no spike layers are found below 2.25 dBar)(parametric)
%     max_iter <1x1 integer> maximum number of intereation to run Hampel Filter
%         set to 0 to run the Hampel filter only once
%         default: 100
%     
%
% OUTPUTS:
%     spikes <NX1 boolean> return true for every indices containing a positive spike
%     neg_spikes <NX1 boolean> return true for every indices containing a negative spike
%     spikes_p <NX1 boolean> pressure of each positive spikes
%     neg_spikes_p <NX1 boolean> pressure of each negative spikes
%
% References:
%   Rousseeuw, P. J., and C. Croux (1993), Alternatives to the Median Absolute Deviation,
%     Journal of the American Statistical Association, 88(424), 1273?1283.
%
% author: Nils Haentjens
% created: May 1, 2018

DEBUG = false;

% Check input
if nargin < 2 || isempty(p); p_flag = false; else; p_flag = true; end
if nargin < 3 || isempty(xerr); xerr_flag = false; xerr = NaN; else; xerr_flag = true; end
if nargin < 4 || isempty(pr); pr_flag = false; else; pr_flag = true; end % ignored if ~p_flag
if nargin < 5 || isempty(max_iter); max_iter = 100; end

if size(x,2) ~= 1; error('x is not a row vector (Nx1).'); end
if p_flag && any(size(p) ~= size(x)); error('Input vectors are different size.'); end

% Check output
if nargout > 2 || (nargout == 2 && ~p_flag); neg_spike_flag = true; else; neg_spike_flag = false; end

% Live debug
if DEBUG
  fig(1);
  subplot(1,3,2); hold('on');
  plot(x,p,'.-'); 
end

% Prepare input
spikes = false(size(x));
if p_flag
  % ignore nan of x and p for processing
  sel = find(~isnan(x) & ~isnan(p) & ~isinf(x) & ~isinf(p));
  if length(sel) < 3
    warning('GET_SPIKES:EMPTY_INPUT', 'get_spikes: Not enough valid values.');
    if nargout == 2 && ~p_flag; varargout = {spikes}; else; varargout = {[], spikes, []}; end
    return;
  end
  p = p(sel);
  % ignore low sampling resolution of profile
  delta = abs([diff(p(1:2)); (p(3:end) - p(1:end-2)) / 2; diff(p(end-1:end))]); % # / dBar
  if ~pr_flag; pr = 3 * median(delta); end
  subsel = delta < pr;
  sel = sel(subsel);
  if DEBUG
    subplot(1,3,1);
    plot(delta, p, '.-');
    set(gca, 'ydir', 'reverse');
  end
  % keep only relevant sampling
  x = x(sel); p = p(subsel);
else
  % only x
  sel = find(~isnan(x) & ~isinf(x));
  x = x(sel);
end

% Apply Hampel filter (same as median filter except onlx change values of spikes)
xd = hampel(x, 15); % > 5 is good
% Recursive Hampel filter
xdo = x; i = 0;
while any(xd ~= xdo) && i < max_iter
  xdo = xd; 
  xd = hampel(xd, 15);
  i = i + 1;
end
if i == max_iter; warning('GET_SPIKES:MAX_ITER', 'get_spikes: Maximum iteration (%d) reached.', max_iter); end

if ~xerr_flag
  % Compute Scaled Median Absolute Deviation (MAD)
  smad = -1/(sqrt(2)*erfcinv(3/2)) * median(abs(x-median(x)));
  % Noisy signal (no significant spike)
  if smad == 0
    % format output & return
    if nargout == 2 && ~p_flag; varargout = {spikes}; else; varargout = {[], spikes, []}; end
    return;
   end
  % Set spike detection threshold to 3 scaled MAD
  xerr = 3 * smad;
end

% Get spikes
x_delta = x - xd;
spikes(sel) = x_delta > xerr;

% Live debug
if DEBUG
  subplot(1,3,2);
  plot(xd, p, '.-');
  plot(x(spikes(sel)), p(spikes(sel)), 'o');
  set(gca, 'ydir', 'reverse');
  y_lim = ylim();
  subplot(1,3,3); hold('on');
  plot(x_delta, p, '.-');
  plot(xerr .* ones(2,1), ylim(), 'k', 'LineWidth', 2);
  ylim(y_lim);
  set(gca, 'ydir', 'reverse');
  xlim([0 0.01]);
%   set(gca, 'xscale', 'log');
  drawnow();
  waitforbuttonpress();
%   pause(1);
end

% Get negative spikes
if neg_spike_flag
  neg_spikes = false(size(spikes));
  neg_spikes(sel) = -x_delta > xerr;
end

% Set output
if p_flag
  switch nargout
    case 2
      varargout = {p(spikes(sel))};
    case 3
      varargout = {p(spikes(sel)), neg_spikes};
    case 4
      varargout = {p(spikes(sel)), neg_spikes, p(neg_spikes(sel))};
  end
else
  varargout = {neg_spikes};
end

end