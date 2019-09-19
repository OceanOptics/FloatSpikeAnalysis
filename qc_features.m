function qc = qc_features(features, cfg, verbose)
% ClEAN_SPIKES return indexes of spike clusters which potentially correspond to mesopelagic organism

if nargin < 2 || isempty(cfg)
  cfg.N_THRESHOLD = 3; 
  cfg.MIN_DEPTH = 10; 
  cfg.MAX_DEPTH = 1000; 
  cfg.MAX_THICKNESS = 350;
  cfg.MIN_REL_INTENSITY = 1;
end

if nargin < 3; verbose = false; end

qc = true(height(features),1);
for i=1:height(features)
  if features.n(i) < cfg.N_THRESHOLD
    if verbose; fprintf('%2d: N_THRESHOLD test FAILED\n', i); end
    qc(i) = false; continue;
  end
  
  if features.p(i) < cfg.MIN_DEPTH
    if verbose; fprintf('%2d: MIN_DEPTH test FAILED\n', i);end
    qc(i) = false; continue;
  end
  if features.p(i) > cfg.MAX_DEPTH
    if verbose; fprintf('%2d: MAX_DEPTH test FAILED\n', i); end
    qc(i) = false; continue;
  end
  
  if features.p_deep(i) - features.p_shallow(i) > cfg.MAX_THICKNESS
    if verbose; fprintf('%2d: MAX_THICKNESS test FAILED\n', i); end
    qc(i) = false; continue;
  end
  
  if features.intensity_norm(i) < cfg.MIN_REL_INTENSITY
    if verbose; fprintf('%2d: MIN_REL_INTENSITY test FAILED\n', i); end
    qc(i) = false; continue;
  end
end


end