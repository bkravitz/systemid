function mn = globalMean(val, cosLAT, timedim)
% function from Chris Jones
% return area-weighted global mean of val

if isempty(val)
    mn = NaN;
    return
end

if ~exist('cosLAT', 'var')
    % uniform weighting:
    mn = nanmean(val); 
    return
end

if ~exist('timedim', 'var')
    assert(all(size(val) == size(cosLAT)), 'globalMean:dimensionMismatch')
    sel = ~isnan(val) & ~isnan(cosLAT);
    mn = sum(cosLAT(sel).*val(sel))./sum(cosLAT(sel));
else
    % permute timeDim to first, the recursively call globalMean:
    perm = [timedim:max(ndims(val), timedim) 1:timedim-1];
    y = permute(val, perm);
    m = size(y,1);
    mn = NaN(m,1); % column vector
    for nt=m:-1:1
        mn(nt) = globalMean(squeeze(y(nt,:,:)), cosLAT);
    end
end
end