function [det, pkv] = sbd_peak_search(phi, eta, W)
% Threshold-gated local peak searcher, sample by sample (same rule as PEAK_SEARCH.sv)
%   - a search starts when phi > eta (and the previous search region has dropped below eta)
%   - the candidate is replaced only by a strictly larger value that is also above eta
%   - the candidate is declared W samples after it was found without a larger value,
%     or at the last sample of the record (flush)
%   - after a declaration a new search needs phi <= eta first (one detection per lobe)
M = numel(phi);
det = [];
pkv = [];
run = false; hold = false;
pv = 0; pi_ = 0; cnt = 0;
for m = 1:M
    above = phi(m) > eta;
    if ~run
        if above && ~hold
            run = true; pv = phi(m); pi_ = m; cnt = 0;
        end
        if ~above, hold = false; end
    else
        if above && phi(m) > pv
            pv = phi(m); pi_ = m; cnt = 0;
        else
            cnt = cnt + 1;
        end
    end
    last = m == M;
    if run && (cnt == W || last)
        det(end + 1, 1) = pi_; %#ok<AGROW>
        pkv(end + 1, 1) = pv; %#ok<AGROW>
        run = false;
        hold = above;
    end
end
end
