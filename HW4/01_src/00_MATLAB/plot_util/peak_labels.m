function peak_labels(mi, v, col, fmt, gap, dy)
% Mark (mi, v) points and label them "(m, value)" without overlap:
% a point whose right neighbour is within gap samples is labelled on its left,
% the other point of the pair on its right; isolated points are labelled on the left
if nargin < 4, fmt = '%.4f'; end
if nargin < 5, gap = 40; end
if nargin < 6, dy = 0; end
[mi, ord] = sort(mi(:));
v = v(:);
v = v(ord);
xl = xlim;
plot(mi, v, 'o', 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 6, 'HandleVisibility', 'off');
for k = 1:numel(mi)
    near_prev = k > 1 && mi(k) - mi(k-1) < gap;
    near_left_edge = mi(k) - xl(1) < 0.15*diff(xl);
    if near_prev || (near_left_edge && ~(k < numel(mi) && mi(k+1) - mi(k) < gap))
        ha = 'left'; dx = 0.008*diff(xl);
    else
        ha = 'right'; dx = -0.008*diff(xl);
    end
    text(mi(k) + dx, v(k) + dy, sprintf(['$(%d,\\,' fmt ')$'], mi(k), v(k)), 'Color', col, 'FontSize', 11, ...
        'HorizontalAlignment', ha, 'VerticalAlignment', 'bottom');
end
end
