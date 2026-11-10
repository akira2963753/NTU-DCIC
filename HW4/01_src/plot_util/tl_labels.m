function tl_labels(tl, ttl, xl)
% Title and shared x label of a tiled layout, in LaTeX
if ~isempty(ttl), title(tl, ttl, 'Interpreter', 'latex', 'FontSize', 15); end
if ~isempty(xl), xlabel(tl, xl, 'Interpreter', 'latex', 'FontSize', 14); end
end
