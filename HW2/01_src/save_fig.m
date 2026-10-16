function save_fig(fh, name)
% Export a figure to 00_report/fig as vector PDF (for LaTeX) and PNG (for preview)
figdir = fullfile(fileparts(mfilename('fullpath')), '..', '00_report', 'fig');
if ~exist(figdir, 'dir'), mkdir(figdir); end

% Render numeric tick labels in math mode so negatives get a true minus sign
for ax = findobj(fh, 'Type', 'axes').'
    ax.XTickLabel = math_labels(ax.XTickLabel);
    ax.YTickLabel = math_labels(ax.YTickLabel);
end

exportgraphics(fh, fullfile(figdir, [name '.pdf']), 'ContentType', 'vector');
exportgraphics(fh, fullfile(figdir, [name '.png']), 'Resolution', 300);
end

function lbl = math_labels(lbl)
lbl = cellstr(lbl);
plain = ~contains(lbl, '$') & ~cellfun(@isempty, lbl);
lbl(plain) = strcat('$', lbl(plain), '$');
end
