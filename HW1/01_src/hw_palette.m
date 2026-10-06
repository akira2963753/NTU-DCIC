function p = hw_palette()
% Shared plot colors for HW1 figures
p.blue = hex2rgb('#2a78d6');      % I / real part
p.orange = hex2rgb('#eb6834');    % Q / imaginary part
p.aqua = hex2rgb('#1baf7a');      % passband
p.red = hex2rgb('#e34948');       % highlight
p.ink = hex2rgb('#0b0b0b');
p.ink2 = hex2rgb('#52514e');
p.muted = hex2rgb('#898781');
p.grid = hex2rgb('#e1e0d9');
end

function c = hex2rgb(h)
c = sscanf(h(2:end), '%2x%2x%2x', [1 3])/255;
end
