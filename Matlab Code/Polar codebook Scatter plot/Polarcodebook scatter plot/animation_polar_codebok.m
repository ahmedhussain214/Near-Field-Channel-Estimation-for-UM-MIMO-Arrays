clc;
clear all;
close all;

% Parameters
fc = 100e9; % Carrier frequency
c = 3e8;
lambda = c/fc;
d = lambda/2;
N = 256;
D = (N-1)*d;
RD = (2*D^2)/lambda;

theta0 = linspace(-pi/2, pi/2, N/8);
label = [];

% Compute sampling points
for i = 1:length(theta0)
    F = 2*D;
    theta = theta0(i);
    rr = [];

    while F <= (RD/10)*(cos(theta)^2)
        r = RD*cos(theta)^2 * F * (1/(RD*cos(theta)^2 - 10*F) - 1/(RD*cos(theta)^2 + 10*F));
        rr = [rr F];
        F = F + r;
    end

t =(RD/10)*(cos(theta).^2);
if t>2*D
    v = t ;
else
    v =[];
end
rr =[rr   v];

    label1 = [ones(1, length(rr)) * theta; rr];
    label = [label label1];
end

% Restructure data for animation: group by unique distances
all_distances = unique(round(label(2, :), 6)); % Use rounded values to avoid numerical precision issues
sorted_angles = sort(unique(label(1, :)));     % Ensure angle consistency

% Create a mapping of distance -> list of (theta, r)
plot_sequence = [];
for r_idx = 1:length(all_distances)
    r_val = all_distances(r_idx);
    for a_idx = 1:length(sorted_angles)
        angle = sorted_angles(a_idx);
        % Find the index of this angle-distance pair
        match = find(abs(label(1,:) - angle) < 1e-6 & abs(label(2,:) - r_val) < 1e-6, 1);
        if ~isempty(match)
            plot_sequence = [plot_sequence, match];
        end
    end
end

% Assign one color per angle
num_angles = length(sorted_angles);
colors = jet(num_angles);
angle_color_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
for i = 1:num_angles
    angle_color_map(sorted_angles(i)) = colors(i, :);
end
point_colors = zeros(size(label, 2), 3);
for i = 1:size(label, 2)
    point_colors(i, :) = angle_color_map(label(1, i));
end

% Plot setup
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
fig = figure('Position', [100, 100, 800, 600]);
ax = gca;
colormap(ax, colors);
xlabel('Angle', 'Interpreter', 'Latex');
ylabel('Distance', 'Interpreter', 'Latex');
subtitle({'Codewords spaced by beamwidth in angle,', ...
          'and \textbf{beamdepth} in range', 'provide 3dB coverage and low correlation'}, ...
          'Interpreter', 'latex');

yline(RD, '--', {'Rayleigh distance'}, 'fontsize', 14, 'Interpreter', 'Latex');
yline(RD/10, '--', {'Rayleigh distance /10'}, 'fontsize', 14, 'Interpreter', 'Latex');
yticks([2*D floor(RD/10) floor(RD)]);
yticklabels({'$2D$', '$\frac{r_{RD}}{10}$', '$r_{RD}$'});
yscale log;
ylim([2*D RD/10]);
xlim([-pi/2.4 pi/2.4]);
xticks([-pi/2.4 0 pi/2.4]);
xticklabels({'$-\pi/2$', '$0$', '$\pi/2$'});
set(gca, 'fontsize', 18);
hold on;

% Animation
filename = 'polar_codebook_distance_first.gif';
scatter_handle = scatter([], [], 45, [], 'filled');
frame_delay_regular = 0.0001;
frame_delay_first_frame = 0.00003;

for i = 1:length(plot_sequence)
    idx = plot_sequence(i);
    set(scatter_handle, ...
        'XData', label(1, plot_sequence(1:i)), ...
        'YData', label(2, plot_sequence(1:i)), ...
        'CData', point_colors(plot_sequence(1:i), :));
    drawnow;

    frame = getframe(fig);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);

    if i == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', frame_delay_first_frame);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', frame_delay_regular);
    end
end
