function show_tiled_fitted_segments(shell_segments, fits)

num_segments = length(shell_segments);
num_rows = floor(sqrt(num_segments));
num_cols = ceil(num_segments / num_rows);
segment_size = size(shell_segments{1}, 1);

tiled_segments = zeros(num_rows * segment_size, num_cols * segment_size);

% Draw segments
for i=1:num_segments
	rows = (1 + (mod(i - 1, num_rows)) * segment_size):((mod(i - 1, num_rows) + 1) * segment_size);
	cols = (1 + floor((i - 1) / num_rows) * segment_size):(floor((i - 1) / num_rows + 1) * segment_size);
	tiled_segments(rows, cols) = shell_segments{i};
end
imshow(tiled_segments, [])

% Draw fit overlays
for i=1:length(fits)
	fit = fits{i};
	x_centre_fit = fit(1);
	y_centre_fit = fit(2);
	radius_fit = fit(3);

	% Correct for difference in coordinate basis and shift into right segment
	x_centre_fit = (x_centre_fit + floor((i - 1) / num_rows) * segment_size + segment_size / 2);
	y_centre_fit = (-y_centre_fit + mod(i - 1, num_rows) * segment_size + segment_size / 2);
	rectangle('Position', [x_centre_fit - radius_fit, y_centre_fit - radius_fit, radius_fit * 2, radius_fit * 2], 'Curvature', [1, 1], 'EdgeColor', 'red');
end

end
