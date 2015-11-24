function shell_images = segment_shells(image_data, centres, segment_half_size)

shell_images = zeros(length(centres) * (2 * segment_half_size + 1), 2 * segment_half_size + 1);

for i=1:length(centres)
	row_centre = floor(centres(i, 2));
	rows = (row_centre - segment_half_size):(row_centre + segment_half_size);
	col_centre = floor(centres(i, 1));
	cols = (col_centre - segment_half_size):(col_centre + segment_half_size);
	shell_image = image_data(rows, cols);
	shell_images((1 + (i - 1) * size(shell_image, 1)) : (i * size(shell_image, 1)), :) = shell_image;
end

end
