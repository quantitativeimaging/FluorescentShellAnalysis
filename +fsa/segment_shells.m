function shell_images = segment_shells(image_data, centres, segment_half_size)

shell_images = cell(length(centres), 1);

for i=1:length(centres)
	row_centre = floor(centres(i, 2));
	rows = (row_centre - segment_half_size):(row_centre + segment_half_size);
	col_centre = floor(centres(i, 1));
	cols = (col_centre - segment_half_size):(col_centre + segment_half_size);
	shell_image = image_data(rows, cols);
	shell_images(i) = {shell_image};
end

end
