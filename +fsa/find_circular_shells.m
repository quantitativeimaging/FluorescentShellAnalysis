function [centres, radii, metric] = find_circular_shells(image_data, radius_lower, radius_upper, segment_half_size, edge_border, ShowPlot);

[centres, radii, metric] = imfindcircles(image_data, [radius_lower radius_upper]);


% Remove candidates near edge
A = [centres, radii, metric];
A(A(:, 1) < segment_half_size + edge_border + 1, :) = [];
A(A(:, 1) > size(image_data, 2) - (segment_half_size + edge_border) - 1, :) = [];
A(A(:, 2) < segment_half_size + edge_border + 1, :) = [];
A(A(:, 2) > size(image_data, 1) - (segment_half_size + edge_border) - 1, :) = [];


if (ShowPlot)
	figure
	imshow(image_data, []);
	colormap(gray)
	truesize;
	% hold on
	% scatter(centres(:,1),centres(:,2), pi * radii.^2, 'co', 'lineWidth', 2)
	hold on
	for (i=1:length(centres))
		x = centres(i, 1) - segment_half_size;
		y = centres(i, 2) - segment_half_size;
		rectangle('Position', [x, y, segment_half_size*2, segment_half_size*2], 'EdgeColor', 'r')
	end
	hold on
	rectangle('Position', [edge_border, edge_border, size(image_data, 2) - 2*edge_border, size(image_data, 1) - 2*edge_border], 'EdgeColor', 'g')
	axis equal
end

end
