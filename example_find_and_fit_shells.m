image_data = imread('testData_megaterium_spheres.tif');

% Set the following flag to 1 to see each segment as it is fitted
SHOW_ALL_FITS = 0;

% Parameters for shell finding
radius_lower = 5;
radius_upper = 15;
segment_half_size = 20;
edge_border = 0;

% Find and display circular shells
[centres, radii, metric] = fsa.find_circular_shells(image_data, radius_lower, radius_upper, segment_half_size, edge_border, true);
title('Identified circular shells')

% Tile segmented shells in figure
shell_segments = fsa.segment_shells(image_data, centres, segment_half_size);
tiled_segments = fsa.tile_segments(shell_segments);
figure
imshow(tiled_segments, [])
title('All segmented shell images')

% Fit all segmented shells and display to the user one by one, if flag set
fits = cell(length(shell_segments), 1);
if (SHOW_ALL_FITS)
	figure
	for i=1:length(shell_segments)
		actual_image = shell_segments{i};
		background = median(actual_image(actual_image < mean(actual_image(:))));
		actual_image = double(actual_image - background);

		% Initial guess for shell parameters
		x_centre = 0; y_centre = 0; radius = 12; psf_sigma = sqrt(2); height = max(actual_image(:));

		% Fit shell to spore segment
		[x_centre_fit, y_centre_fit, radius_fit, psf_sigma_fit, height_fit] = fsa.fit_sphere_thin(x_centre, y_centre, radius, psf_sigma, height, actual_image);

		fits{i} = [x_centre_fit, y_centre_fit, radius_fit, psf_sigma_fit, height_fit];

		% Display segmented spore, with red fit overlay
		subplot(1, 3, 1)
		imshow(actual_image, [])
		hold on
		rectangle('Position', [size(actual_image, 2) / 2 + x_centre_fit - radius_fit,  size(actual_image, 1) / 2 - y_centre_fit - radius_fit, radius_fit * 2, radius_fit * 2], 'Curvature', [1, 1], 'EdgeColor', 'red');
		xlabel('Segmented + radial fit overlay')

		% Display initial guess
		subplot(1, 3, 2)
		imshow(fsa.image_sphere_thin(x_centre, y_centre, radius, psf_sigma, height, actual_image), [])
		title(sprintf('Final fit params: x = %f, y = %f, r = %f, s = %f, h = %f', x_centre_fit, y_centre_fit, radius_fit, psf_sigma_fit, height_fit))
		xlabel('Initial guess')

		% Display final fit
		subplot(1, 3, 3)
		imshow(fsa.image_sphere_thin(x_centre_fit, y_centre_fit, radius_fit, psf_sigma_fit, height_fit, actual_image), [])
		xlabel('Final fit')

		input(sprintf('Displaying segment %d of %d. Press return to continue...', i, length(shell_segments)))
	end
else
	parfor i=1:length(shell_segments)
		actual_image = shell_segments{i};
		background = median(actual_image(actual_image < mean(actual_image(:))));
		actual_image = double(actual_image - background);

		% Initial guess for shell parameters
		x_centre = 0; y_centre = 0; radius = 12; psf_sigma = sqrt(2); height = max(actual_image(:));

		% Fit shell to spore segment
		[x_centre_fit, y_centre_fit, radius_fit, psf_sigma_fit, height_fit] = fsa.fit_sphere_thin(x_centre, y_centre, radius, psf_sigma, height, actual_image);

		fits{i} = [x_centre_fit, y_centre_fit, radius_fit, psf_sigma_fit, height_fit];
	end
end

figure
fsa.show_tiled_fitted_segments(shell_segments, fits)
title('All segmented shell images with fit overlays')
