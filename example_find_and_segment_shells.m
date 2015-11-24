image_data = imread('testData_megaterium_spheres.tif');

% Parameters for shell finding
radius_lower = 5;
radius_upper = 15;
segment_half_size = 40;
edge_border = 0;

% Find and display circular shells
[centres, radii, metric] = fsa.find_circular_shells(image_data, radius_lower, radius_upper, segment_half_size, edge_border, true);
title('Identified circular shells')

% Tile segmented shells in figure
shell_segments = fsa.segment_shells(image_data, centres, segment_half_size);
tiled_segments = fsa.tile_segments(shell_segments);
figure
imshow(tiled_segments, [])
