x_centre = 5;
y_centre = 7;
radius = 3;
psf_sigma = 1;
height = 1;
X(:, 1) = 0:10;
X(:, 2) = 0:10;

tolerance = 1E-7;

sphere_thin = fsa.image_sphere_thin(x_centre, y_centre, radius, psf_sigma, height, X);

%% Test that thin sphere image has right intensities
assert(all(abs([1.7782397997146e-08;1.95543222989977e-05;0.00311874507453789;0.0756635854607892;0.312091279411532;0.303263466529731;0.201062512910099;0.303263466529731;0.312091279411532;0.0756635854607892;0.00311874507453789] - sphere_thin) <= tolerance), 'Thin sphere output incorrect.');
