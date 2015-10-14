x_centre  = 5;
y_centre  = 7;
radius    = 3;
psf_sigma = 1;
height    = 1;
thickness = 2; % EJR's code uses half-thickness
X(:, 1)   = 0:10;
X(:, 2)   = 0:10;

% For floating point comparison
tolerance = 1E-7;

%% Test that thin sphere image has right intensities
sphere_thin = fsa.image_sphere_thin(x_centre, y_centre, radius, psf_sigma, height, X);
assert(all(abs([1.7782397997146e-08;1.95543222989977e-05;0.00311874507453789;0.0756635854607892;0.312091279411532;0.303263466529731;0.201062512910099;0.303263466529731;0.312091279411532;0.0756635854607892;0.00311874507453789] - sphere_thin) <= tolerance), 'Thin sphere output incorrect.');


%% Test that thick sphere image has right intensities
sphere_thick = fsa.image_sphere_thick(x_centre, y_centre, radius, psf_sigma, height, thickness, X);
assert(all(abs([4.62568635601214e-06;0.00171481448673769;0.103947583982774;1.18364166944416;3.26728012668752;3.25706406712521;2.51312391692102;3.25706406712521;3.26728012668752;1.18364166944416;0.103947583982774] - sphere_thick) <= tolerance), 'Thick sphere output incorrect.');
