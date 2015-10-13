function I = image_sphere_thin(x_centre, y_centre, radius, psf_sigma, height, X)
% IMAGE_SPHERE_THIN - Return radial intensities of image of a thin spherical shell
%
% Input:
% 	x_centre  - x coordinate of shell centre.
% 	y_centre  - y coordinate of shell centre.
% 	radius    - Shell radius.
% 	psf_sigma - Sigma of Gaussian point spread function.
% 	height    - Height of image intensity.
% 	X         - Array of (x, y) coordinates.
%
% Output:
% 	I - Vector of radial image intensities.

% Radial position
r = sqrt((X(:, 1) - x_centre).^2 + (X(:, 2) - y_centre).^2);
psf_variance = psf_sigma^2;

I = height * (exp(-(r - radius).^2 / (2 * psf_sigma)) - exp(-(r + radius).^2 / (2 * psf_sigma))) ./ r;

% For r = 0, a singular point, use the limiting case
I(r == 0) = (2 * radius * height / psf_variance) * exp(-(radius^2) / (2 * psf_variance));

end
