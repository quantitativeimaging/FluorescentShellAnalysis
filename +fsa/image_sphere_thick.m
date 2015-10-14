function I = image_sphere_thick(x_centre, y_centre, radius, psf_sigma, height, thickness, X)
% IMAGE_SPHERE_THICK - Return radial intensities of image of a thick spherical shell
%
% Input:
% 	x_centre  - x coordinate of shell centre.
% 	y_centre  - y coordinate of shell centre.
% 	radius    - Shell radius.
% 	psf_sigma - Sigma of Gaussian point spread function.
% 	height    - Height of image intensity.
% 	thickness - Thickness of shell.
% 	X         - Array of (x, y) coordinates.
%
% Output:
% 	I - Vector of pixel values at each X (x, y) position.

% Radial position
r = sqrt((X(:, 1) - x_centre).^2 + (X(:, 2) - y_centre).^2);

psf_variance = psf_sigma^2;
radius_inner = radius - thickness / 2;
radius_outer = radius + thickness / 2;

I = (height * psf_sigma ./ r) .* ( ...
	  (sqrt(2 * pi) .* r .* erf((radius_outer - r) / (sqrt(2) * psf_sigma))) ...
	+ (sqrt(2 * pi) .* r .* erf((radius_outer + r) / (sqrt(2) * psf_sigma))) ...
	+ (2 * psf_sigma .* exp(-(radius_outer + r).^2 / (2 * psf_variance))) ...
	- (2 * psf_sigma .* exp(-(radius_outer - r).^2 / (2 * psf_variance))) ...
	- (sqrt(2 * pi) .* r .* erf((radius_inner - r) / (sqrt(2) * psf_sigma))) ...
	- (sqrt(2 * pi) .* r .* erf((radius_inner + r) / (sqrt(2) * psf_sigma))) ...
	- (2 * psf_sigma .* exp(-(radius_inner + r).^2 / (2 * psf_variance))) ...
	+ (2 * psf_sigma .* exp(-(radius_inner - r).^2 / (2 * psf_variance))) ...
	);

end
