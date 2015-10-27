function I = cross_section_cylinder_thin(x_centre, radius, psf_sigma, height, x)

r = x - x_centre;

I = height * sqrt(2) * exp(-(radius^2 + r.^2) / (2 * psf_sigma^2)) * pi^(3/2) * psf_sigma .* besseli(0, radius .* r / psf_sigma^2);

end
