function I = image_cylinder_thin(x_centre, y_centre, cylinder_angle, cylinder_length, radius, psf_sigma, height, imagemat)

image_width = size(imagemat, 2);
image_height = size(imagemat, 1);

image_centre_x = image_width / 2;
image_centre_y = image_height / 2;

for x = 1:image_width
	for y = 1:image_height
		x_pos = x - image_centre_x;
		y_pos = y - image_centre_y;

		% Rotate by angle
		x_pos = x_pos * cos(cylinder_angle) + y_pos * -sin(cylinder_angle);
		y_pos = x_pos * sin(cylinder_angle) + y_pos * cos(cylinder_angle);

		% Shift centre
		x_pos = x_pos - x_centre * cos(cylinder_angle) - y_centre * sin(cylinder_angle);
		y_pos = y_pos - x_centre * sin(cylinder_angle) + y_centre * cos(cylinder_angle);

		imagemat(y, x) = fsa.cross_section_cylinder_thin(0, radius, psf_sigma, height, x_pos);
	end
end

I = imagemat;

end
