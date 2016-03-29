actual_image = double(imread('one_ellipsoid.tif'));

actual_image = actual_image(6:32, 5:31);

background = median(actual_image(actual_image < mean(actual_image(:)))); % 3/15. Good. 
actual_image = actual_image - background;

bw_image = actual_image;
threshold = 35 - background;
bw_image(actual_image > threshold) = 1;
bw_image(actual_image <= threshold) = 0;

% bw_image = actual_image > mean(actual_image(:));

stats = regionprops(bw_image, 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation');

x_shift = 0;
y_shift = 0;
orientation = -1;
% orientation = stats.Orientation * pi / 180;
% semiminor_axis = 4;
% semiminor_axis = stats.MinorAxisLength / 2;
semiminor_axis = 6;
psf_variance = 6;
height = max(actual_image(:));
% eccentricity = 1.3;
% eccentricity = stats.MajorAxisLength / stats.MinorAxisLength - 1;
eccentricity = 0.2;
equatoriality = 0;

ori_image = fsa.image_ellipsoid_biased(x_shift, y_shift, orientation, semiminor_axis, psf_variance, height, eccentricity, equatoriality, actual_image);

[x_shift, y_shift, orientation, semiminor_axis, psf_variance, height, eccentricity, equatoriality] = fsa.fit_ellipsoid(x_shift, y_shift, orientation, semiminor_axis, psf_variance, height, eccentricity, equatoriality, actual_image)

fit_image = fsa.image_ellipsoid_biased(x_shift, y_shift, orientation, semiminor_axis, psf_variance, height, eccentricity, equatoriality, actual_image);


figure
subplot(2, 2, 1)
imshow(actual_image, [])
title('Actual')
subplot(2, 2, 2)
imshow(bw_image, [])
title('Thresholded')
subplot(2, 2, 3)
imshow(ori_image, [])
title('Initial')
subplot(2, 2, 4)
imshow(fit_image, [])
title('Final')