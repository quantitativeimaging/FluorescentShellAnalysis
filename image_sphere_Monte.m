function I = image_sphere_Monte(b0, X)% Monte Carlo spherical shell image

xcen   = b0(1); % X-coordinate of image centre
ycen   = b0(2); % Y-coordinate of image centre
crcrad = b0(3); % radius of spherical shell
var    = b0(4); % sigma squared (optical PSF variance on an axis)
height = b0(5); % height
%ellip = beta(6); % Ellipticity, which is not used in this model

rr = sqrt( (X(:,1)-xcen).^2 + (X(:,2)-ycen).^2 ); % radial position
a  = crcrad;
twoSS = 2*var;

nPoints = 2500;  % Number of fluorophores to simulate

phi      = 2*pi*rand(nPoints,1);
costheta = 2*rand(nPoints,1) -1;
sintheta = sqrt(1 - costheta.^2);

surfX = a*sintheta.*cos(phi) +xcen;   % Project to surface. Isotropic. 
surfY = a*sintheta.*sin(phi) +ycen;
surfZ = a*costheta;

I = zeros(size(X,1), 1);
for lp = 1:nPoints
    dispsSq = ( (surfX(lp) - X(:,1) ).^2 + ...
                (surfY(lp) - X(:,2) ).^2 + ....
                (surfZ(lp)).^2 );
              
    ints  = exp(-(dispsSq)/twoSS);
    I = I + ints;
end

I = I * height / max( I(:) );

I(isnan(I)) = 0;    % Prevent NaN returns
I(isinf(I)) = 1E10; % Prevent inf returns (Better to bound ss?)
end