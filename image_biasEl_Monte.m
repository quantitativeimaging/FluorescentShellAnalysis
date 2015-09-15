function I = image_unifEl_Monte(b0, X)% Monte Carlo uniform shell image

flagShowModelImages = 0;
flagSaveModelIms    = 0; % Need to pre-set  count = 1  in base workspace;

xcen   = b0(1); % X-coordinate of image centre
ycen   = b0(2); % Y-coordinate of image centre
crcrad = b0(3); % semi-minor axis of prolate ellipsoid of revolution
var    = b0(4); % sigma squared (optical PSF variance on an axis)
height = b0(5); % height
ellip  = b0(6)+1; % Ellipse aspect ratio - 1 (c/a - 1) or shape factor -1
psi    = b0(7); % Azimuthal orientation
Q      = b0(8); % Equatoriality. A bias parameter for protein density

rr = sqrt( (X(:,1)-xcen).^2 + (X(:,2)-ycen).^2 ); % radial position
a  = crcrad;
twoSS = 2*abs(var);

nPoints = 3000;  % Number of fluorophores to simulate (before rejection)
if(length(b0) > 8)
   nPoints = b0(9); 
end

phi      = 2*pi*rand(nPoints,1);
costheta = 2*rand(nPoints,1) -1;
sintheta = sqrt(1-costheta.^2);

semiMajor = a * ellip; % semi-minor * aspect ratio
semiMinor = a;

surfXp = semiMajor.*sintheta.*cos(phi); % Let Xp be long-axis co-ord
surfYp = semiMinor.*sintheta.*sin(phi); % Let Yp be short-axis co-ord
surfZ = semiMinor*costheta; % Zp is short-axis co-ord, long axis in Zp=0.

% Now discard some points randomly to get uniform sampling on surface
g = sqrt( (surfYp./semiMinor).^2 + (surfZ./semiMinor).^2 + ...
          ((surfXp.*semiMinor).^2)./(semiMajor^4) ); % accept ratio. Ref Chen.
% Bias g by the equatoriality. 
sinTsq = (surfYp.^2 + surfZ.^2)./sqrt(surfXp.^2 + surfYp.^2 + surfZ.^2);
sinT   = sqrt(sinTsq);
g = g.*(1 + Q.*sinT);
      
pAccept = g./max(g(:));
chance  = rand(length(pAccept),1); % Or =0; to accept all
surfXp  = surfXp(chance < pAccept);
surfYp  = surfYp(chance < pAccept);
surfZ   = surfZ(chance < pAccept);

surfX = surfXp.*cos(psi) + surfYp.*sin(psi)  + xcen; % Rotate about Z-axis
surfY = -surfXp.*sin(psi) + surfYp.*cos(psi) + ycen;

nPointsAccepted = length(surfX);

I = zeros(size(X,1), 1);
for lp = 1:nPointsAccepted
    dispsSq = ( (surfX(lp) - X(:,1) ).^2 + ...
                (surfY(lp) - X(:,2) ).^2 + ....
                (surfZ(lp)).^2 );
              
    ints  = exp(-(dispsSq)/twoSS);
    I = I + ints;
end

I = I * height / max( I(:) );


if(flagShowModelImages)
  im = zeros(41);
  for lp = 1:length(I)
    im(X(lp,1),X(lp,2)) = I(lp);
  end
  figure(9)
  imagesc(im') % Note transpose for image orientation!
  % figure(10)
  % scatter3(surfX, surfY, surfZ)
end
  
if(flagSaveModelIms)
   count = evalin('base', 'count');
   imwrite(double(im')./max(im(:)), ['C:\Users\user\Documents\Projects\2014_Spores\MATLAB\out\out',int2str(count), '.tif'])
   count = count+1;
   assignin('base','count', count);
end

I(isnan(I)) = 0;    % Prevent NaN returns
I(isinf(I)) = height/2; % Prevent inf returns (Better to bound ss?)
end