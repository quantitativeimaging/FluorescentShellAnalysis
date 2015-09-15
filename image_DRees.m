function I = image_DRees(beta, X)%UNTITLED2 Summary of this function goes here
%   Returns I = intensities in image of a thin spherical fluorescent shell

flagCompareDataModel = 0;

xcen   = beta(1); % X-coordinate of image centre
ycen   = beta(2); % Y-coordinate of image centre
crcrad = beta(3); % radius of spherical shell
var    = beta(4); % sigma squared (optical PSF variance on an axis)
height = beta(5); % height
%ellip = beta(6); % Ellipticity, which is not used in this model

rr = sqrt( (X(:,1)-xcen).^2 + (X(:,2)-ycen).^2 ); % radial position
a = crcrad;
ss = var;        

% Equation for integral intensity from sphere:
I = height * (exp(-(rr-a).^2/(2*ss)) - exp(-(rr+a).^2/(2*ss)) )./rr;

% In case rr = 0 exactly, a singular point, use the limiting case of I(r):
I(rr==0) = (2*a*height/(ss))*exp(-(a^2)/(2*ss));

I(isnan(I)) = 0;         % Prevent NaN returns (Should be none now)
I(isinf(I)) = height*10; % Prevent inf returns (Should be none)

if(flagCompareDataModel)
   figure(11)
   plot(rr*74, evalin('base', 'listI'), 'b');
   hold on
     plot(rr*74, I, 'r');
   hold off
   
  % set(gca,'yTick',[0 1000])
  %  AX2 = legend('Data', 'Model', 'FontSize', 14);
  % LEG = findobj(AX2,'type','text');
  % set(LEG,'FontSize',10)
  % legend boxoff
  xlim([0 1500])
  ylim([0 1200])
  set(gca,'FontSize',16,'fontweight','normal');
  xlabel('r / nm')
  ylabel(' pixel value', 'FontSize', 16);
  set(gcf,'color','w')
  % set(10,'Position',[100,100,330,280]); % 720 px wide, 600 high
end

end