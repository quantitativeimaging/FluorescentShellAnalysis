function beta = fitSphereParams(X, listI, b0)
% Estimate best-fit parameters of a spherical shell to fit image data
% This method tries to explore and shrink a region of confidence 
% based around an initial guess of parameters for the mixed model. 
% It is an heuristic method 
% The results should be checked against spore images by the user

% xcen   = b0(1); % X-coordinate of image centre
% ycen   = b0(2); % Y-coordinate of image centre
% crcrad = b0(3); % radius of spherical shell
% var    = b0(4); % sigma squared (optical PSF variance on an axis)
% height = b0(5); % height
% ellip = beta(6); % Ellipticity, which is not used in this model

% I     = image_sphere_Monte(b0, X);
% sumSq = sum((I - listI).^2);       % Quantifies misfit at initial guess

radX   = 0.4;       % parameter adjustments to consider
radY   = 0.4;
radR   = 0.2*b0(3);
radVar = 0.5*b0(4);
radHt  = 0.1*b0(5);

I2     = image_sphere_Monte(b0 + [0,0, 0,+0.2,0,0], X);
sumSq2 = sum((I2 - listI).^2);

numberIts = 25;
shift     = 0.95;     % The range 0.9 to 0.95 seems reasonable
listParams= zeros(numberIts*2,6);

for lpIts = 1:numberIts
    I     = image_sphere_Monte(b0, X);
    sumSq = sum((I - listI).^2);       % Quantifies misfit at initial guess

    % Check for sphere radius improvement
    IradHi = image_sphere_Monte(b0 + [0,0,radR,0,0,0], X);
    IradLo = image_sphere_Monte(b0 - [0,0,radR,0,0,0], X);
    ssRadH = sum((IradHi - listI).^2);
    ssRadL = sum((IradLo - listI).^2);
    if(ssRadH < sumSq && ssRadH < ssRadL)
        b0(3) = b0(3) + radR/2;
    elseif(ssRadL < sumSq && ssRadL < ssRadH)
        b0(3) = b0(3) - radR/2;
    end
    radR = shift*radR;
    
    % Check for blur radius (point spread function) improvement 
    I     = image_sphere_Monte(b0, X);
    sumSq = sum((I - listI).^2);
    IvarHi = image_sphere_Monte(b0 + [0,0,0,radVar,0,0], X);
    IvarLo = image_sphere_Monte(b0 - [0,0,0,radVar,0,0], X);
    ssVarH = sum((IvarHi - listI).^2);
    ssVarL = sum((IvarLo - listI).^2);
    if(ssVarH < sumSq && ssVarH < ssVarL)
        b0(4) = b0(4) + radVar/2;
    elseif(ssVarL < sumSq && ssVarL < ssVarH)
        b0(4) = b0(4) - radVar/2;
    end
    radVar = shift*radVar;
    
    % Check for brightness (signal height) improvement 
    I     = image_sphere_Monte(b0, X);
    sumSq = sum((I - listI).^2);
    IhtHi = image_sphere_Monte(b0 + [0,0,0,0,radHt,0], X);
    IhtLo = image_sphere_Monte(b0 - [0,0,0,0,radHt,0], X);
    ssHtH = sum((IhtHi - listI).^2);
    ssHtL = sum((IhtLo - listI).^2);
    if(ssHtH < sumSq && ssHtH < ssHtL)
        b0(5) = b0(5) + radHt/2;
    elseif(ssHtL < sumSq && ssHtL < ssHtH)
         b0(5) = b0(5) - radHt/2;
    end
    radHt = radHt*shift;
    
    % Check for centre co-ordinate improvement (X-direction)
    I     = image_sphere_Monte(b0, X);
    sumSq = sum((I - listI).^2);
    IxcHi = image_sphere_Monte(b0 + [radX,0,0,0,0,0], X);
    IxcLo = image_sphere_Monte(b0 - [radX,0,0,0,0,0], X);
    ssXcH = sum((IxcHi - listI).^2);
    ssXcL = sum((IxcLo - listI).^2);
    if(ssXcH < sumSq && ssXcH < ssXcL)
        b0(1) = b0(1) + radX/2;
    elseif(ssXcL < sumSq && ssXcL < ssXcH)
         b0(1) = b0(1) - radX/2;
    end
    radX = radX*shift;

    % Check for centre co-ordinate improvement (Y-direction)
    I     = image_sphere_Monte(b0, X);
    sumSq = sum((I - listI).^2);
    IycHi = image_sphere_Monte(b0 + [0,radY,0,0,0,0], X);
    IycLo = image_sphere_Monte(b0 - [0,radY,0,0,0,0], X);
    
    ssYcH = sum((IycHi - listI).^2);
    ssYcL = sum((IycLo - listI).^2);
    if(ssYcH < sumSq && ssYcH < ssYcL)
        b0(2) = b0(2) + radY/2;
    elseif(ssYcL < sumSq && ssYcL < ssYcH)
         b0(2) = b0(2) - radX/2;
    end
    radY = radY*shift;
    
    listParams(lpIts,:) = b0;
    
    figure(7)
    rr = sqrt((X(:,1)-b0(1)).^2 + (X(:,2)-b0(2)).^2);
    plot(rr, listI)
    hold on
      plot(rr,I,'g')
    hold off
    legend('Data','Fit');
    
end

% Further iterature to refine radius.
for lpIts = (numberIts+1): (2*numberIts)
    
    I     = image_sphere_Monte(b0, X);
    sumSq = sum((I - listI).^2);       % Quantifies misfit at initial guess

    % Check for sphere radius improvement
    IradHi = image_sphere_Monte(b0 + [0,0,radR,0,0,0], X);
    IradLo = image_sphere_Monte(b0 - [0,0,radR,0,0,0], X);
    ssRadH = sum((IradHi - listI).^2);
    ssRadL = sum((IradLo - listI).^2);
    if(ssRadH < sumSq && ssRadH < ssRadL)
        b0(3) = b0(3) + radR/2;
    elseif(ssRadL < sumSq && ssRadL < ssRadH)
        b0(3) = b0(3) - radR/2;
    end
    radR = shift*radR;
    
    listParams(lpIts,:) = b0;
    
    figure(7)
    rr = sqrt((X(:,1)-b0(1)).^2 + (X(:,2)-b0(2)).^2);
    plot(rr, listI)
    hold on
      plot(rr,I,'g')
    hold off
    legend('Data','Fit');
end
%

% Save an estimate of near-final fit quality to the base workspace
relSumSq = sumSq / sum(listI.^2);
assignin('base', 'relSumSq', relSumSq);

figure(8)
plot(listParams(:,3), 'b');
hold on
 plot(listParams(:,4), 'g');
 %plot(listParams(:,5)), 'r';
 plot(listParams(:,1), 'r');
 plot(listParams(:,2), 'k');
hold off
legend('radius', 'var', 'xCen', 'yCen' )
xlabel('fit iterations')

beta = b0;
end