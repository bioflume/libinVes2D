classdef curve 
% This class implements that basic calculus on the curve.
% The basic data structure is a matrix X in which the columns 
% represent periodic C^{\infty} closed curves with N points, 
% X(1:n,j) is the x-coordinate of the j_th curve and X(n+1:N,j) 
% is the y-coordinate of the j_th curve; here n=N/2
% X coordinates do not have to be periodic, but the curvature,
% normals, etc that they compute will be garbage.  This allows
% us to store target points for tracers or the pressure using
% this class and then using near-singular integration is easy
% to implement

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y]=getXY(o,X)
% [x,y] = getXY(X) get the [x,y] component of curves X
N = size(X,1)/2;
x = X(1:N,:);
y = X(N+1:end,:);

end % getXY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V = setXY(o,x,y)
% V = setXY(x,y) set the [x,y] component of vector V on the curve
N = size(x,1);
V=zeros(2*N,size(x,2));
V(1:N,:) = x;
V(N+1:end,:) = y;

end % setXY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dx,Dy]=getDXY(o,X)
% [Dx,Dy]=getDXY(X), compute the derivatives of each component of X 
% these are the derivatives with respect the parameterization 
% not arclength
[x,y] = o.getXY(X);
N = size(x,1);
nv = size(x,2);
IK = fft1.modes(N,nv);
Dx = fft1.diffFT(x,IK);
Dy = fft1.diffFT(y,IK);

end % getDXY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [jacobian,tangent,curvature] = diffProp(o,X)
% [jacobian,tangent,curvature] = diffProp(X) returns differential
% properties of the curve each column of the matrix X. Each column of 
% X should be a closed curve defined in plane. The tangent is the 
% normalized tangent.
%
% EXAMPLE:
%    n = 128; nv = 3;
%    X = boundary(n,'nv',nv,'curly');
%    c = curve;
%    [k t s] = c.diffProp(X);

n = size(X,1);
nv = size(X,2);
N = n/2; 
IK = fft1.modes(N,nv);

% get the x y components
[Dx,Dy] = o.getDXY(X);

jacobian = sqrt( Dx.^2 + Dy.^2 ); 

if nargout>1  % if user requires tangent
  tangent = o.setXY( Dx./jacobian, Dy./jacobian);
end

if nargout>2  % if user requires curvature
  DDx = curve.arcDeriv(Dx,1,ones(N,nv),IK);
  DDy = curve.arcDeriv(Dy,1,ones(N,nv),IK);
  curvature = (Dx.*DDy - Dy.*DDx)./(jacobian.^3);
end
% curvature of the curve

end % diffProp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [reducedArea,area,length] = geomProp(o,X)
% [reducedArea area length] = geomProp(X) calculate the length, area 
% and the reduced volume of domains inclose by columns of X. 
% Reduced volume is defined as 4*pi*A/L. 
% EXAMPLE:
%   X = boundary(64,'nv',3,'curly');
%   c = curve(X);
%   [rv A L] = c.geomProp(X);

[x,y] = o.getXY(X);
N = size(x,1);
[Dx,Dy] = o.getDXY(X);

length = sum( sqrt(Dx.^2 + Dy.^2) )*2*pi/N;
area = sum( x.*Dy - y.*Dx)*pi/N;

reducedArea = 4*pi*area./length.^2;

end % geomProp



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,nv] = initConfig(o,N,varargin)       
% [X,nv] = initConfig(n,varargin) returns N coordinates of boundary
% points.
% X = BOUNDARY(N,OPTIONS) can be used to call for different
% configuration.  Available OPTIONS are
%
%   'nv'          - (followed by) number of vesicles (in 
%                   lineajr configuration),
%   'angle'       - inclination angle of the vesicle(s) form the vertical
%                   position, 
%   'curly'       - returns a somehow wiggly vesicle.
%   'couette'     - the vesicles inside the default geometry of 
%                   couette flow (confined). 
%   'scale'       - multiplies the size of the output boundary
%   'choke'       - returns a choked domain.  Usually a solid boundary
%
%   'couette'     - returns a domain for a couette apparatus.
%   'figureEight' - returns a domain that resembles a pinched ellipse
%   'tube'        - returns a domain that is an ellongated ellipse

% EXAMPLE: X = boundary(64,'nv',3,'theta',pi/6);
%

options = varargin;
X = [];

if(any(strcmp(options,'nv')))
  nv = options{find(strcmp(options,'nv'))+1};
else
  nv = 1;
end
% number of vesicles

if(any(strcmp(options,'angle')))
  theta = options{find(strcmp(options,'angle'))+1};
else
  theta = zeros(nv,1);
end
% rotate the vesicle by a certain angle

if(any(strcmp(options,'center')))
  cen = options{find(strcmp(options,'center'))+1};
else
  cen = [(0:nv-1);zeros(1,nv)];
end
% pick the center of the vesicles

if(any(strcmp(options,'reducedArea')))
  ra = options{find(strcmp(options,'reducedArea'))+1};
else 
  ra = [];
end
% desired reduced area

if(any(strcmp(options,'scale')))
  scale = options{find(strcmp(options,'scale'))+1};
else
  scale = 1;
end
% scale the size of the vesicle (reduced area is invariant under
% scale)


t = (0:N-1)'*2*pi/N;
% Discretization in parameter space
  
if any(strcmp(options,'curly'))
  a = 1; b = 3*a; c = 0.85; 
  r = 0.5*sqrt( (a*cos(t-theta)).^2 + (b*sin(t-theta)).^2) + ...
      .07*cos(12*(t-theta));
  % radius of curly vesicle
  X0 = scale*[c*r.*cos(t); r.*sin(t)];

elseif any(strcmp(options,'crazy'))
  mollify = exp(1./((t/pi-1).^2-1));
  mollify(1) = 0;
  mollify(end) = 0;
  % mollifier supported on [0,2*pi] and centered at pi

  index = ceil(N/8);
  mollify1 = [mollify(index+1:end);mollify(1:index)];
  % shift the center of the mollifier

  r1 = .8*cos(16*t) + .16*sin(20*t);
  r1 = r1.*mollify1;
  % radius with additional high frequencies 'wiggles'

  index = ceil(N/2);
  mollify2 = [mollify(index+1:end);mollify(1:index)];
  % shift the center of the mollifier
  r2 = 2.4*sin(7*t) + .64*cos(10*t);
  r2 = r2.*mollify2;
  % radius with additional high frequencies 'wiggles'

  X = [8*cos(t) + r1.*cos(t) + r2.*cos(t);...
      sin(t) + r1.*sin(t) + r2.*sin(t)];
  % A somewhat complicated geometry

elseif any(strcmp(options,'star'))
  folds = options{2};
  radius = 1 + 0.98*cos(folds*t);
  X = [radius.*cos(t);radius.*sin(t)];
  % a star that comes very close to intersecting itself
  % at the origin

elseif any(strcmp(options,'openStar'))
  if any(strcmp(options,'amplitude'))
    amp = options{find(strcmp(options,'amplitude'))+1};
  else
    amp = 5e-2;
  end
  folds = options{2};
  radius = 1 + 0.5*cos(folds*t) + amp*cos(30*t);
  X = [radius.*cos(t);radius.*sin(t)];
  % a star that comes very close to intersecting itself
  % at the origin

elseif any(strcmp(options,'wigglyStar'))
%  radius = 1 + .8*cos(2*t);
%  X = [radius.*cos(t);radius.*sin(t)] + [0.01*cos(80*t).*cos(t);0.01*sin(80*t).*sin(t)];
  % a star with only two folds, but a high frequency
  % added on top
  radius = 1 + 0.5*cos(3*t) + 0.05*cos(30*t);
  X = [radius.*cos(t);radius.*sin(t)];

elseif any(strcmp(options,'spikes'))
  % This shape is not implemented very well.  The curve
  % is too crazy
  N = numel(t);
  uprate = 5;
  t1 = t;
  t2 = (uprate*N-1:-1:0)'*2*pi/(uprate*N);
  % want to upsample the top part of the vesicle since it
  % has much sharper features

  x = [t1;t2];
  x = x(1:(uprate+1):end);
  y = [0.3*cos(10*t1)+0.7;ones(size(t2))+0.1*(cos(2*t2)+1.1)];
  y = y(1:(uprate+1):end);

  width = 0.2;
  cen = [0.5 2 3 3.8 5];
  height = [5 9 13 9 5];
  for k = 1:numel(cen)
    s = find(x > cen(k)-width & x < cen(k)+width & y > 1.01);
    y(s) = y(s) + height(k)*cos(pi/2/width*(x(s)-cen(k)));
  end

  z = x+1i*y;
  z = z - 1/2*(max(x) + min(x)) - 1i/2*(max(y) + min(y));
  z1 = z;
  zh = fftshift(fft(z));
  modes = (-N/2+1:N/2)';
  zh(abs(modes) > 100) = 0;
  z = ifft(ifftshift(zh));
  x = real(z); y = imag(z);
  X = [x;y];


elseif any(strcmp(options,'choke'))
  a = 10; b = 3; c = 0.6; order = 8;
  % parameters for the boundary
  Nsides = ceil(0.5*b/(2*a+2*b)*N);
  Ntop = (N-4*Nsides)/2;
  t1 = linspace(0,0.2*pi,Nsides+1); t1 = t1(1:end-1)';
  t2 = linspace(0.2*pi,pi-0.2*pi,Ntop+1); t2 = t2(1:end-1)';
  t3 = linspace(pi-0.2*pi,pi+0.2*pi,2*Nsides+1); t3 = t3(1:end-1)';
  t4 = linspace(pi+0.2*pi,2*pi-0.2*pi,Ntop+1); t4 = t4(1:end-1)';
  t5 = linspace(2*pi-0.2*pi,2*pi,Nsides+1); t5 = t5(1:end-1)';
  t = [t1;t2;t3;t4;t5];
  % Parameterize t so that geometry is closer to equi-spaced in arclength
%  t = (0:N-1)'*2*pi/N;
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  ind = abs(x) < pi;
  y(ind) = y(ind).*(1-c*cos(x(ind)))/(1+c);
  X0 = [x;y];
  % choked domain.  a and b control the length and height.  c
  % controls the width of the gap, and order controls the
  % regularity

elseif any(strcmp(options,'choke2'))
  a = 40; b = 3; c = 0.6; order = 8;
  % parameters for the boundary
  Nsides = ceil(0.5*b/(2*a+2*b)*N);
  Ntop = (N-4*Nsides)/2;
  t1 = linspace(0,0.2*pi,Nsides+1); t1 = t1(1:end-1)';
  t2 = linspace(0.2*pi,pi-0.2*pi,Ntop+1); t2 = t2(1:end-1)';
  t3 = linspace(pi-0.2*pi,pi+0.2*pi,2*Nsides+1); t3 = t3(1:end-1)';
  t4 = linspace(pi+0.2*pi,2*pi-0.2*pi,Ntop+1); t4 = t4(1:end-1)';
  t5 = linspace(2*pi-0.2*pi,2*pi,Nsides+1); t5 = t5(1:end-1)';
  t = [t1;t2;t3;t4;t5];
  % Parameterize t so that geometry is closer to equi-spaced in arclength
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  ind = abs(x) < pi;
  y(ind) = y(ind).*(1-c*cos(x(ind)))/(1+c);
  X0 = [x;y];
  % choked domain.  a and b control the length and height.  c
  % controls the width of the gap, and order controls the
  % regularity

elseif any(strcmp(options,'doublechoke'))
  a = 10; b = 3; c = 0.6; order = 8;
  shift = pi/2 + 0.1;
  % parameters for the boundary
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  ind = abs(x-shift) < pi/2;
  y(ind) = y(ind).*(1-c*cos(2*(x(ind)-shift)))/(1+c);
  ind = abs(x+shift) < pi/2;
  y(ind) = y(ind).*(1-c*cos(2*(x(ind)+shift)))/(1+c);
  X0 = [x;y];
  % choked domain.  a and b control the length and height.  c
  % controls the width of the gap, and order controls the
  % regularity.  shift controls the distance between the two
  % regions where the domain is restricted

elseif any(strcmp(options,'couette'))
  x = [20*cos(t)+cen(1,1) 10*cos(-t)+cen(1,2)];
  y = [20*sin(t)+cen(2,1) 10*sin(-t)+cen(2,2)];
  X = o.setXY(x,y);
  % annular domain

elseif any(strcmp(options,'doubleCouette'))
  x = [20*cos(t)+cen(1,1) 5*cos(-t)+cen(1,2) 5*cos(-t)+cen(1,3)];
  y = [20*sin(t)+cen(2,1) 5*sin(-t)+cen(2,2) 5*sin(-t)+cen(2,3)];
  X = o.setXY(x,y);
  % annular domain of genus 2

elseif any(strcmp(options,'quadCouette'))
  x = [20*cos(t)+cen(1,1) 5*cos(-t)+cen(1,2) 5*cos(-t)+cen(1,3) ...
                          5*cos(-t)+cen(1,4) 5*cos(-t)+cen(1,5)];
  y = [20*sin(t)+cen(2,1) 5*sin(-t)+cen(2,2) 5*sin(-t)+cen(2,3) ...
                          5*sin(-t)+cen(2,4) 5*sin(-t)+cen(2,5)];
  X = o.setXY(x,y);
  % annular domain of genus 4

elseif any(strcmp(options,'doubleFlower'))
  r = 17 + 2*cos(7*t);
  x = [r.*cos(t)+cen(1,1) 5*cos(-t)+cen(1,2) 5*cos(-t)+cen(1,3)];
  y = [r.*sin(t)+cen(2,1) 5*sin(-t)+cen(2,2) 5*sin(-t)+cen(2,3)];
  X = o.setXY(x,y);
  % annular domain of genus 2

elseif any(strcmp(options,'cylinder'))
  x = [20*scale*cos(t)+cen(1,1)];
  y = [20*scale*sin(t)+cen(2,1)];
  X = o.setXY(x,y);
  % single cylinder

elseif any(strcmp(options,'figureEight'))
  r = 1 + .7*cos(2*t);
  x = r.*cos(t); x = 2*x/max(x);
  y = r.*sin(t); y = y/max(y); 
  X = o.setXY(x,y);
  % pinched elliptical cylinder

elseif any(strcmp(options,'tube'))
  a = 10; b = 3; order = 8;
  % parameters for the boundary
  Nsides = ceil(0.5*b/(2*a+2*b)*N);
  Ntop = (N-4*Nsides)/2;
  t1 = linspace(0,0.2*pi,Nsides+1); t1 = t1(1:end-1)';
  t2 = linspace(0.2*pi,pi-0.2*pi,Ntop+1); t2 = t2(1:end-1)';
  t3 = linspace(pi-0.2*pi,pi+0.2*pi,2*Nsides+1); t3 = t3(1:end-1)';
  t4 = linspace(pi+0.2*pi,2*pi-0.2*pi,Ntop+1); t4 = t4(1:end-1)';
  t5 = linspace(2*pi-0.2*pi,2*pi,Nsides+1); t5 = t5(1:end-1)';
  t = [t1;t2;t3;t4;t5];
  % Parameterize t so that geometry is closer to 
  % equispaced in arclength
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  X0 = [x;y];
  % rounded off cylinder.  a and b control the length and height 
  % and order controls the regularity

else
  if ~isempty(ra)
    X0 = o.ellipse(N,ra);
    % build a vesicle of reduced area ra with N points
  else
    a = 1; b = 4.8*a; c = 1; 
%    r = 0.5*sqrt( (a*cos(t-theta)).^2 + (b*sin(t-theta)).^2);
    r = 0.5*sqrt( (a*cos(t)).^2 + (b*sin(t)).^2);
    X0 = [c*r.*cos(t); r.*sin(t)];
    % this shape has reduced area very close to 0.65
  end
  oc = curve;
  [ra,area,length] = oc.geomProp(X0);
  scale = 2*scale/sqrt(area/pi);
  X0 = scale*X0;
  % build reference elliptical vesicle

  if any(strcmp(options,'volFrac'))
    ind = find(strcmp(options,'volFrac'));
    volFrac = options{ind+1};
    Xwalls = options{ind+2};
    fmm = options{ind+3};
    X = o.fillDomain(X0,volFrac,Xwalls,fmm);
    nv = size(X,2);
  end


  if size(cen,2) ~= nv
    b = 3*max(X0(1:N));
    cen = [b*(0:nv-1);zeros(1,nv)];
  end
  % centers if there are multiple vesicles.
end 
% end of building reference vesicles.  Only need to rotate
% and shift them as desired

if isempty(X)
  % if X has not been defined, we only have X0 which is 
  % the reference vesicle which not needs to be
  % rotated and translated
  X = zeros(2*N,nv);
  for k=1:nv
    X(1:N,k) = cos(theta(k)) * X0(1:N) - ...
      sin(theta(k)) * X0(N+1:2*N) + cen(1,k);
    X(N+1:2*N,k) = sin(theta(k)) * X0(1:N) +  ...
      cos(theta(k)) * X0(N+1:2*N) + cen(2,k);
  end
  % Rotate vesicles as requested 
end


end % initConfig



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X0 = ellipse(o,N,ra)
% X0 = o.ellipse(N,ra) finds the ellipse (a*cos(theta),sin(theta)) so
% that the reduced area is ra.  Uses N points.  Parameter a is found 
% by using bisection method

t = (0:N-1)'*2*pi/N;
a = (1 - sqrt(1-ra^2))/ra;
% initial guess using approximation length = sqrt(2)*pi*sqrt(a^2+1)


X0 = [a*cos(t);sin(t)];

[raNew,~,~] = o.geomProp(X0);

cond = abs(raNew - ra)/ra < 1e-4;
maxiter = 10;
iter = 0;
while (~cond && iter < maxiter);
  iter = iter + 1;
  if (raNew > ra)
    a = 0.9 * a;
  else
    a = 1.05 * a;
  end
  % update the major axis
  X0 = [cos(t);a*sin(t)];
  % Compute new possible configuration
  [raNew,~,~] = o.geomProp(X0);
  % compute new reduced area
  cond = abs(raNew - ra) < 1e-2;
  % check if the residual is small enough
end
% iteration quits if reduced area is achieved within 1% or 
% maxiter iterations have been performed


end % ellipse



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = fillDomain(o,Xref,volFrac,Xwalls,fmm)
% fillDomain(X0,volFrac,Xwalls) uses a Monte-Carlo method to fill 
% the domain enclosed by Xwalls with dialations and translations
% of Xref.  Tries to obtain volume fraction volFrac

[xwalls,ywalls] = o.getXY(Xwalls);
% solid walls
[x0,y0] = o.getXY(Xref);
% reference solution
N = numel(x0);
%figure(1); clf; hold on
%plot(xwalls,ywalls,'k')
nvbd = size(xwalls,2);
% number of solid walls
[~,area,~] = o.geomProp(Xwalls);
areaGeom = sum(area);
% Total area occupied by the physical domain
[~,areaVes,~] = o.geomProp(Xref);
% area of reference vesicle

nv = ceil(volFrac*areaGeom/areaVes);
% total number of vesicles to achieve volume fraction
% volFrac

walls = capsules(Xwalls,[],[],0,0);
% build object for walls
radx = 1/2*(max(xwalls(:,1)) - min(xwalls(:,1)));
rady = 1/2*(max(ywalls(:,1)) - min(ywalls(:,1)));
% need radx and rady to decide where to randomly
% sample the geometry.  This way we do not pick 
% a center that is  outside of the solid walls 
% on a regular basis

X = zeros(2*N,nv);
k = 1;
% counter for the number of successfully placed vesicles
while k <= nv 
  cx = 2*radx*(rand-1/2);
  cy = 2*rady*(rand-1/2);
  phi = 2*pi*rand;
  % center
  xpot = cx + x0*cos(phi) + y0*sin(phi);
  ypot = cy - x0*sin(phi) + y0*cos(phi);
  % potential vesicle location

  accept = true; % tentatively accept the vesicle

  vesicle = capsules([xpot;ypot],[],[],[],[]);
  % create capsule with the potential vesicle
  [~,icollisionWall] = o.collision(walls,vesicle,[],[],fmm,0);
%  figure(1); plot(xpot,ypot,'r--')
%  pause
  if icollisionWall
    accept = false;
    % at least one of the vesicles's points is outside of one
    % of the solid wall components
  end

  if accept 
%    figure(1); plot(xpot,ypot,'r')
%    pause
    X(:,k) = [xpot;ypot];
    % if vesicle is not outside of physical walls, accept it as
    % a potential new vesicle.  It will be kept as long as it intersects
    % no other vesicles

    vesicle = capsules(X(:,1:k),[],[],0,0);
    % create an object with the current configuration of
    % vesicles

    icollisionVes = o.collision(...
        vesicle,[],[],[],fmm,0);
    % see if vesicles have crossed using collision detection code
    % Can be used with or without the fmm
    if icollisionVes
      X(:,k) = 0;
      % if they've crossed, reject it
    else
      k = k + 1;
      fprintf('%d Vesicles left to fill the domain to the desired volume fraction\n', nv-k+1)
      % if they haven't crossed, increase the total number of vesicles
    end
  end

end % while k <= nv

end % fillDomain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [icollisionVes,icollisionWall] = ...
    collision(o,vesicle,walls,NearV2V,NearV2W,fmm,near)
% [icollisionVes,icollisionWall] = ...
%    collision(vesicle,walls,NearV2V,NearV2W,fmm,near)
% uses double-layer potential for Laplace's equation to determing if
% vesicles have crossed (stored in icollisionVes) or if solid wall 
% point is inside a vesicle (icollisionWall).
% Can be done with or without the fmm.  vesicle stores all the
% vesicles, walls stores the solid walls (if they exist), NearV2V
% and NearV2W are near-singular integration structures, fmm and near
% are flgas for using the FMM and near-singular integration

%figure(2); clf; hold on
%plot(vesicle.X(1:64,:),vesicle.X(65:end,:),'k')
%figure(1)
%pause
[x,y] = o.getXY(vesicle.X);
[nx,ny] = o.getXY(vesicle.normal);

f = [ones(vesicle.N,vesicle.nv);zeros(vesicle.N,vesicle.nv)];
% Density function is constant.  Pad second half of it with zero
op = poten(vesicle.N);
% load object for doing near-singular integration and evaluating
% laplace double-layer potential

if ~fmm
  kernel = @op.exactLaplaceDL;
else
  kernel = @op.exactLaplaceDLfmm;
end
% kernel for laplace's double layer potential.  Only difference
% is if FMM is used or not

if near
  DLP = zeros(2*vesicle.N,2*vesicle.N,vesicle.nv);
  for k=1:vesicle.nv
    DLP(1:vesicle.N,1:vesicle.N,k) = 0;
    % can cheat here because we know that the double-layer
    % potential applied to our function f will always be 0
    % This won't work if we are considering density functions
    % that are not one everywhere
  end

  Fdlp = op.nearSingInt(vesicle,f,DLP,...
          NearV2V,kernel,vesicle,true);
  % Apply the double-layer laplace potential with constant
  % boundary condition.  Skips self-vesicle term.  This means
  % that if the vesicles have crossed, the Fdlp will achieve
  % a value of 1.  If they have not crossed, it will be 0
  if ~isempty(walls)
    FDLPwall = op.nearSingInt(vesicle,f,DLP,...
            NearV2W,kernel,walls,false);
  end
else
  Fdlp = kernel(vesicle,f);
  if ~isempty(walls)
    [~,FDLPwall] = kernel(vesicle,f,walls.X,(1:vesicle.nv));
    % vesicles are sources and solid walls are targets
    % If vesicles have not crossed a solid wall, values is 0
%    figure(2); clf
%    plot(FDLPwall)
%    figure(1)
  end
end
Fdlp = Fdlp(1:vesicle.N,:);
% this layer-potential is not vector valued, so we throw away the
% second half.

if near
  buffer = 1e-4;
else
  buffer = 1e-1;
  % without near-singular, it is less robust.  This is only being used
  % in fillDomain
end
% can't set buffer too big because near singular integration does not
% assign a value of 1 when near points cross.  This is because I did not
% swtich the direction of the normal for this case.  So, the lagrange
% interpolation points may look something like [0 1 1 1 1 ...] instead
% of [1 1 1 1 1 ...].  The local interpolant is affected by this and
% a value somewhere between 0 and 1 will be returned
% This can cause problems when the vesicles first cross
icollisionVes = any(abs(Fdlp(:)) > buffer);

if ~isempty(walls)
  FDLPwall = FDLPwall(1:walls.N,:);
  icollisionWall = any(abs(FDLPwall(:)-1) > buffer);
%  disp(icollisionWall)
else
  icollisionWall = false;
end

%if nargout == 2
%  figure(3); clf; hold on
%  plot(walls.X(1:walls.N,:),walls.X(walls.N+1:2*walls.N,:),'o')
%  plot(vesicle.X(1:vesicle.N,:),vesicle.X(vesicle.N+1:2*vesicle.N,:),...
%      'r-o')
%  figure(4); clf;
%  plot(abs(FDLPwall(:)-1),'b-o')
%  pause
%end
%icollisionWall = any(x.^2 + y.^2 > 10^2);


end % collision


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = correctAreaAndLength(o,X,a0,l0)
% Xnew = correctAreaAndLength(X,a0,l0) changes the shape of the vesicle
% by finding the shape Xnew that is closest to X in the L2 sense and
% has the same area and length as the original shape

N = size(X,1)/2;
minFun = @(z) 1/N*min(sum((z - X).^2));

%active-set | interior-point | interior-point-convex | levenberg-marquardt | ...
%                           sqp | trust-region-dogleg | trust-region-reflective

options = optimset('Algorithm','sqp','TolCon',1e-8,'TolFun',1e-2,...
    'display','off');
Xnew = fmincon(minFun,X,[],[],[],[],[],[],@(z) o.nonlcon(z,a0,l0),options);

end % correctAreaAndLength


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cIn,cEx] = nonlcon(o,X,a0,l0)
% [cIn,cEx] = nonlcon(X,a0,l0) is the non-linear constraints required
% by fmincon

[~,a,l] = o.geomProp(X);

cIn = [];
% new inequalities in the constraint
cEx = [(a-a0)./a0 (l-l0)./l];
% want to keep the area and length the same

end % nonlcon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta,arcLength] = redistributeArcLength(o,x,y)
% theta = redistributeArcLength(o,x,y) finds a discretization of paramter space
% theta so that the resulting geometry will be equispaced in arclength

N = numel(x);
t = (0:N-1)'*2*pi/N;
[~,~,len] = o.geomProp([x;y]);
% find total perimeter
[Dx,Dy] = o.getDXY([x;y]);
% find derivative
arc = sqrt(Dx.^2 + Dy.^2);
arch = fft(arc);
modes = -1i./[(0:N/2-1) 0 (-N/2+1:-1)]';
modes(1) = 0;
modes(N/2+1) = 0;
arcLength = real(ifft(modes.*arch) - sum(modes.*arch/N) + arch(1)*t/N);
% arclength from 0 to t(n) in parameter space

z1 = [arcLength(end-9:end)-len;arcLength;arcLength(1:10)+len];
z2 = [t(end-9:end)-2*pi;t;t(1:10)+2*pi];
% put in some overlap to account for periodicity

theta = [interp1(z1,z2,(0:N-1)'*len/N,'spline')];

end % redistributeArcLength




end % methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static)

function f = arcDeriv(f,m,sa,IK)
% df = arcDeriv(f,m,s,IK) is the arclength derivative of order m.
% f is a matrix of scalar functions (each function is a column)
% f is assumed to have an arbitrary parametrization
% sa = d s/ d a, where a is the aribtrary parameterization
% IK is the fourier modes which is saved and used to accelerate 
% this routine

for j=1:m
  f = sa.*ifft(IK.*fft(f));
end
% compute the mth order derivative
f = real(f);
%norm(f)

end % arcDeriv

end % methods (Static)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




end % classdef
