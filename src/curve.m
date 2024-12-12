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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y]=getXY(o,X)
% [x,y] = getXY(X) get the [x,y] component of curves X
N = size(X,1)/2;
x = X(1:N,:);
y = X(N+1:end,:);

end % getXY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V = setXY(o,x,y)
% V = setXY(x,y) set the [x,y] component of vector V on the curve
N = size(x,1);
V=zeros(2*N,size(x,2));
V(1:N,:) = x;
V(N+1:end,:) = y;

end % setXY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

N = size(X,1)/2;
nv = size(X,2);
IK = fft1.modes(N,nv);

% get the x y components
[Dx,Dy] = o.getDXY(X);

jacobian = sqrt(Dx.^2 + Dy.^2); 

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

length = sum(sqrt(Dx.^2 + Dy.^2))*2*pi/N;
area = sum(x.*Dy - y.*Dx)*pi/N;

reducedArea = 4*pi*area./length.^2;

end % geomProp



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  x = scale*c*r.*cos(t);
  y = scale*r.*sin(t);
  tNew = o.redistributeArcLength(x,y);
  r = 0.5*sqrt( (a*cos(tNew-theta)).^2 + (b*sin(tNew-theta)).^2) + ...
      .07*cos(12*(tNew-theta));
  x = scale*c*r.*cos(tNew);
  y = scale*r.*sin(tNew);
  X = [x;y];
  % radius of curly vesicle

elseif any(strcmp(options,'cShape'))
  x = [-(1.5+sin(t)).*(cos(0.99*pi*cos(t)))];
  y = [+(1.5+sin(t)).*(sin(0.99*pi*cos(t)))];
  X0 = [x;y];

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
  radius = 1 + 0.2*cos(folds*t) + amp*cos(30*t);
  x = radius.*cos(t);
  y = radius.*sin(t);
  tNew = o.redistributeArcLength(x,y);
  % interpolate inverse function to find new parameter spacing

  radius = 1 + 0.2*cos(folds*tNew) + amp*cos(30*tNew);
  x = radius.*cos(tNew);
  y = radius.*sin(tNew);
  % re-define geometry equi-spaced in arclength

  X = [x;y];
  % a star that comes very close to intersecting itself at the origin

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
  %a = 10; b = 3; c = 0.6; order = 8; %original setup
  %a = 10; b = 3; c = 0; order = 8; %like long tube
  %a = 1.5; b = 3; c = 0; order = 8; %like a thin box
  %a = 10; b = 3; c = 0.78; order = 8;
  %a = 30; b = 3; c = 0.77; order = 8;
  a = 20; b = 3; c = 0.78; order = 8;
  %a = 25; b = 3; c = 0.78; order = 8;
  % parameters for the boundary
  Nsides = ceil(0.5*b/(2*a+2*b)*N);
  Ntop = (N-4*Nsides)/2;
  t1 = linspace(0,0.2*pi,Nsides+1); t1 = t1(1:end-1)';
  t2 = linspace(0.2*pi,pi-0.2*pi,Ntop+1); t2 = t2(1:end-1)';
  t3 = linspace(pi-0.2*pi,pi+0.2*pi,2*Nsides+1); t3 = t3(1:end-1)';
  t4 = linspace(pi+0.2*pi,2*pi-0.2*pi,Ntop+1); t4 = t4(1:end-1)';
  t5 = linspace(2*pi-0.2*pi,2*pi,Nsides+1); t5 = t5(1:end-1)';
  t = [t1;t2;t3;t4;t5];
  % Parameterize t so that geometry is closer to equi-spaced in
  % arclength
%  t = (0:N-1)'*2*pi/N;
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  ind = abs(x) < pi;
  y(ind) = y(ind).*(1-c*cos(x(ind)))/(1+c);
  X0 = [x;y];
  % choked domain.  a and b control the length and height.  c
  % controls the width of the gap, and order controls the
  % regularity
elseif any(strcmp(options,'chokeerf'))
  %a = 25; b = 3; c = 0.1; order = 8; stif = 0.04;
  %a = 25; b = 3; c = 0.125; order = 8; stif = 0.04;
  %a = 25; b = 3; c = 0.065; order = 8; stif = 0.01;
  %a = 25; b = 3; c = 0.065; order = 8; stif = 0.04;lscale = 1; % default test
  %a = 25; b = 3; c = 0.08; order = 8; stif = 0.04; % larger choke for test
  %a = 25; b = 3; c = 0.065; order = 8; stif = 0.16; % choke part more smooth
  %a = 25; b = 3; c = 0.10; order = 8; stif = 0.16; % largest choke part more smooth
  %a = 25; b = 3; c = 0.08; order = 8; stif = 0.16; % choke part more smooth
  a = 20; b = 3; c = 0.065; order = 8; stif = 0.18;lscale = 1.5; 
  if any(strcmp(options,'stiff'))
    ind = find(strcmp(options,'stiff'));
    stif = options{ind+1};
  end
  if any(strcmp(options,'cwide'))
    ind = find(strcmp(options,'cwide'));
    c = options{ind+1};
  end
  if any(strcmp(options,'lscale'))
    ind = find(strcmp(options,'lscale'));
    lscale = options{ind+1};
  end
  % parameters for the boundary
  Nsides = ceil(0.5*b/(2*a+2*b)*N);
  Ntop = (N-4*Nsides)/2;
  t1 = linspace(0,0.2*pi,Nsides+1); t1 = t1(1:end-1)';
  t2 = linspace(0.2*pi,pi-0.2*pi,Ntop+1); t2 = t2(1:end-1)';
  t3 = linspace(pi-0.2*pi,pi+0.2*pi,2*Nsides+1); t3 = t3(1:end-1)';
  t4 = linspace(pi+0.2*pi,2*pi-0.2*pi,Ntop+1); t4 = t4(1:end-1)';
  t5 = linspace(2*pi-0.2*pi,2*pi,Nsides+1); t5 = t5(1:end-1)';
  t = [t1;t2;t3;t4;t5];
  % Parameterize t so that geometry is closer to equi-spaced in
  % arclength
  %  t = (0:N-1)'*2*pi/N;
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  	
  ind = abs(x) < (lscale*pi);
  y(ind) = y(ind).*(-(0.5-c)*(erf(cos(x(ind)/lscale)/sqrt(stif))) + 0.5+c);
  X0 = [x;y];
elseif any(strcmp(options,'stenosis'))
  a = 25; b = 3; c = 0.065; order = 8; stif = 0.04;
  % parameters for the boundary
  Nsides = ceil(0.5*b/(2*a+2*b)*N);
  Ntop = (N-4*Nsides)/2;
  t1 = linspace(0,0.2*pi,Nsides+1); t1 = t1(1:end-1)';
  t2 = linspace(0.2*pi,pi-0.2*pi,Ntop+1); t2 = t2(1:end-1)';
  t3 = linspace(pi-0.2*pi,pi+0.2*pi,2*Nsides+1); t3 = t3(1:end-1)';
  t4 = linspace(pi+0.2*pi,2*pi-0.2*pi,Ntop+1); t4 = t4(1:end-1)';
  t5 = linspace(2*pi-0.2*pi,2*pi,Nsides+1); t5 = t5(1:end-1)';
  t = [t1;t2;t3;t4;t5];
  % Parameterize t so that geometry is closer to equi-spaced in
  % arclength
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);

  indx = abs(x) < pi;
  indy = y < 0;
  ind = and(indx,indy);
  y(ind) = 2*y(ind).*(-(0.5-c)*(erf(cos(x(ind))/sqrt(stif))) + c);
  X0 = [x;y];
elseif any(strcmp(options,'tube'))
  a = 10; b = 3; c = 0; order = 8; %like long tube
  % parameters for the boundary
  Nsides = ceil(0.5*b/(2*a+2*b)*N);
  Ntop = (N-4*Nsides)/2;
  t1 = linspace(0,0.2*pi,Nsides+1); t1 = t1(1:end-1)';
  t2 = linspace(0.2*pi,pi-0.2*pi,Ntop+1); t2 = t2(1:end-1)';
  t3 = linspace(pi-0.2*pi,pi+0.2*pi,2*Nsides+1); t3 = t3(1:end-1)';
  t4 = linspace(pi+0.2*pi,2*pi-0.2*pi,Ntop+1); t4 = t4(1:end-1)';
  t5 = linspace(2*pi-0.2*pi,2*pi,Nsides+1); t5 = t5(1:end-1)';
  t = [t1;t2;t3;t4;t5];
  % Parameterize t so that geometry is closer to equi-spaced in
  % arclength
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  ind = abs(x) < pi;
  y(ind) = y(ind).*(1-c*cos(x(ind)))/(1+c);
  X0 = [x;y];
  % choked domain.  a and b control the length and height.  c
  % controls the width of the gap, and order controls the
  % regularity
elseif any(strcmp(options,'box'))
  a = 1.6; b = 3; c = 0; order = 8; %like a thin box
  %a = 1.5; b = 3; c = 0; order = 8; %like a thin box
  %a = 1.6; b = 3.2; c = 0; order = 8; %like a thin box
  %a = 12; b = 3; c = 0; order = 8; %like a thin box
  % parameters for the boundary
  Nsides = ceil(0.5*b/(2*a+2*b)*N);
  Ntop = (N-4*Nsides)/2;
  t1 = linspace(0,0.2*pi,Nsides+1); t1 = t1(1:end-1)';
  t2 = linspace(0.2*pi,pi-0.2*pi,Ntop+1); t2 = t2(1:end-1)';
  t3 = linspace(pi-0.2*pi,pi+0.2*pi,2*Nsides+1); t3 = t3(1:end-1)';
  t4 = linspace(pi+0.2*pi,2*pi-0.2*pi,Ntop+1); t4 = t4(1:end-1)';
  t5 = linspace(2*pi-0.2*pi,2*pi,Nsides+1); t5 = t5(1:end-1)';
  t = [t1;t2;t3;t4;t5];
  % Parameterize t so that geometry is closer to equi-spaced in
  % arclength
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  ind = abs(x) < pi;
  y(ind) = y(ind).*(1-c*cos(x(ind)))/(1+c);
  X0 = [x;y];
  
  %xtmp = cos(pi/2)*x - sin(pi/2)*y;
  %ytmp = sin(pi/2)*x + cos(pi/2)*y;
  %X0 = [xtmp;ytmp];
  % choked domain.  a and b control the length and height.  c
  % controls the width of the gap, and order controls the
  % regularity
elseif any(strcmp(options,'choke2'))
  a = 40; b = 3; c = 0.6; order = 8;
  % parameters for the boundary
%  Nsides = ceil(0.5*b/(2*a+2*b)*N);
%  Ntop = (N-4*Nsides)/2;
%  t1 = linspace(0,0.2*pi,Nsides+1); t1 = t1(1:end-1)';
%  t2 = linspace(0.2*pi,pi-0.2*pi,Ntop+1); t2 = t2(1:end-1)';
%  t3 = linspace(pi-0.2*pi,pi+0.2*pi,2*Nsides+1); t3 = t3(1:end-1)';
%  t4 = linspace(pi+0.2*pi,2*pi-0.2*pi,Ntop+1); t4 = t4(1:end-1)';
%  t5 = linspace(2*pi-0.2*pi,2*pi,Nsides+1); t5 = t5(1:end-1)';
%  t = [t1;t2;t3;t4;t5];
  t = (0:N-1)'*2*pi/N;
  % Parameterize t so that geometry is closer to equi-spaced in arclength
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
  ind = abs(x) < pi;
  y(ind) = y(ind).*(1-c*cos(x(ind)))/(1+c);

  tNew = redistributeArcLength(o,x,y);
  r = (cos(tNew).^order + sin(tNew).^order).^(-1/order);
  x = a*r.*cos(tNew); y = b*r.*sin(tNew);
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
  
elseif any(strcmp(options,'tube1'))
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
  
elseif any(strcmp(options,'diffuser'))
  a = 10; b = 1; c = 4; order = 8;
  % parameters for the boundary
%  Nsides = ceil(0.5*b/(2*a+2*b)*N);
%  Ntop = (N-4*Nsides)/2;
%  t1 = linspace(0,0.2*pi,Nsides+1); t1 = t1(1:end-1)';
%  t2 = linspace(0.2*pi,pi-0.2*pi,Ntop+1); t2 = t2(1:end-1)';
%  t3 = linspace(pi-0.2*pi,pi+0.2*pi,2*Nsides+1); t3 = t3(1:end-1)';
%  t4 = linspace(pi+0.2*pi,2*pi-0.2*pi,Ntop+1); t4 = t4(1:end-1)';
%  t5 = linspace(2*pi-0.2*pi,2*pi,Nsides+1); t5 = t5(1:end-1)';
%  t = [t1;t2;t3;t4;t5];
  t = (0:N-1)'*2*pi/N;
  % Parameterize t so that geometry is closer to equi-spaced in
  % arclength
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  x = a*r.*cos(t); y = b*r.*sin(t);
%  ind = abs(x) < pi;
%  y(ind) = y(ind) + c/2*sign(y(ind)).*exp(-c./(pi^2 - x(ind).^2));
  ind = (abs(x) < 0.5 & y < 0);
  c = 0.5;
  y(ind) = y(ind) + 10*c/2.*exp(-c./(0.5^2 - x(ind).^2));

  X = [x;y];
  % opposite of chocked domain.  a and b control the length and
  % height.  c controls the width of the gap, and order controls the
  % regularity

elseif any(strcmp(options,'microfluidic'));
  t = (0:N-1)'*2*pi/N;

  a = 2; b = 2; order = 12;
  % parameters for the outer boundary
  r = (cos(t).^order + sin(t).^order).^(-1/order);
  Xout = [a*r.*cos(t);b*r.*sin(t)];

  a = 0.85; b = 0.85; order = 12;
  % parameters for the inner boundary
  r = (cos(-t).^order + sin(-t).^order).^(-1/order);
  Xin1 = [a*r.*cos(-t)-0.95;b*r.*sin(-t)-0.95];
  Xin2 = [a*r.*cos(-t)-0.95;b*r.*sin(-t)+0.95];
  Xin3 = [a*r.*cos(-t)+0.95;b*r.*sin(-t)-0.95];
  Xin4 = [a*r.*cos(-t)+0.95;b*r.*sin(-t)+0.95];

  X = [Xout Xin1 Xin2 Xin3 Xin4];

else
  if ~isempty(ra)
    X0 = o.ellipse(N,ra);
    % build a vesicle of reduced area ra with N points
  else
    a = 1; b = 4.8*a; c = 1; 
    r = 0.5*sqrt( (a*cos(t)).^2 + (b*sin(t)).^2);
    x = c*r.*cos(t); y = r.*sin(t);
    tNew = o.redistributeArcLength(x,y);
    r = 0.5*sqrt( (a*cos(tNew)).^2 + (b*sin(tNew)).^2);
    x = c*r.*cos(tNew); y = r.*sin(tNew);
    X0 = [x;y];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
area = abs(area);
areaGeom = sum(area);
% Total area occupied by the physical domain
[~,areaVes,~] = o.geomProp(Xref);
% area of reference vesicle

nv = ceil(volFrac*areaGeom/areaVes);
% total number of vesicles to achieve volume fraction
% volFrac

walls = capsules(Xwalls,[],[],0,0,false);
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
xcen = [];
ycen = [];
angcen = [];
while k <= nv 
  cx = 2*radx*(rand-1/2);
  cy = 2*rady*(rand-1/2);
  phi = 2*pi*rand;
  % center
  xpot = cx + x0*cos(phi) - y0*sin(phi);
  ypot = cy + x0*sin(phi) + y0*cos(phi);
  % potential vesicle location

  accept = false; % tentatively accept the vesicle

  vesicle = capsules([xpot;ypot],[],[],[],[],false);
  % create capsule with the potential vesicle
  [~,icollisionWall] = vesicle.collision(walls,vesicle,[],[],fmm,0);
  if icollisionWall
    accept = true;
    % at least one of the vesicles's points is outside of one
    % of the solid wall components
  end

  if accept 
    X(:,k) = [xpot;ypot];
    % if vesicle is not outside of physical walls, accept it as
    % a potential new vesicle.  It will be kept as long as it intersects
    % no other vesicles

    vesicle = capsules(X(:,1:k),[],[],0,0,false);
    % create an object with the current configuration of
    % vesicles

    icollisionVes = vesicle.collision(...
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
      xcen = [xcen;cx];
      ycen = [ycen;cy];
      angcen = [angcen;phi];
    end
  end

end % while k <= nv
xcen
ycen
angcen
end % fillDomain


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = addAndRemove(o,X,walls,near,fmm)
% X = addAndRemove(X,walls,near,fmm) is used to fudge in periodicity.
% When a vesicle passes a certain threshold (hard-coded into this
% function), it randomly places it at the start of the computational
% domain with the same shape.  If a collision with another vesicle
% results, it tries again.  If too many attempts are made, it simply
% continues with its current position and tries again at the next
% iteration

near = false;
[x,y] = o.getXY(X);
[N,nv] = size(x);

ind = find(max(x) > 8);
if numel(ind) > 0
  [~,s] = max(max(x(:,ind)));
  % find the index of the point farthest to the right
  x = x(:,[1:ind(s)-1 (ind(s)+1:nv)]);
  y = y(:,[1:ind(s)-1 (ind(s)+1:nv)]);
  Xminus = o.setXY(x,y);
  % all points except the one that is farthest to the right

  icollisionVes = true;
  maxTrys = 10;
  % maximum number of attempts before attempt to build in peridocity is
  % postponed
  ntrys = 0;
  while (icollisionVes && ntrys < maxTrys)
    ntrys = ntrys + 1;
%    center = [-8.5;0.2*rand(1)];
    center = [-8;0.7*2*(rand(1)-0.5)];
    % random center
    Xadd = X(:,ind(s));
    Xadd(1:end/2) = Xadd(1:end/2) - mean(Xadd(1:end/2));
    Xadd(1+end/2:end) = Xadd(1+end/2:end) - mean(Xadd(1+end/2:end));
    Xadd(1:end/2) = Xadd(1:end/2) + center(1);
    Xadd(1+end/2:end) = Xadd(1+end/2:end) + center(2);
    if max(abs(Xadd(1+end/2:end,:)) > 9.5e-1)
      icollisionWall = true;
    else
      icollisionWall = false;
    end
    % solid wall is at about 0.95 at the point -8.5
    if ~icollisionWall
      % parameterization of new vesicle
      vesicle = capsules([Xminus Xadd],[],[],0,0,o.antiAlias);
      % class corresponding to the new configuration
      if near
        [NearV2V,NearV2W] = vesicle.getZone(walls,3);
      else
        NearV2V = [];
        NearV2W = [];
      end
      icollisionVes = o.collision(vesicle,walls,...
          NearV2V,NearV2W,fmm,near);
      % check if the vesicles have crossed
    end
  end

  if ntrys < maxTrys
    X = vesicle.X;
  end
  % if it's tried to many times, simply move on and try again at the
  % next iteration
end

end % addAndRemove

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xnew = correctAreaAndLength(o,X,a0,l0)
% Xnew = correctAreaAndLength(X,a0,l0) changes the shape of the vesicle
% by finding the shape Xnew that is closest to X in the L2 sense and
% has the same area and length as the original shape

N = size(X,1)/2;

options = optimset('Algorithm','sqp','TolCon',1e-8,'TolFun',1e-4,...
    'display','off');

Xnew = zeros(size(X));
for k = 1:size(Xnew,2);
  minFun = @(z) 1/N*min(sum((z - X(:,k)).^2));
  Xnew(:,k) = fmincon(minFun,X(:,k),[],[],[],[],[],[],...
      @(z) o.nonlcon(z,a0(k),l0(k)),options);
end
% Looping over vesicles, correct the area and length of each vesicle

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
function [X,u,sigma] = redistributeParameterize(o,X,u,sigma)
% [X,u,sigma] = resdistributeParameterize(o,X,u,sigma) resdistributes
% the vesicle shape eqiuspaced in arclength and adjusts the tension and
% velocity according to the new parameterization

N = size(X,1)/2;
nv = size(X,2);
modes = [(0:N/2-1) (-N/2:-1)];
jac = o.diffProp(X);
for k = 1:nv
  if norm(jac(:,k) - mean(jac(:,k)),inf) > 0
    theta = o.redistributeArcLength(X(1:end/2,k),...
        X(end/2+1:end,k));
    zX = X(1:end/2,k) + 1i*X(end/2+1:end,k);
    zu = u(1:end/2,k) + 1i*u(end/2+1:end,k);
    zXh = fft(zX)/N;
    zuh = fft(zu)/N;
    sigmah = fft(sigma(:,k))/N;
    zX = zeros(N,1);
    zu = zeros(N,1);
    sigma(:,k) = zeros(N,1);
    for j = 1:N
      zX = zX + zXh(j)*exp(1i*modes(j)*theta);
      zu = zu + zuh(j)*exp(1i*modes(j)*theta);
      sigma(:,k) = sigma(:,k) + sigmah(j)*exp(1i*modes(j)*theta);
    end
    sigma = real(sigma);
    % redistribute the vesicle positions and tension so that it is
    % equispaced in arclength

    X(:,k) = o.setXY(real(zX),imag(zX));
    u(:,k) = o.setXY(real(zu),imag(zu));
  end
end

end % redistributeParameterize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta,arcLength] = redistributeArcLength(o,x,y)
% theta = redistributeArcLength(o,x,y) finds a discretization of
% paramter space theta so that the resulting geometry will be
% equispaced in arclength

uprate = 1;
N = numel(x);
Nup = uprate*N;
t = (0:Nup-1)'*2*pi/Nup;
x = interpft(x,Nup);
y = interpft(y,Nup);
[~,~,len] = o.geomProp([x;y]);
% find total perimeter
[Dx,Dy] = o.getDXY([x;y]);
% find derivative
arc = sqrt(Dx.^2 + Dy.^2);
arch = fft(arc);
modes = -1i./[(0:Nup/2-1) 0 (-Nup/2+1:-1)]';
modes(1) = 0;
modes(Nup/2+1) = 0;
arcLength = real(ifft(modes.*arch) - sum(modes.*arch/Nup) + ...
    arch(1)*t/Nup);
% arclength from 0 to t(n) in parameter space

z1 = [arcLength(end-6:end)-len;arcLength;arcLength(1:7)+len];
z2 = [t(end-6:end)-2*pi;t;t(1:7)+2*pi];
% put in some overlap to account for periodicity

theta = [interp1(z1,z2,(0:N-1)'*len/N,'spline')];

end % redistributeArcLength

end % methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = arcDeriv(f,m,sa,IK)
% f = arcDeriv(f,m,s,IK,col) is the arclength derivative of order m.
% f is a matrix of scalar functions (each function is a column)
% f is assumed to have an arbitrary parametrization
% sa = d s/ d a, where a is the aribtrary parameterization
% IK is the fourier modes which is saved and used to accelerate 
% this routine

for j=1:m
  %f = sa.*ifft(IK.*fft(f));
  f = sa.*real(ifft(IK.*fft(f)));
end
f = real(f);

end % arcDeriv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zCoarse = restrict(z,Nfine,Ncoarse)
% zCoarse = restrict(z,Nfine,Ncoarse) restricts the periodic function
% z.  z has a column vector containing arbitrarly many periodic
% functions of size Nfine.  The output, zCoarse, has the same number of
% periodic copies at a grid with Ncoarse points.  NOTE: Nfine/Ncoarse
% must be a power of 2 if method == 'local'

method = 'spectral';
%method = 'local';

if Nfine == Ncoarse
  zCoarse = z;
  % if coarse and fine grids are the same, nothing to do
else
  nSecs = size(z,1)/Nfine;
  % number of functions that need to be restricted.  This is the number
  % of periodic functions that are stacked in each row

  if strcmp(method,'spectral')
    zFine = zeros(Nfine*nSecs,size(z,2));
    for j = 1:nSecs
      istart = (j-1)*Nfine + 1;
      iend = istart + Nfine - 1;
      istart2 = (j-1)*Ncoarse + 1;
      iend2 = istart2 + Ncoarse - 1;
      zh = fft(z(istart:iend,:));
      % take fft of a matrix of periodic columns
      zh = [zh(1:Ncoarse/2,:);zh(end-Ncoarse/2+1:end,:)]*Ncoarse/Nfine;
      % zero all the high frequencies and scale appropriatly
      zCoarse(istart2:iend2,:) = real(ifft(zh));
      % move to physical space and take real part as we are assuming
      % input is real
    end
  end
  % spectral restriction

  if strcmp(method,'local')
    nRestrict = log2(Nfine/Ncoarse);
    % number of halvings to do
    for k = 1:nRestrict
      Nfine = Nfine/2;
      zCoarse = zeros(size(z,1)/2,size(z,2));
      % create grid half the size
      for j = 1:nSecs
        istart = (j-1)*Nfine + 1;
        iend = istart + Nfine - 1;
        zCoarse(istart:iend,:) = 1/2*z(2*(istart:iend)-1,:);
        zCoarse(istart:iend,:) = zCoarse(istart:iend,:) + ...
          1/4*z(2*(istart:iend),:);
        zCoarse(istart+1:iend,:) = zCoarse(istart+1:iend,:) + ...
          1/4*z(2*(istart+1:iend)-2,:);
        zCoarse(istart,:) = zCoarse(istart,:) + 1/4*z(2*iend,:);
        % usual (1/4,1/2,1/4) restriction with periodicity built in
      end
      z = zCoarse;
    end
  end
  % local restriction
end

end % restrict

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zFine = prolong(z,Ncoarse,Nfine)
% zFine = prolong(z,Ncoarse,Nfine) prolongs the periodic function z.  z
% has a column vector containing arbitrarly many periodic functions of
% size Ncoarse.  The output, zfine, has the same number of periodic
% copies at a grid with Nfine points.  NOTE: Nfine/Ncoarse must be a
% power of 2

method = 'spectral';
%method = 'local';

if Nfine == Ncoarse
  zFine = z;
% if coarse and fine grids are the same, nothing to do
else
  nSecs = size(z,1)/Ncoarse;
  % number of functions that need to be prolonged.  This is the number of
  % periodic functions that are stacked in each row

  if strcmp(method,'spectral')
    zFine = zeros(Nfine*nSecs,size(z,2));
    for j = 1:nSecs
      istart = (j-1)*Ncoarse + 1;
      iend = istart + Ncoarse - 1;
      istart2 = (j-1)*Nfine + 1;
      iend2 = istart2 + Nfine - 1;

      zFine(istart2:iend2,:) = interpft(z(istart:iend,:),Nfine);
    end
  end
  % spectral prolongation

  if strcmp(method,'local')
    nProlong = log2(Nfine/Ncoarse);
    % number of doublings to do

    for k = 1:nProlong
      Ncoarse = Ncoarse*2;
      zFine = zeros(size(z,1)*2,size(z,2));
      for j = 1:nSecs
        istart = (j-1)*Ncoarse + 1;
        iend = istart + Ncoarse - 1;
        zFine(istart:2:iend,:) = z((istart+1)/2:iend/2,:);
        zFine(istart+1:2:iend,:) = 1/2*z((istart+1)/2:iend/2,:);
        zFine(istart+1:2:iend-2,:) = zFine(istart+1:2:iend-2,:) + ...
            1/2*z(((istart+1)/2:(iend-2)/2)+1,:);
        zFine(iend,:) = zFine(iend,:) + 1/2*z((istart+1)/2,:);
        % usual (1/2,1,1/2) prolongation with periodicity built in
      end
      z = zFine;
    end
  end
  % local prolongation
end

end % prolong

end % methods (Static)

end % classdef
