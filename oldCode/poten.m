classdef poten 
% this class defines single and double layers for various kernels 
% (stokes, laplace) on 2D periodic curves.  Also defines the 
% integrals required for pressure and stress.
% Defines the matricies that map a density function defined on the
% boundary of a curve to the layer potential evaluated on the 
% curve, and defines the operators that takes a density function 
% defined on the boundary of a curve and returns the layer 
% potential at arbitrary target points.
% This class also has the main routine that evaluates layer
% potentials using near-singular integration.
    
properties
  N; % points per curve
  qw; % quadrature weights for logarithmic singularity
  qp; % quadrature points for logarithmic singularity (Alpert's rule)
  interpMat;  
  % interpolation matrix used for near-singular integration
  % This matrix replaces the need to use polyfit
end % properties

methods 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = poten(N)
% o = poten(N): constructor; N is the number of points per curve.
% initialize class.

o.N = N;
o.interpMat = o.lagrangeInterp;
% load in the interpolation matrix which is precomputed
% with 7 interpolation points
accuracyOrder = 8;
%qw = o.quadratureSold(accuracyOrder);
qw = o.quadratureS(accuracyOrder);
o.qp = qw(:,2:end);
o.qw = qw(:,1);

end % poten: constructor
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = Lagrange(o,order);
% A = Lagrange(order) generates that matrix that takes order
% function values and returns the coefficients of the polynomial
% of order-1.  Points are assumed to be equally spaced in [0,1] 

A = zeros(order);
for j=1:order
  A(:,j) = ((0:order-1)/(order-1)).^(order-j);
end
A = inv(A);

end % Lagrange


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function qw = quadratureS(o,q);
% qw = quadratureS(q) generates the quadrature rules for a function
% with o.N points and a logarithmic singularity at the origin.  q
% controls the accuracy.  All rules are from Alpert 1999.  This is
% Bryan's reformulation which uses Alpert's quadrature rules as
% described in Section 7, but with Fourier interpolation

[v,u,a] = o.getWeights(q);
% get the weights coming from Table 8 of Alpert's 1999 paper

h = 2*pi/o.N;
n = o.N - 2*a + 1;

of = fft1;
A1 = of.sinterpS(o.N,v*h);
A2 = of.sinterpS(o.N,2*pi-flipud(v*h));
yt = h*(a:n-1+a)';
% regular points away from the singularity
wt = [h*u; h*ones(length(yt),1); h*flipud(u)]/4/pi;
% quadrature points away from singularity

B = sparse(length(yt),o.N);
pos = 1 + (a:n-1+a)';

for k = 1:length(yt)
  B(k, pos(k)) = 1;
end
A = [sparse(A1); B; sparse(A2)];
qw = [wt, A];

end % quadratureS



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function qw = quadratureSold(o,q);
% qw = quadratureS(q) generates the quadrature rules for a function
% with m points and a logarithmic singularity at the origin.  q 
% controls the accuracy.  All rules are from Alpert 1999.

[v,u,a] = o.getWeightsOld(q,2);
[x,w] = o.regularWeights(a);

m = o.N/2;
n = m-2*a+1;
h = 1/m;

evalpots1 = v*h;
evalpots3 = 1-x*h;

ys = pi*[evalpots1;evalpots3];
wt = pi*h*[u;w];

wt = [wt; flipud(wt)]; ys = [ys; 2*pi-flipud(ys)];
of = fft1;
A = of.sinterpS(2*m, ys); 
h = pi/m; 
yt = [a*h:h:(m - a)*h]';
%regular points away from singularity
wt = [wt; h*ones(2*length(yt),1)]/4/pi;
%quadrature weights at regular points away from singularity
lyt = 2*length(yt); 
B = sparse(lyt, 2*m); 
pos = 1+[(a:m-a)'; (m+a:2*m-a)'];

for k = 1:lyt
  B(k, pos(k)) = 1;
end
A = [sparse(A); B];
qw = [wt, A];


end % quadratureS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v,u,a] = getWeightsOld(o,q,ker)
% [v,u,a] = getWeightsOld(q,ker) loads quadrature rules for different
% types of singularities.  All rules come from Bradley Alpert's
% papers.  We are interested in nodesLogx.dat;

switch ker
case 0
  xp = load('nodesr.dat');
  lenth = [2;4;8];
  par = [2;4;7];
case 1
  xp = load('nodes_sqrtx.dat');
  lenth = [4;8;16];
  par = [3;5;10];
case 2
  xp = load('nodesLog.dat');
  lenth = [3;7;15];
  par = [2;5;10];
end

switch q
case 4
  v = xp(1:lenth(1), 1);
  u = xp(1:lenth(1), 2);
  a = par(1);
case 8
  v = xp(1+lenth(1):lenth(1)+lenth(2),1);
  u = xp(1+lenth(1):lenth(1)+lenth(2),2);
  a = par(2);
case 16
  v = xp(1+lenth(2)+lenth(1):sum(lenth),1);
  u = xp(1+lenth(2)+lenth(1):sum(lenth),2);
  a = par(3);
end


end % getWeightsOld

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v,u,a] = getWeights(o,q)
% [v,u,a] = getWeights(q) loads quadrature rules for different types of
% singularities.  All rules come from Bradley Alpert's papers.  We are
% interested in nodesLogx.dat

xp = o.nodesLog;
lenth = [3;7;15];
par = [2;5;10];

switch q
case 4
  v = xp(1:lenth(1), 1);
  u = xp(1:lenth(1), 2);
  a = par(1);
case 8
  v = xp(1+lenth(1):lenth(1)+lenth(2),1);
  u = xp(1+lenth(1):lenth(1)+lenth(2),2);
  a = par(2);
case 16
  v = xp(1+lenth(2)+lenth(1):sum(lenth),1);
  u = xp(1+lenth(2)+lenth(1):sum(lenth),2);
  a = par(3);
end


end % getWeights

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,w] = regularWeights(o,a)
% [x,w] = regularWeights(a) gets quadrature rules
% for regular functions (Alpert's quadrature rules).  These are
% needed for integrating away from singularities

par = [2;3;4;5;7;10];
for i = 1:length(par)
  if par(i)==a
    key = i;
  end
end

lenth = [2;3;4;6;8;12];
xp = load('nodesRegular.dat');
if key==1
  starting = 1;
else
  starting = sum(lenth(1:key-1))+1; 
end

x = xp( starting:starting+lenth(key)-1,1);
w = xp(starting:starting+lenth(key)-1,2);

end % regularWeights


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G = stokesSLmatrix(o,vesicle)
% G = stokesSLmatrix(vesicle) generate single layer potential for
% Stokes vesicle is a data structure defined as in the curve class
% G is (2N,2N,nv) array where N is the number of points per curve and 
% nv is the number of curves in X 

X = vesicle.X; % position of vesicle
oc = curve;
[x,y] = oc.getXY(X);
% seperate x and y coordinates of vesicle
sa = vesicle.sa; % Jacobian

G = zeros(2*o.N,2*o.N,vesicle.nv);
% initalize single-layer potential to zero
for k=1:vesicle.nv  % Loop over curves
  for j=1:o.N % Loop over targets
    ind = 1 + mod(j-1 + (0:o.N-1),o.N);
    xin = o.qp*x(ind,k);
    yin = o.qp*y(ind,k);

    rho = (xin-x(j,k)).^2 + (yin-y(j,k)).^2;
    br = 1./rho;

    logpart = -1/2 * o.qw .* log(rho);
    
    G(j, ind,   k) = (logpart+...
        (o.qw.*( ( x(j,k) - xin).^2 .* br)))'*o.qp;
    G(j, o.N+ind, k) = ...
        (o.qw.*(x(j,k)-xin).*(y(j,k)-yin) .* br)'*o.qp;
    G(j+o.N, 1:o.N ,k) = G(j,o.N+1:2*o.N,k);
    G(j+o.N,o.N+ind,k) = (logpart+...
        (o.qw.*( ( y(j,k) - yin).^2 .* br)))'*o.qp;
  end % j
  sak = repmat(sa(:,k)',2*o.N,2);
  
  G(:,:,k) = G(:,:,k).*sak;
  % multiply both components of G by the Jacobian
  % don't need to divide by 4*pi as it is store in o.qw
end % k


end % stokesSLmatrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Gtrap = stokesSLtrapmatrix(o,vesicle)
% Gtrap = stokesSLtrapmatrix(vesicle) generates double-layer potential
% for Stokes vesicle is a data structure defined as in the capsules
% class Gtrap is (2N,2N,nv) array where N is the number of points per
% curve and nv is the number of curves in X 

Nup = vesicle.N*ceil(sqrt(vesicle.N));
% upsample to N^(3/2).  
if Nup < 1024
  X = vesicle.X; % Vesicle positions
  Xup = [interpft(X(1:vesicle.N,:),Nup);...
         interpft(X(vesicle.N+1:2*vesicle.N,:),Nup)];

  oc = curve;
  [x,y] = oc.getXY(Xup);

  Gtrap = zeros(2*Nup,2*Nup,vesicle.nv);
  for k=1:vesicle.nv  % Loop over curves
    for j=1:Nup % Loop over targets
      rho2 = (x(j,k)-x(:,k)).^2 + (y(j,k)-y(:,k)).^2;
      rho2(j) = 1;
      br = 1./rho2;
      % Set diagonal term to one to avoid dividing by zero
      logpart = -1/2*log(rho2);

      Gtrap(j,1:Nup,k) = logpart + (x(j,k) - x(:,k)).^2 .* br;
      Gtrap(j,Nup+1:2*Nup,k) = (x(j,k)-x(:,k)).*(y(j,k)-y(:,k)) .* br;
      Gtrap(j+Nup,1:Nup,k) = Gtrap(j,Nup+1:2*Nup,k);
      Gtrap(j+Nup,1+Nup:2*Nup,k) = logpart + (y(j,k) - y(:,k)).^2 .* br;
    end % j
  end % k

  Gtrap = Gtrap/4/pi;
  % constant that multiplies for viscosity contrast.  For solid walls, 
  % simply set viscCont to 0 to obtain the right scaling
else
  Gtrap = [];
end


end % stokesSLtrapmatrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = stokesDLmatrix(o,vesicle)
% D = stokesDLmatrix(vesicle), generate double-layer potential for 
% Stokes vesicle is a data structure defined as in the capsules class
% D is (2N,2N,nv) array where N is the number of points per curve and 
% nv is the number of curves in X 

oc = curve;
[x,y] = oc.getXY(vesicle.X);
% Vesicle positions

normal = [vesicle.xt(vesicle.N+1:2*vesicle.N,:);...
         -vesicle.xt(1:vesicle.N,:)]; % Normal vector

D = zeros(2*vesicle.N,2*vesicle.N,vesicle.nv);
for k=1:vesicle.nv  % Loop over curves
  for j=1:vesicle.N % Loop over targets
    rho2 = (x(j,k)-x(:,k)).^2 + (y(j,k)-y(:,k)).^2;
    rho2(j) = 1;
    % Set diagonal term to one to avoid dividing by zero

    coeff = ((x(j,k) - x(:,k)).*normal(1:vesicle.N,k) + ...
        (y(j,k) - y(:,k)).*normal(vesicle.N+1:2*vesicle.N,k)).*...
        vesicle.sa(:,k)./rho2.^2/pi;
    % part of kernel of the double-layer potential

    D(j,:,k) = 2*pi/vesicle.N*[coeff.*(x(j,k) - x(:,k)).^2; ...
      coeff.*(x(j,k) - x(:,k)).*(y(j,k) - y(:,k))]';
    D(j+vesicle.N,:,k) = 2*pi/vesicle.N*[...
      coeff.*(y(j,k) - y(:,k)).*(x(j,k)-x(:,k)); ...
      coeff.*(y(j,k) - y(:,k)).^2]';
    % Build double-layer potential matrix D

    rc = [j j+vesicle.N];
    D(rc,rc,k) = -2*pi/vesicle.N*vesicle.sa(j,k)*...
      vesicle.cur(j,k)/2/pi*...
      [vesicle.xt(j,k)^2 vesicle.xt(j,k)*vesicle.xt(j+vesicle.N,k);...
      vesicle.xt(j+vesicle.N,k)*vesicle.xt(j,k) vesicle.xt(j+vesicle.N,k)^2];
    % Diagonal term requires the limiting value.  Otherwise, above
    % formula would divide by zero
  end % j
  D(:,:,k) = D(:,:,k)*(1-vesicle.viscCont(k));
  % constant that multiplies for viscosity contrast.  For solid walls,
  % simply set viscCont to 0 to obtain the right scaling
end % k

end % stokesDLmatrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N0 = stokesN0matrix(o,vesicle)
% N0 = stokesN0matrix(vesicle) generates the the integral 
% operator with kernel normal(x) \otimes normal(y) which removes
% the rank one defficiency of the double-layer potential.  Need
% this operator for solid walls

N = vesicle.N; % Number of points per vesicle
nv = vesicle.nv; % Number of vesicles
X = vesicle.X; % Vesicle positions
oc = curve;
[x,y] = oc.getXY(X);

normal = vesicle.normal; % Normal vector
sa = vesicle.sa; % Jacobian


N0 = zeros(2*N,2*N,nv);
for k = 1:1 % Loop over curves
  % Only want to form N0 for the outer boundary.
  % Rotlets and Stokeslets take care of null space for inner 
  % boundaries
  for j = 1:2*N % Loop over targets
    N0(j,:,k) = normal(j,k)*normal(:,k).*...
        [sa(:,k);sa(:,k)]*2*pi/N;
    % compute the columns of the modification to the double-layer
    % potential.
  end
end
% Use N0 if solving (-1/2 + DLP)\eta = f where f has no flux through
% the boundary.  By solving (-1/2 + DLP + N0)\eta = f, we guarantee
% that \eta also has no flux through the boundary.  This is not
% required, but it means we're enforcing one addition condition on eta
% which removes the rank one kernel.  DLP is the double-layer potential
% for stokes equation

end % stokesN0matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G = laplaceSLmatrix(o,vesicle)
% G = laplaceSLmatrix(vesicle), generate single layer potential for
% laplace.  vesicle is a data structure defined as in the capsules class
% G is (N,N,nv) array where N is the number of points per curve and 
% nv is the number of curves in X

N = vesicle.N; % Number of points per vesicle
nv = vesicle.nv; % Number of vesicles
X = vesicle.X; % Vesicle positions
oc = curve;
[x,y] = oc.getXY(X);

sa = vesicle.sa; % Jacobian

G = zeros(o.N,o.N,vesicle.nv);
% initalize single-layer potential to zero
for k=1:vesicle.nv  % Loop over curves
  for j=1:o.N % Loop over targets
    ind = 1 + mod(j-1 + (0:o.N-1),o.N);
    xin = o.qp*x(ind,k);
    yin = o.qp*y(ind,k);

    rho2 = (xin-x(j,k)).^2 + (yin-y(j,k)).^2;

    logpart = -1/2 * o.qw .* log(rho2);
    G(j,ind,k) = logpart'*o.qp;
    
  end % j
  sak = repmat(sa(:,k)',o.N,1);
  
  G(:,:,k) = G(:,:,k).*sak;
  % multiply G by the Jacobian
end % k


end % stokesSLmatrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = laplaceDLmatrix(o,vesicle)
% D = laplaceDLmatrix(vesicle), generate double layer potential for
% laplace.  vesicle is a data structure defined as in the capsules class
% D is (N,N,nv) array where N is the number of points per curve and 
% nv is the number of curves in X.  Need this for collision detection

oc = curve;
[x,y] = oc.getXY(vesicle.X);
% Vesicle positions

normal = [vesicle.xt(vesicle.N+1:2*vesicle.N,:);...
         -vesicle.xt(1:vesicle.N,:)]; % Normal vector

D = zeros(vesicle.N,vesicle.N,vesicle.nv);
% initialize double-layer potential to zero
for k=1:vesicle.nv  % Loop over curves
  for j=1:vesicle.N % Loop over targets
    rho2 = (x(:,k)-x(j,k)).^2 + (y(:,k)-y(j,k)).^2;
    rho2(j) = 1;
    % Set diagonal term to one to avoid dividing by zero

    coeff = ((x(:,k) - x(j,k)).*normal(1:vesicle.N,k) + ...
        (y(:,k) - y(j,k)).*normal(vesicle.N+1:2*vesicle.N,k)).*...
        vesicle.sa(:,k)./rho2/2/pi;
    % kernel of double-layer potential for Laplace
    coeff(j) = 0.5*vesicle.cur(j,k)/2/pi*vesicle.sa(j,k);
    % diagonal term is proportional to the curvature
    D(j,:,k) = 2*pi/vesicle.N*coeff;
    % multiply by d\theta term

  end % j
end % k


end % laplaceDLmatrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = pressSLmatrix(o,vesicle)
% P = pressSLmatrix(vesicle), generates the matrix that returns
% the pressure given the traction jump.  Matrix has dimensions
% (N,2*N,nv) where N is the number of points per curve and nv
% is the number of curves in X.  Matrix is not square since
% traction jump is vector-valued whereas the pressure is a 
% scalar-valued function

oc = curve;
[x,y] = oc.getXY(vesicle.X);

P = zeros(vesicle.N,2*vesicle.N,vesicle.nv);
% initialize double-layer potential to zero
for k=1:vesicle.nv  % Loop over curves
  index1 = (1:2:vesicle.N)'; % odd-indexed source points
  for j=2:2:vesicle.N % Loop over targets
    rho2 = (x(j,k) - x(index1,k)).^2 + (y(j,k) - y(index1,k)).^2;
    % distance squared
    coeff = (x(j,k) - x(index1,k)).*vesicle.sa(index1,k)./rho2/2/pi;
    P(j,index1,k) = 4*pi/vesicle.N*coeff;
    % part that multiplies x-component of traction jump
    coeff = (y(j,k) - y(index1,k)).*vesicle.sa(index1,k)./rho2/2/pi;
    P(j,N+index1,k) = 4*pi/vesicle.N*coeff;
    % part that multiplies y-component of traction jump
    % need factor of 4 instead of 2 because we are only using
    % half the terms in the quadrature
  end % j

  index1 = (2:2:vesicle.N)'; % even-indexed source points
  for j=1:2:vesicle.N % Loop over targets
    rho2 = (x(j,k) - x(index1,k)).^2 + (y(j,k) - y(index1,k)).^2;
    % distance squared
    coeff = (x(j,k) - x(index1,k)).*vesicle.sa(index1,k)./rho2/2/pi;
    P(j,index1,k) = 4*pi/vesicle.N*coeff;
    % part that multiplies x-component of traction jump
    coeff = (y(j,k) - y(index1,k)).*vesicle.sa(index1,k)./rho2/2/pi;
    P(j,vesicle.N+index1,k) = 4*pi/vesicle.N*coeff;
    % part that multiplies y-component of traction jump
    % need factor of 4 instead of 2 because we are only using
    % half the terms in the quadrature
  end % j
end % k


end % pressSLmatrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = pressDLmatrix(o,vesicle)
% P = pressDLmatrix(vesicle), generates the matrix that returns
% the pressure given the traction jump.  Matrix has dimensions
% (N,2*N,nv) where N is the number of points per curve and nv
% is the number of curves in X.  Matrix is not square since
% traction jump is vector-valued whereas the pressure is a 
% scalar-valued function

oc = curve;
[x,y] = oc.getXY(vesicle.X);
nx = vesicle.xt(vesicle.N+1:2*vesicle.N,:);
ny = -vesicle.xt(1:vesicle.N,:);
% Normal vector

P = zeros(vesicle.N,2*vesicle.N,vesicle.nv);
% initialize double-layer potential to zero
for k=1:vesicle.nv  % Loop over curves
  index1 = (1:2:vesicle.N)'; % odd-indexed source points
  for j=2:2:vesicle.N % Loop over targets
    rho2 = (x(j,k) - x(index1,k)).^2 + (y(j,k) - y(index1,k)).^2;
    % distance squared
    rx = x(j,k) - x(index1,k);
    ry = y(j,k) - y(index1,k);
    rdotn = rx.*nx(index1,k) + ry.*ny(index1,k);

    coeff = (-nx(index1,k)./rho2 + 2*rdotn./rho2.^2.*rx) .* ...
        vesicle.sa(index1,k)/pi;
    P(j,index1,k) = 4*pi/vesicle.N*coeff;
    % part that multiplies x-component of traction jump
    % need factor of 4 instead of 2 because we are only using
    % half the terms in the quadrature
    P(j,j,k) = P(j,j,k) - sum(coeff)*4*pi/vesicle.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2
    coeff = (-ny(index1,k)./rho2 + 2*rdotn./rho2.^2.*ry) .* ...
        vesicle.sa(index1,k)/pi;
    P(j,vesicle.N+index1,k) = 4*pi/vesicle.N*coeff;
    % part that multiplies y-component of traction jump
    % need factor of 4 instead of 2 because we are only using
    % half the terms in the quadrature
    P(j,vesicle.N+j,k) = P(j,vesicle.N+j,k) - ...
        sum(coeff)*4*pi/vesicle.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2
  end % j

  index1 = (2:2:vesicle.N)'; % even-indexed source points
  for j=1:2:vesicle.N % Loop over targets
    rho2 = (x(j,k) - x(index1,k)).^2 + (y(j,k) - y(index1,k)).^2;
    % distance squared
    rx = x(j,k) - x(index1,k);
    ry = y(j,k) - y(index1,k);
    rdotn = rx.*nx(index1,k) + ry.*ny(index1,k);
    % dot product of r with normal

    coeff = (-nx(index1,k)./rho2 + 2*rdotn./rho2.^2.*rx) .* ...
        vesicle.sa(index1,k)/pi;
    P(j,index1,k) = 4*pi/vesicle.N*coeff;
    % part that multiplies x-component of traction jump
    % need factor of 4 instead of 2 because we are only using
    % half the terms in the quadrature
    P(j,j,k) = P(j,j,k) - sum(coeff)*4*pi/vesicle.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2

    coeff = (-ny(index1,k)./rho2 + 2*rdotn./rho2.^2.*ry) .* ...
        vesicle.sa(index1,k)/pi;
    P(j,vesicle.N+index1,k) = 4*pi/vesicle.N*coeff;
    % part that multiplies y-component of traction jump
    % need factor of 4 instead of 2 because we are only using
    % half the terms in the quadrature
    P(j,vesicle.N+j,k) = P(j,vesicle.N+j,k) - ...
        sum(coeff)*4*pi/vesicle.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2
  end % j
end % k


end % pressDLmatrix




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S1,S2] = stressSLmatrix(o,vesicle)
% [S1,S2] = stressSLmatrix(vesicle), generates the matrix that returns
% the stress tensor of the single-layer potential applied to [1;0] (S1) 
% and to [0;1] (S2) given the traction jump.  Matricies have dimensions 
% (2*N,2*N,nv) where N is the number of points per curve and nv is the 
% number of curves in X.  Matrix is square since traction jump is 
% vector-valued and so is the stress

oc = curve;
[x,y] = oc.getXY(vesicle.X);
% Vesicle positions

sa = vesicle.sa; % Jacobian

S1 = zeros(2*vesicle.N,2*vesicle.N,vesicle.nv);
S2 = zeros(2*vesicle.N,2*veslcie.N,vesicle.nv);
% initialize double-layer potential to zero
for k=1:vesicle.nv  % Loop over curves
  index1 = (1:2:vesicle.N)'; % odd-indexed source points
  for j=2:2:vesicle.N % Loop over targets
    rho4 = ((x(j,k) - x(index1,k)).^2 + (y(j,k) - y(index1,k)).^2).^2;
    % distance squared

    coeff = (x(j,k) - x(index1,k)).^3.*vesicle.sa(index1,k)./rho4/pi;
    S1(j,index1,k) = 4*pi/vesicle.N*coeff;
    % part of first component of stress1 that multiplies first 
    % component of traction jump

    coeff = (x(j,k) - x(index1,k)).^2.*(y(j,k) - y(index1,k)).* ...
        vesicle.sa(index1,k)./rho4/pi;
    S1(j,index1+vesicle.N,k) = 4*pi/vesicle.N*coeff;
    % part of first component of stress1 that multiplies second
    % component of traction jump
    S1(j+vesicle.N,index1,k) = 4*pi/vesicle.N*coeff;
    % part of second component of stress1 that multiplies first 
    % component of traction jump
    S2(j,index1,k) = 4*pi/vesicle.N*coeff;
    % part of first component of stress2 that multiplies first 
    % component of traction jump

    coeff = (x(j,k) - x(index1,k)).*(y(j,k) - y(index1,k)).^2.* ...
        vesicle.sa(index1,k)./rho4/pi;
    S1(j+vesicle.N,index1+vesicle.N,k) = 4*pi/vesicle.N*coeff;
    % part of second component of stress1 that multiplies second 
    % component of traction jump
    S2(j,index1+vesicle.N,k) = 4*pi/vesicle.N*coeff;
    % part of first component of stress2 that multiplies second 
    % component of traction jump
    S2(j+vesicle.N,index1,k) = 4*pi/vesicle.N*coeff;
    % part of second component of stress2 that multiplies first 
    % component of traction jump

    coeff = (y(j,k) - y(index1,k)).^3.*vesicle.sa(index1,k)./rho4/pi;
    S2(j+vesicle.N,index1+vesicle.N,k) = 4*pi/vesicle.N*coeff;
    % part of second component of stress2 that multiplies second 
    % component of traction jump

    % need factors of 4 instead of 2 because we are only using
    % half the terms in the quadrature
  end % j

  index1 = (2:2:vesicle.N)'; % even-indexed source points
  for j=1:2:vesicle.N % Loop over targets
    rho4 = ((x(j,k) - x(index1,k)).^2 + (y(j,k) - y(index1,k)).^2).^2;
    % distance squared

    coeff = (x(j,k) - x(index1,k)).^3.*vesicle.sa(index1,k)./rho4/pi;
    S1(j,index1,k) = 4*pi/vesicle.N*coeff;
    % part of first component of stress1 that multiplies first 
    % component of traction jump

    coeff = (x(j,k) - x(index1,k)).^2.*(y(j,k) - y(index1,k)).* ...
        vesicle.sa(index1,k)./rho4/pi;
    S1(j,index1+vesicle.N,k) = 4*pi/vesicle.N*coeff;
    % part of first component of stress1 that multiplies second
    % component of traction jump
    S1(j+vesicle.N,index1,k) = 4*pi/vesicle.N*coeff;
    % part of second component of stress1 that multiplies first 
    % component of traction jump
    S2(j,index1,k) = 4*pi/vesicle.N*coeff;
    % part of first component of stress2 that multiplies first 
    % component of traction jump

    coeff = (x(j,k) - x(index1,k)).*(y(j,k) - y(index1,k)).^2.* ...
        vesicle.sa(index1,k)./rho4/pi;
    S1(j+vesicle.N,index1+vesicle.N,k) = 4*pi/vesicle.N*coeff;
    % part of second component of stress1 that multiplies second 
    % component of traction jump
    S2(j,index1+vesicle.N,k) = 4*pi/vesicle.N*coeff;
    % part of first component of stress2 that multiplies second 
    % component of traction jump
    S2(j+vesicle.N,index1,k) = 4*pi/vesicle.N*coeff;
    % part of second component of stress2 that multiplies first 
    % component of traction jump

    coeff = (y(j,k) - y(index1,k)).^3.*vesicle.sa(index1,k)./rho4/pi;
    S2(j+vesicle.N,index1+vesicle.N,k) = 4*pi/vesicle.N*coeff;
    % part of second component of stress2 that multiplies second 
    % component of traction jump

    % need factors of 4 instead of 2 because we are only using
    % half the terms in the quadrature
  end % j
end % k


end % stressSLmatrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D1,D2] = stressDLmatrix(o,vesicle)
% [D1,D2] = stressDLmatrix(vesicle), generates the matrix that returns
% the stress tensor due to the double-layer potential applied to [1;0] 
% (D1) and to [0;1] (D2) given the traction jump.  Matricies have 
% dimensions (2*N,2*N,nv) where N is the number of points per curve 
% and nv is the number of curves in X.  Matrix is square since traction 
% jump is vector-valued and so is the stress

normal = [vesicle.xt(vesicle.N+1:2*vesicle.N,:);...
          -vesicle.xt(1:vesicle.N,:)]; % Normal vector
oc = curve;
[x,y] = oc.getXY(vesicle.X);
% Vesicle positions

D1 = zeros(2*vesicle.N,2*vesicle.N,vesicle.nv);
D2 = zeros(2*vesicle.N,2*vesicle.N,vesicle.nv);
% initialize double-layer potential to zero
for k=1:vesicle.nv  % Loop over curves
  index1 = (1:2:vesicle.N)'; % odd-indexed source points
  for j=2:2:vesicle.N % Loop over targets
    rx = x(j,k) - x(index1,k);
    ry = y(j,k) - y(index1,k);
    rho2 = rx.^2 + ry.^2;
    % distance squared
    nx = normal(index1,k);
    ny = normal(index1+vesicle.N,k);
    % normal vector
    rdotn = rx.*nx + ry.*ny;
    % dot product of r with normal

    coeff = 1./rho2.*nx - ...
        8*rx.*rdotn./rho2.^3.*rx.^2 + ...
        rdotn./rho2.^2.*(2*rx) + ...
        1./rho2.^2.*(2*rx.^2.*nx);
    coeff = coeff/pi.*vesicle.sa(index1,k);
    D1(j,index1,k) = 4*pi/vesicle.N*coeff;
    % part of first component of stress1 that multiplies first 
    % component of traction jump
    D1(j,j,k) = D1(j,j,k) - sum(coeff)*4*pi/vesicle.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2

    coeff = 1./rho2.*ny - ...
        8*rx.*rdotn./rho2.^3.*rx.*ry + ...
        1./rho2.^2.*(2*rx.*ry.*nx);
    coeff = coeff/pi.*vesicle.sa(index1,k);
    D1(j,index1+vesicle.N,k) = 4*pi/vesicle.N*coeff;
    % part of first component of stress1 that multiplies second
    % component of traction jump
    D1(j,vesicle.N+j,k) = D1(j,vesicle.N+j,k) - ...
        sum(coeff)*4*pi/vesicle.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2

    coeff = -8*rx.*rdotn./rho2.^3.*rx.*ry + ...
        rdotn./rho2.^2.*ry + ...
        1./rho2.^2.*(rx.^2.*ny+rx.*ry.*nx);
    coeff = coeff/pi.*vesicle.sa(index1,k);
    D1(j+vesicle.N,index1,k) = 4*pi/vesicle.N*coeff;
    % part of second component of stress1 that multiplies first 
    % component of traction jump
    D1(vesicle.N+j,j,k) = D1(vesicle.N+j,j,k) - ...
        sum(coeff)*4*pi/vesicle.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2

    coeff = -8*rx.*rdotn./rho2.^3.*ry.^2 + ...
        rdotn./rho2.^2.*rx + ...
        1./rho2.^2.*(rx.*ry.*ny+ry.^2.*nx);
    coeff = coeff/pi.*vesicle.sa(index1,k);
    D1(j+vesicle.N,index1+vesicle.N,k) = 4*pi/vesicle.N*coeff;
    % part of second component of stress1 that multiplies second 
    % component of traction jump
    % Have built the stress tensor applied to [1;0] at half the points
    D1(vesicle.N+j,vesicle.N+j,k) = ...
        D1(vesicle.N+j,vesicle.vesicle.N+j,k) - ...
        sum(coeff)*4*pi/vesicle.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2


    coeff = -8*ry.*rdotn./rho2.^3.*rx.^2 + ...
        rdotn./rho2.^2.*ry + ...
        1./rho2.^2.*(rx.^2.*ny+rx.*ry.*nx);
    coeff = coeff/pi.*vesicle.sa(index1,k);
    D2(j,index1,k) = 4*pi/vesicle.N*coeff;
    % part of first component of stress1 that multiplies first 
    % component of traction jump
    D2(j,j,k) = D2(j,j,k) - sum(coeff)*4*pi/vesicle.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2

    coeff = -8*ry.*rdotn./rho2.^3.*rx.*ry + ...
        rdotn./rho2.^2.*rx + ...
        1./rho2.^2.*(rx.*ry.*ny+ry.^2.*nx);
    coeff = coeff/pi.*vesicle.sa(index1,k);
    D2(j,index1+vesicle.N,k) = 4*pi/vesicle.N*coeff;
    % part of first component of stress1 that multiplies second
    % component of traction jump
    D2(j,vesicle.N+j,k) = D2(j,vesicle.N+j,k) - ...
        sum(coeff)*4*pi/vesicle.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2

    coeff = 1./rho2.*nx - ...
        8*ry.*rdotn./rho2.^3.*rx.*ry + ...
        1./rho2.^2.*(2*rx.*ry.*ny);
    coeff = coeff/pi.*vesicle.sa(index1,k);
    D2(j+vesicle.N,index1,k) = 4*pi/vesicle.N*coeff;
    % part of second component of stress1 that multiplies first 
    % component of traction jump
    D2(vesicle.N+j,j,k) = D2(vesicle.N+j,j,k) - ...
        sum(coeff)*4*pi/vesicle.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2

    coeff = 1./rho2.*ny - ...
        8*ry.*rdotn./rho2.^3.*ry.^2 + ...
        rdotn./rho2.^2.*(2*ry) + ...
        1./rho2.^2.*(2*ry.^2.*ny);
    coeff = coeff/pi.*vesicle.sa(index1,k);
    D2(j+vesicle.N,index1+vesicle.N,k) = 4*pi/vesicle.N*coeff;
    % part of second component of stress1 that multiplies second 
    % component of traction jump
    D2(vesicle.N+j,vesicle.N+j,k) = ...
        D2(vesicle.N+j,vesicle.N+j,k) - ...
        sum(coeff)*4*pi/vesicle.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2
    % Have built the stress tensor applied to [0;1] at half the points


    % need factors of 4 instead of 2 because we are only using
    % half the terms in the quadrature
  end % j




  index1 = (2:2:vesicle.N)'; % even-indexed source points
  for j=1:2:vesicle.N % Loop over targets
    rx = x(j,k) - x(index1,k);
    ry = y(j,k) - y(index1,k);
    rho2 = rx.^2 + ry.^2;
    % distance squared
    nx = normal(index1,k);
    ny = normal(index1+vesicle.N,k);
    % normal vector
    rdotn = rx.*nx + ry.*ny;
    % dot product of r with normal

    coeff = 1./rho2.*nx - ...
        8*rx.*rdotn./rho2.^3.*rx.^2 + ...
        rdotn./rho2.^2.*(2*rx) + ...
        1./rho2.^2.*(2*rx.^2.*nx);
    coeff = coeff/pi.*vesicle.sa(index1,k);
    D1(j,index1,k) = 4*pi/vesicle.N*coeff;
    % part of first component of stress1 that multiplies first 
    % component of traction jump
    D1(j,j,k) = D1(j,j,k) - sum(coeff)*4*pi/vesicle.N;

    coeff = 1./rho2.*ny - ...
        8*rx.*rdotn./rho2.^3.*rx.*ry + ...
        1./rho2.^2.*(2*rx.*ry.*nx);
    coeff = coeff/pi.*vesicle.sa(index1,k);
    D1(j,index1+vesicle.N,k) = 4*pi/vesicle.N*coeff;
    % part of first component of stress1 that multiplies second
    % component of traction jump
    D1(j,vesicle.N+j,k) = D1(j,vesicle.N+j,k) - ...
        sum(coeff)*4*pi/vesicle.N;

    coeff = -8*rx.*rdotn./rho2.^3.*rx.*ry + ...
        rdotn./rho2.^2.*ry + ...
        1./rho2.^2.*(rx.^2.*ny+rx.*ry.*nx);
    coeff = coeff/pi.*vesicle.sa(index1,k);
    D1(j+vesicle.N,index1,k) = 4*pi/vesicle.N*coeff;
    % part of second component of stress1 that multiplies first 
    % component of traction jump
    D1(vesicle.N+j,j,k) = D1(vesicle.N+j,j,k) - ...
        sum(coeff)*4*pi/vesicle.N;

    coeff = -8*rx.*rdotn./rho2.^3.*ry.^2 + ...
        rdotn./rho2.^2.*rx + ...
        1./rho2.^2.*(rx.*ry.*ny+ry.^2.*nx);
    coeff = coeff/pi.*vesicle.sa(index1,k);
    D1(j+vesicle.N,index1+vesicle.N,k) = 4*pi/vesicle.N*coeff;
    % part of second component of stress1 that multiplies second 
    % component of traction jump
    D1(vesicle.N+j,vesicle.N+j,k) = ...
        D1(vesicle.N+j,vesicle.N+j,k) - ...
        sum(coeff)*4*pi/vesicle.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2
    % Have built the stress tensor applied to [0;1] at half the points
    % Have built the stress tensor applied to [1;0] at the other 
    % half of the points


    coeff = -8*ry.*rdotn./rho2.^3.*rx.^2 + ...
        rdotn./rho2.^2.*ry + ...
        1./rho2.^2.*(rx.^2.*ny+rx.*ry.*nx);
    coeff = coeff/pi.*vesicle.sa(index1,k);
    D2(j,index1,k) = 4*pi/vesicle.N*coeff;
    % part of first component of stress1 that multiplies first 
    % component of traction jump
    D2(j,j,k) = D2(j,j,k) - sum(coeff)*4*pi/vesicle.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2

    coeff = -8*ry.*rdotn./rho2.^3.*rx.*ry + ...
        rdotn./rho2.^2.*rx + ...
        1./rho2.^2.*(rx.*ry.*ny+ry.^2.*nx);
    coeff = coeff/pi.*vesicle.sa(index1,k);
    D2(j,index1+vesicle.N,k) = 4*pi/vesicle.N*coeff;
    % part of first component of stress1 that multiplies second
    % component of traction jump
    D2(j,vesicle.N+j,k) = D2(j,vesicle.N+j,k) - ...
        sum(coeff)*4*pi/vesicle.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2

    coeff = 1./rho2.*nx - ...
        8*ry.*rdotn./rho2.^3.*rx.*ry + ...
        1./rho2.^2.*(2*rx.*ry.*ny);
    coeff = coeff/pi.*vesicle.sa(index1,k);
    D2(j+vesicle.N,index1,k) = 4*pi/vesicle.N*coeff;
    % part of second component of stress1 that multiplies first 
    % component of traction jump
    D2(vesicle.N+j,j,k) = D2(vesicle.N+j,j,k) - ...
        sum(coeff)*4*pi/vesicle.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2

    coeff = 1./rho2.*ny - ...
        8*ry.*rdotn./rho2.^3.*ry.^2 + ...
        rdotn./rho2.^2.*(2*ry) + ...
        1./rho2.^2.*(2*ry.^2.*ny);
    coeff = coeff/pi.*vesicle.sa(index1,k);
    D2(j+vesicle.N,index1+vesicle.N,k) = 4*pi/vesicle.N*coeff;
    % part of second component of stress1 that multiplies second 
    % component of traction jump
    D2(vesicle.N+j,vesicle.N+j,k) = ...
        D2(vesicle.N+j,vesicle.N+j,k) - ...
        sum(coeff)*4*pi/vesicle.N;
    % need to subtract off the input at the target location
    % so that the singularity is only 1/r instead of 1/r^2
    % Have built the stress tensor applied to [0;1] at the other 
    % half of the points

    % need factors of 4 instead of 2 because we are only using
    % half the terms in the quadrature
  end % j
end % k

end % stressDLmatrix



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function uc = nearIntStockeval(o,vesicleSou,f, ...
%    kernel,vesicleTar,tEquals)
%% Near Singular Integration with stockeval    
%nvsou = vesicleSou.nv;
%nvtar = vesicleTar.nv;
%N = vesicleTar.N;
%uc = zeros(2*N,nvtar);
%if tEquals
%  t.x = zeros(N*(nvtar-1),1);
%else   
%  t.x = zeros(N*nvtar,1);
%  idxs = 1:nvtar;
%end
%for i = 1:nvsou
%  if tEquals
%    idxs = [(1:i-1) (i+1:nvsou)];
%  end
%  t.x(:) = vesicleTar.x(:,idxs);
%  SouCur.x = vesicleSou.x(:,i);
%  SouCur.xp = vesicleSou.xp(:,i);
%  SouCur.sp = vesicleSou.sa(:,i);
%  SouCur.nx = vesicleSou.nx(:,i);
%  SouCur.tang = vesicleSou.tang(:,i);
%  SouCur.t = vesicleSou.t;
%  SouCur.x = vesicleSou.x(:,i);
%  SouCur.cw = vesicleSou.cw(:,i);
%  SouCur.w = vesicleSou.w(:,i);
%  SouCur.a = vesicleSou.a(i);
%  fcur = f(1:N,i) + 1i*f(N+1:end,i);
%  uccur = kernel(t.x, SouCur, fcur, 'e');
%  for j = 1:numel(idxs)
%    uc(1:N,idxs(j)) = uc(1:N,idxs(j)) + ...
%       real( uccur((j-1)*N+1:j*N) );
%    uc(N+1:2*N,idxs(j)) = uc(N+1:2*N,idxs(j)) + ...
%       imag( uccur((j-1)*N+1:j*N) );
%  end
%end     
%
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LP = nearSingStokesSLP(o,vesicleSou,f,vesicleTar,...
    Gbarnett,tEqualS,ifmm)
% function LP = nearSingStokesSLP(o,vesicleSou,f,vesicleTar,...
%     Gbarnett,tEqualS,ifmm) 
% uses the near-singular laplace layer potentials to form the SLP for
% Stokes equation.  vesicleSou is an object of sources, f is the
% density function, vesicleTar is an object of targets, Gbarnett is a
% matrix that takes a density function and returns the boundary data of
% a complex-valued function that is analytic and whose real part is the
% desired single-layer potential, tEqualS is a flag stating if the
% targets is equal to the sources or not, and ifmm is a flag for using
% the fmm.  The output LP is the single-layer potential

Xsou = vesicleSou.X;
Nsou = size(Xsou,1)/2;
nvSou = size(Xsou,2);
Xtar = vesicleTar.X;
Ntar = size(Xtar,1)/2;
nvTar = size(Xtar,2);

uStokes = zeros(Ntar,nvTar);
vStokes = zeros(Ntar,nvTar);
% initialize two components of SLP

for ksou = 1:nvSou
  if tEqualS
    K = [(1:ksou-1) (ksou+1:nvTar)];
    % skip diagonal vesicle
  else % sources ~= targets
    K = (1:nvTar);
    % consider all vesicles
  end
  geom.x = Xsou(1:Nsou,ksou) + 1i*Xsou(Nsou+1:2*Nsou,ksou);
  geom.der = curve.arcDeriv(real(geom.x),1,ones(Nsou,1),...
        vesicleSou.IK(:,ksou)) + ...
      1i*curve.arcDeriv(imag(geom.x),1,ones(Nsou,1),...
        vesicleSou.IK(:,ksou));
  geom.sa = vesicleSou.sa(:,ksou);
  geom.nor = -1i*geom.der./geom.sa;
  geom.center = vesicleSou.center(1,ksou) + ...
             1i*vesicleSou.center(2,ksou);
  sigma1 = f(1:Nsou,ksou);
  sigma2 = f(Nsou+1:2*Nsou,ksou);
  side = 'e';

  zTar = Xtar(1:Ntar,K) + 1i*Xtar(Ntar+1:2*Ntar,K);
  zTar = zTar(:);

  u = zeros(size(zTar));
  v = zeros(size(zTar));
  [vSLP,vxSLP,vySLP] = o.laplaceSLPclose(geom,sigma1,zTar,side,...
      Gbarnett(:,:,ksou),ifmm);
  u = u + 0.5*vSLP;
  u = u - 0.5*real(zTar).*vxSLP;
  v = v - 0.5*real(zTar).*vySLP;
  [vSLP,vxSLP,vySLP] = o.laplaceSLPclose(geom,sigma2,zTar,side,...
      Gbarnett(:,:,ksou),ifmm);
  v = v + 0.5*vSLP;
  u = u - 0.5*imag(zTar).*vxSLP;
  v = v - 0.5*imag(zTar).*vySLP;
  souDOTsig = Xsou(1:Nsou,ksou).*sigma1 + ...
              Xsou(Nsou+1:2*Nsou,ksou).*sigma2;
  [~,vxSLP,vySLP] = o.laplaceSLPclose(geom,souDOTsig,zTar,side,...
      Gbarnett(:,:,ksou),ifmm);
  u = u + 0.5*vxSLP;
  v = v + 0.5*vySLP;

  for ktar = 1:numel(K)
    uStokes(:,K(ktar)) = uStokes(:,K(ktar)) + u((ktar-1)*Ntar+1:ktar*Ntar);
    vStokes(:,K(ktar)) = vStokes(:,K(ktar)) + v((ktar-1)*Ntar+1:ktar*Ntar);
  end
end

LP = [uStokes;vStokes];


end % nearSingStokesSLP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uSLP,uxSLP,uySLP] = laplaceSLPclose(o,...
    geom,sigma,zTar,side,Gbarnett,ifmm)
% function [uSLP,uxSLP,uySLP] = laplaceSLPclose(o,...
%    geom,sigma,zTar,side,Gbarnett)
% evalutes the potential and the gradient of the SLP using the
% barycentric ideas outlined in Barnett and Veerapaneni paper.  geom is
% structure with the geometry in complex form, sigma is the complex
% valued density function, zTar is the collection of target locations,
% side is either 'e' or 'i' which tells the function if the targets are
% inside or outside the geometry, Gbarnett is a precomputed matrix that
% takes a density function and maps it to the complex-valued boundary
% data of a holomorphic function whose real part is the desired layer
% potential, and ifmm is a flag to decide if the FMM is used or not

%d = repmat(geom.x,[1 N]) - repmat(geom.x.',[N 1]);
%x = exp(1i*theta);
%S = -log(d) + log(repmat(x,[1 N]) - repmat(x.',[N,1]));
%for k = 1:N
%  S(k,k) = log(1i*exp(1i*theta(k))./geom.der(k));
%end
%for i = 1:numel(S) - 1;
%  p = imag(S(i+1) - S(i));
%  S(i+1) = S(i+1) - 2*pi*1i*round(p/(2*pi));
%end
%
%m = 1:N/2-1;
%if side == 'i'
%  R = ifft([zeros(1,N/2) 1/N 1./(m(end:-1:1))]);
%else
%  R = ifft([0 1./m 1/N zeros(1,N/2-1)]);
%end
%R = R(:);
%R = toeplitz([R(1); R(end:-1:2)],R);
%
%S = S/N + R;
%vbd = S * (sigma.*geom.sa);
vbd = Gbarnett * (sigma.*geom.sa);
% compute the boundary data of the holomorphic function whose real part
% is the single-layer potential with density sigma
N = numel(sigma);
theta = (0:N-1)'*2*pi/N;

if side == 'e'
  sawlog = -1i*theta + log(geom.center - geom.x);
  for k = 1:numel(sawlog) - 1
    p = imag(sawlog(k+1) - sawlog(k));
    sawlog(k+1) = sawlog(k+1) - 2*1i*pi*round(p/(2*pi));
  end

  totCharge = sum(geom.sa.*sigma) * 1/N;
  vbd = vbd + totCharge * sawlog;

  cw = 1i*geom.nor.*geom.sa * 2*pi/N;
  vInf = sum(1./(geom.x - geom.center).*vbd.*cw)/(2*1i*pi);
  vbd = vbd - vInf;
end
% adjustment for exterior problem as outlined in Section 4.2.2 of
% Barnett and Veerapaneni

if ifmm
  [v,Dv] = o.cauchyCmplxFMM(zTar(:),geom,vbd,side);
else
  [v,Dv] = o.cauchyCmplxDirect(zTar(:),geom,vbd,side);
end
% compute the required cauchy integral

if side == 'e'
  vbd = vbd - totCharge*log(abs(geom.center - geom.x));
  v = v - totCharge*log(abs(geom.center - zTar(:))) + real(vInf);
  Dv = Dv + totCharge./(geom.center - zTar(:));
end
% correct the boundary data (for debugging only), and the LP and its
% gradient

uSLP = real(v);
uxSLP = real(Dv);
uySLP = -imag(Dv);

end % laplaceSLPclose

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uDLP,uxDLP,uyDLP,vbd] = laplaceDLPclose(o,...
    geom,sigma,Dsigma,zTar,side,ifmm)
% function [uDLP,uxDLP,uyDLP,vbd] = laplaceDLPclose(o,...
%     geom,sigma,Dsigma,zTar,side,ifmm)
% computes the double-layer potential of Laplace using method outlined
% by Barnett and Veerapaneni.  geom is the geometry in complex form,
% sigma and Dsigma are the density function and its derivative, zTar is
% the set of target points, and side is either 'e' or 'i' depending if
% the target points are interior or exterior to the geometry.  ifmm is
% a flag to signal the function to use or not use the FMM
% uDLP, uxDLP, uyDLP are the DLP and its gradient, and vbd is the
% complex-valued boundary data of uDLP

N = numel(geom.x);

vbd = zeros(N,1);
if ifmm
  disp('haven''t don''t this peice of the code yet')
else
  for k = 1:N
    ind = [(1:k-1) (k+1:N)];
    vbd(k) = sum((sigma(ind) - sigma(k))./...
        (geom.x(ind) - geom.x(k)).*geom.der(ind));
  end
end

vbd = 1i/2/pi * vbd * 2*pi/N + 1i/N * Dsigma;

if side == 'i'
  vbd = vbd - sigma;
end
% compute complex boundary value of the anlaytic continuation of the

% harmonic function we are seeking (the double-layer potential)

if ifmm
  [v,Dv] = o.cauchyCmplxFMM(zTar(:),geom,vbd,side);
else
  [v,Dv] = o.cauchyCmplxDirect(zTar(:),geom,vbd,side);
end
% compute the required cauchy integral

uDLP = real(v);
uxDLP = real(Dv);
uyDLP = -imag(Dv);

end % laplaceDLPclose

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v,Dv] = cauchyCmplxDirect(o,zTar,geom,vbd,side)
% [v,Dv] = cauchyCmplxDirect(zTar,geom,vbd,side) evaluates the cauchy
% integral theorem, and its complex derivative where vbd is the
% boundary data of a holomorphic function.  zTar is the set of target
% points, geom is a structure for the geometry, vbd is the boundary
% data, and side is 'i' or 'e' for interior or exterior, respectively

Nsou = numel(geom.x);
% number of source points
Ntar = numel(zTar);
% number of target points

v = zeros(Ntar,1);
% initialize function values

if side == 'i'
  for k = 1:Ntar
    I0 = sum(vbd./(geom.x - zTar(k)).*geom.der);
    J0 = sum(1./(geom.x - zTar(k)).*geom.der);
    v(k) = I0/J0;
  end
  % Cauchy integral for function values at interior points
else
  for k = 1:Ntar
    I0 = sum(vbd./(geom.x - zTar(k)).*geom.der);
    J0 = sum(1./(geom.x - geom.center)./(geom.x - zTar(k)).*geom.der);
    v(k) = 1/(zTar(k) - geom.center) * I0/J0;
  end
  % Cauchy integral for function values at exterior points
end

Dv = zeros(Ntar,1);
% initialize derivative values

if side == 'i'
  for k = 1:Ntar
    I0 = sum((vbd - v(k))./(geom.x - zTar(k)).^2.*geom.der);
    J0 = sum(1./(geom.x - zTar(k)).*geom.der);
    Dv(k) = I0/J0;
  end
  % Cauchy integral for derivative at interior points
else
  for k = 1:Ntar
    I0 = sum((vbd - v(k))./(geom.x - zTar(k)).^2.*geom.der);
    J0 = sum(1./(geom.x - geom.center)./(geom.x - zTar(k)).*geom.der);
    Dv(k) = 1/(zTar(k) - geom.center) * I0/J0;
  end
  % Cauchy integral for derivative at exterior points
end
% TODO: Have not included fix where term vbd - v(k) has a large error.
% This can be fixed as described in Barnett, Wu, Veerapaneni in equation
% (3.11)




end % cauchyCmplxDirect


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v,Dv] = cauchyCmplxFMM(o,zTar,geom,vbd,side)
% [v,Dv] = cauchyCmplxFMM(zTar,geom,vbd,side) evaluates the cauchy
% integral theorem, and its complex derivative where vbd is the
% boundary data of a holomorphic function.  zTar is the set of target
% points, geom is a structure for the geometry, vbd is the boundary
% data, and side is 'i' or 'e' for interior or exterior, respectively

Nsou = numel(geom.x);
% number of source points
Ntar = numel(zTar);
% number of target points

x = [real(geom.x); real(zTar)];
y = [imag(geom.x); imag(zTar)];
% stack source points followed by target points

I0 = zeros(Ntar,1);
DI0 = zeros(Ntar,1);
J0 = zeros(Ntar,1);
DJ0 = zeros(Ntar,1);

sig = [vbd.*geom.der;zeros(Ntar,1)];
[rI0,iI0,rDI0,iDI0] = laplaceLPfmm(real(sig),imag(sig),x,y);
I0 = I0 + rI0(Nsou+1:end) + 1i*iI0(Nsou+1:end);
DI0 = DI0 + rDI0(Nsou+1:end) + 1i*iDI0(Nsou+1:end);

sig = [geom.der;zeros(Ntar,1)];
[rJ0,iJ0,rDJ0,iDJ0] = laplaceLPfmm(real(sig),imag(sig),x,y);
DJ0save = rDJ0(Nsou+1:end) + 1i*iDJ0(Nsou+1:end);
if side == 'i';
  J0 = J0 + rJ0(Nsou+1:end) + 1i*iJ0(Nsou+1:end);
  DJ0 = DJ0 + rJ0(Nsou+1:end) + 1i*iJ0(Nsou+1:end);
end

if side == 'e'
  sig = [geom.der./(geom.x - geom.center);zeros(Ntar,1)];
  [rJ0,iJ0,~,~] = laplaceLPfmm(real(sig),imag(sig),x,y);
  J0 = J0 + rJ0(Nsou+1:end) + 1i*iJ0(Nsou+1:end);
  DJ0 = DJ0 + rJ0(Nsou+1:end) + 1i*iJ0(Nsou+1:end);
end

v = I0./J0;
if side == 'e'
  v = v./(zTar - geom.center);
end

Dv = (DI0 - v.*DJ0save)./DJ0;
if side == 'e'
  Dv = Dv./(zTar - geom.center);
end
% TODO: Have not included fix where term vbd - v(k) has a large error.
% This can be fixed as described in Barnett, Wu, Veerapaneni in equation
% (3.11)


end % cauchyCmplxFMM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G = laplaceSLcomplexMatrix(o,vesicle)
% G = laplaceSLcomplexMatrix(vesicle) generates complex boudnary
% condition which is the analytic continuation of the single-layer
% potential of a given input density function.  vesicle is a data
% structure defined as in the capsules class G is (N,N,nv) array where
% N is the number of points per curve and nv is the number of curves in
% X
% NOTE THAT THIS MATRIX CORRESPONDS TO THE LIMITING TERM AS A POINT
% APPROACHES THE BOUNDARY FROM THE EXTERIOR.  SO FAR, THIS IS THE ONLY
% CASE WE NEED TO CONSIDER FOR THE SLP SINCE VESICLES AND SOLID WALLS
% ARE ALWAYS EXTERIOR TO THE OTHER VESICLES

oc = curve;
[x,y] = oc.getXY(vesicle.X); % Vesicle positions

G = zeros(o.N,o.N,vesicle.nv);
% initalize single-layer potential to zero
theta = (0:o.N-1)'*2*pi/o.N;

for kves=1:vesicle.nv  % Loop over curves
  z = x(:,kves) + 1i*y(:,kves);
  zDer = curve.arcDeriv(real(z),1,ones(o.N,1),vesicle.IK(:,kves)) + ...
      1i*curve.arcDeriv(imag(z),1,ones(o.N,1),vesicle.IK(:,kves));

  zexp = exp(1i*theta);
  d = repmat(z,[1,o.N]) - repmat(z.',[o.N,1]);

  S = -log(d) + log(...
      repmat(zexp,[1 o.N]) - repmat(zexp.',[o.N 1]));
  for k = 1:o.N
    S(k,k) = log(1i*zexp(k)./zDer(k));
  end

  p = S(:); p = imag(p(2:end) - p(1))';
  S(2:end) = S(2:end) - 2*pi*1i*round(p/(2*pi));
  % vectorized method of removing the branch cut

%  for k = 1:o.N^2 - 1
%    p = imag(S(k+1) - S(k));
%    S(k+1) = S(k+1) - 2*pi*1i*round(p/(2*pi));
%  end
% old way of removing the branch cut.  The for loop is over N^{2}
% entries, so this is a killer
  m = 1:o.N/2-1;
  R = ifft([0 1./m 1/o.N zeros(1,o.N/2-1)]);
  % This vector R is different if we are taking the limit from the
  % interior rather than the exterior
  R = R(:);
  R = toeplitz([R(1); R(end:-1:2)],R);

  S = S/o.N + R;
  G(:,:,kves) = S;

end % k

end % laplaceSLcomplexMatrix


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LP = nearSingInt(o,vesicleSou,f,selfMat,...
    NearStruct,kernel,vesicleTar,tEqualS,trash)
% LP = nearSingInt(vesicle,f,selfMat,...
% zone,dist,nearest,icp,argnear,vesicleTar,tequalS) computes a 
% layer potential due to f at all points in vesicleTar.X.  If
% tEqualS==true, then the vesicleTar == vesicleSou and the
% self-vesicle interaction is skipped.
% selfMat is the diagonal of the potential needed to compute
% the layer potential of each vesicle indepenedent of all others
% if using near-singular integration
% zone,dist,nearest,icp,argnear are required by near-singular
% integration (they keep everything sorted and precomputed)
% Everything is in the 2*N x nv format
% Can pass a final argument if desired so that plots of the 
% near-singular integration algorithm are displayed

dist = NearStruct.dist;
zone = NearStruct.zone;
nearest = NearStruct.nearest;
icp = NearStruct.icp;
argnear = NearStruct.argnear;

Xsou = vesicleSou.X; % source positions
Nsou = size(Xsou,1)/2; % number of source points
nvSou = size(Xsou,2); % number of source 'vesicles'
Xtar = vesicleTar.X; % target positions
Ntar = size(Xtar,1)/2; % number of target points
nvTar = size(Xtar,2); % number of target 'vesicles'

h = vesicleSou.length/Nsou; % arclength term

Nup = Nsou*2^ceil(1/2*log2(Nsou));
Nup = Nsou*ceil(sqrt(Nsou));
% upsample to N^(3/2).  

vself = selfMat(f);
%vself = zeros(2*Nsou,nvSou);
%for k=1:nvSou
%  vself(:,k) = selfMat(:,:,k)*f(:,k);
%end
% Compute velocity due to each vesicle independent of others.
% This is needed when using near-singular integration since
% we require a value for the layer-potential on the vesicle of 
% sources 

Xup = [interpft(Xsou(1:Nsou,:),Nup);...
       interpft(Xsou(Nsou+1:2*Nsou,:),Nup)];
fup = [interpft(f(1:Nsou,:),Nup);...
       interpft(f(Nsou+1:2*Nsou,:),Nup)];
% upsample positions, traction jump

vesicleUp = capsules(Xup,[],[],vesicleSou.kappa,vesicleSou.viscCont);
% Build an object with the upsampled vesicle

interpOrder = size(o.interpMat,1);
% lagrange interpolation order
p = ceil((interpOrder+1)/2);
% want half of the lagrange interpolation points to the left of the
% closest point and the other half to the right
vel = zeros(2*Ntar,nvTar,nvSou);
% allocate space for storing velocity at intermediate points
% needed by near-singular integration

if tEqualS % sources == targets
  if nvSou > 1
    if strfind(char(kernel),'fmm');
      farField = kernel(vesicleUp,fup);
      farField = farField(1:Nup/Ntar:end,:);
      % evaluate layer potential at all targets except ignore the
      % diagonal term
    else
      for k = 1:nvSou
        K = [(1:k-1) (k+1:nvSou)];
        [~,farField(:,k)] = kernel(vesicleUp,fup,Xtar(:,k),K);
      end
%      % This is a huge savings if we are using a direct method rather
%      % than the fmm to evaluate the layer potential.  The speedup is
%      % more than N^{1/2}, where N is the resolution of the vesicles
%      that % we are computing on
    end
  else
    farField = zeros(2*Ntar,nvTar);
  end

else % sources ~= targets
  [~,farField] = kernel(vesicleUp,fup,Xtar,1:nvSou);
  % evaluate layer potential due to all 'vesicles' at all points
  % in Xtar;
end
% Use upsampled trapezoid rule to compute layer potential

nearField = zeros(2*Ntar,nvTar);
% Initialize potential at near points to zero

beta = 1.1;
% small buffer to make sure Lagrange interpolation points are
% not in the near zone
for k1 = 1:nvSou
  if tEqualS % sources == targets
    K = [(1:k1-1) (k1+1:nvTar)];
    % skip diagonal vesicle
  else % sources ~= targets
    K = (1:nvTar);
    % consider all vesicles
  end
  for k2 = K
    J = find(zone{k1}(:,k2) == 1);
    % set of points on vesicle k2 close to vesicle k1
    if (numel(J) ~= 0)
      indcp = icp{k1}(J,k2);
      % closest point on vesicle k1 to each point on vesicle k2 that is
      % close to k1
      for j = 1:numel(J)
        pn = mod((indcp(j)-p+1:indcp(j)-p+interpOrder)' - 1,Nsou) + 1;
        % index of points to the left and right of the closest point
        v = filter(1,[1 -full(argnear{k1}(J(j),k2))],...
          o.interpMat*vself(pn,k1));
        vel(J(j),k2,k1) = v(end);
        % x-component of the velocity at the closest point
        v = filter(1,[1 -full(argnear{k1}(J(j),k2))],...
          o.interpMat*vself(pn+Nsou,k1));
        vel(J(j)+Ntar,k2,k1) = v(end);
        % y-component of the velocity at the closest point
      end
%     compute values of velocity at required intermediate points
%     using local interpolant

      [~,potTar] = kernel(vesicleUp,fup,...
          [Xtar(J,k2);Xtar(J+Ntar,k2)],k1);
      % Need to subtract off contribution due to vesicle k1 since
      % its layer potential will be evaulted using Lagrange
      % interpolant of nearby points
      nearField(J,k2) =  - potTar(1:numel(J));
      nearField(J+Ntar,k2) =  - potTar(numel(J)+1:end);

      XLag = zeros(2*numel(J),interpOrder - 1);
      for i = 1:numel(J)
        nx = (Xtar(J(i),k2) - nearest{k1}(J(i),k2))/...
            dist{k1}(J(i),k2);
        ny = (Xtar(J(i)+Ntar,k2) - nearest{k1}(J(i)+Ntar,k2))/...
            dist{k1}(J(i),k2);

        XLag(i,:) = nearest{k1}(J(i),k2) + ...
            beta*h*nx*(1:interpOrder-1);
        XLag(i+numel(J),:) = nearest{k1}(J(i)+Ntar,k2) + ...
            beta*h*ny*(1:interpOrder-1);
        % Lagrange interpolation points coming off of vesicle k1
        % All points are behind Xtar(J(i),k2) and are sufficiently
        % far from vesicle k1 so that the Nup-trapezoid rule gives
        % sufficient accuracy
      end
      [~,lagrangePts] = kernel(vesicleUp,fup,XLag,k1);
      % evaluate velocity at the lagrange interpolation points
      x = XLag(1,:);
      y = XLag(2,:);

      for i = 1:numel(J)
        Px = o.interpMat*[vel(J(i),k2,k1) ...
            lagrangePts(i,:)]';
        Py = o.interpMat*[vel(J(i)+Ntar,k2,k1) ...
            lagrangePts(i+numel(J),:)]';
        % Build polynomial interpolant along the one-dimensional
        % points coming out of the vesicle
        dscaled = full(dist{k1}(J(i),k2)/(beta*h*(interpOrder-1)));
        % Point where interpolant needs to be evaluated

        v = filter(1,[1 -dscaled],Px);
        v1 = v;
        nearField(J(i),k2) = nearField(J(i),k2) + ...
            v(end);
        v = filter(1,[1 -dscaled],Py);
        v2 = v;
        nearField(J(i)+Ntar,k2) = nearField(J(i)+Ntar,k2) + ...
            v(end);
        % Assign higher-order results coming from Lagrange 
        % integration to velocity at near point.  Filter is faster
        % than polyval
        
        if (nargin == 9)
          figure(2); clf; hold on;
          plot(Xsou(1:Nsou,:),Xsou(Nsou+1:end,:),'r.')
          plot(Xtar(1:Ntar,:),Xtar(Ntar+1:end,:),'k.')
          plot(Xtar(J,k2),Xtar(Ntar+J,k2),'b.')
          plot(XLag(1:numel(J),:),XLag(numel(J)+1:end,:),'kx')
          plot(XLag(i,:),XLag(numel(J)+i,:),'gx')
          axis equal

          figure(1); clf; hold on
          plot((0:interpOrder-1)*beta*h,...
              [vel(J(i),k2,k1) lagrangePts(i,:)],'g-o')
          plot((0:interpOrder-1)*beta*h,...
              [vel(J(i)+Ntar,k2,k1) lagrangePts(i+numel(J),:)],'r--o')
          pause
        end
        % DEBUG: PASS IN A DUMMY VARIABLE INTO THIS ROUTINE AND THEN
        % YOU CAN SEE THE INTERPOLATION POINTS AND CHECK THE SMOOTHNESS
        % OF THE INTERPOLANT

      end % i
    end % numel(J) ~= 0
    % Evaluate layer potential at Lagrange interpolation
    % points if there are any
  end % k2
end % k1

%if tEqualS % vesicleSou == vesicleTar
%  LP = farField(1:Nup/Ntar:end,:) + nearField;
%else % vesicleSou ~= vesicleTar
%  LP = farField + nearField;
%end
LP = farField + nearField;
% Add kernel due to far points and near points.  Far points were
% upsampled if source==vesicle so need to truncate here.  We are 
% only using Ntar target points.  Note that it is only the sources 
% that were upsampled


end % nearSingInt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SLP = exactStokesSLdiag(o,vesicle,G,f)
% SLP = exactStokesSLdiag(vesicle,G(f) computes the diagonal term of
% the single-layer potential due to f around vesicle.  Source and
% target points are the same.  This uses Alpert's quadrature formula.
% The goal is to get rid of forming matricies, which is why this
% function is required

SLP = zeros(size(f));
for k = 1:vesicle.nv
  SLP(:,k) = G(:,:,k) * f(:,k);
end


end % exactStokesSLdiag




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DLP = exactStokesDLdiag(o,vesicle,D,f)
% DLP = exactStokesDLdiag(vesicle,f,K) computes the diagonal term of
% the double-layer potential due to f around all vesicles stored in K.
% If K is empty all interaction due to all vesicles is computed  Source
% and target points are the same.  This uses trapezoid rule with the
% curvature at the diagonal.  The goal is to get rid of forming
% matricies, which is why this function is required


DLP = zeros(2*vesicle.N,vesicle.nv);
for k = 1:vesicle.nv
  DLP(:,k) = D(:,:,k) * f(:,k);
end

end % exactStokesDLdiag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N0 = exactStokesN0diag(o,vesicle,f)
% DLP = exactStokesN0diag(vesicle,f) computes the diagonal term of the
% modification of the double-layer potential due to f around outermost
% vesicle.  Source and target points are the same.  This uses trapezoid
% rule

oc = curve;
den = f.*[vesicle.sa(:,1);vesicle.sa(:,1)];
% multiply by the Jacobian
[x,y] = oc.getXY(vesicle.X(:,1));
[denx,deny] = oc.getXY(den);
% position
nx = vesicle.xt(vesicle.N+1:2*vesicle.N,1);
ny = -vesicle.xt(1:vesicle.N:vesicle.N,1);
% normal vector

N0 = zeros(2*vesicle.N,1);
fdotn = sum(denx.*nx + deny.*ny);
N0(1:vesicle.N) = 2*pi/vesicle.N*fdotn*nx;
N0(1+vesicle.N:end) = 2*pi/vesicle.N*fdotn*ny;


end % exactStokesN0diag


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SLP = exactLaplaceSLdiag(o,vesicle,G,f)
% SLP = exactLaplaceSLdiag(vesicle,G,f) computes the diagonal term of
% the single-layer potential due to f around vesicle.  Source and
% target points are the same.  For now, we just pass the matrix of
% SLP's and loop

SLP = zeros(size(f));
for k = 1:vesicle.nv
  SLP(:,k) = G(:,:,k) * f(:,k);
end
%SLP = [SLP;zeros(vesicle.N,vesicle.nv)];


end % exactLaplaceSLdiag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokesSLP,stokesSLPtar] = ...
    exactStokesSL(o,vesicle,f,Xtar,K1)
% [stokesSLP,stokesSLPtar] = exactStokesSL(vesicle,f,Xtar,K1)
% computes the single-layer potential due to f around all vesicles 
% except itself.  Also can pass a set of target points Xtar and a 
% collection of vesicles K1 and the single-layer potential due to
% vesicles in K1 will be evaluated at Xtar.
% Everything but Xtar is in the 2*N x nv format
% Xtar is in the 2*Ntar x ncol format

if nargin == 5
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  stokesSLPtar = zeros(2*Ntar,ncol);
else
  K1 = [];
  stokesSLPtar = [];
  ncol = 0;
  % if nargin ~= 5, the user does not need  the velocity at arbitrary
  % points
end

den = f.*[vesicle.sa;vesicle.sa]*2*pi/vesicle.N;
% multiply by arclength term

for k = 1:ncol % loop over columns of target points
  for j = 1:Ntar % loop over rows of target points
    diffxy = [Xtar(j,k) - vesicle.X(1:vesicle.N,K1) ; ...
        Xtar(j+Ntar,k) - vesicle.X(vesicle.N+1:2*vesicle.N,K1)];
    dis2 = diffxy(1:vesicle.N,:).^2 + ...
        diffxy(vesicle.N+1:2*vesicle.N,:).^2;
    % distance squared and difference of source and target location

    coeff = 0.5*log(dis2);
    % first part of single-layer potential for Stokes
    stokesSLPtar(j,k) = -sum(sum(...
        coeff.*den(1:vesicle.N,K1)));
    stokesSLPtar(j+Ntar,k) = -sum(sum(...
        coeff.*den(vesicle.N+1:2*vesicle.N,K1)));
    % log part of stokes single-layer potential

    coeff = (diffxy(1:vesicle.N,:).*den(1:vesicle.N,K1) + ...
        diffxy(vesicle.N+1:2*vesicle.N,:).*...
        den(vesicle.N+1:2*vesicle.N,K1))./dis2;
    % second part of single-layer potential for Stokes

    stokesSLPtar(j,k) = stokesSLPtar(j,k) + sum(sum(...
        coeff.*diffxy(1:vesicle.N,:)));
    stokesSLPtar(j+Ntar,k) = stokesSLPtar(j+Ntar,k) + sum(sum(...
        coeff.*diffxy(vesicle.N+1:2*vesicle.N,:)));
    % r \otimes r term of the stokes single-layer potential
  end % j

end % k
% Evaluate single-layer potential at arbitrary target points
stokesSLPtar = 1/(4*pi)*stokesSLPtar;
% 1/4/pi is the coefficient in front of the single-layer potential

stokesSLP = zeros(2*vesicle.N,vesicle.nv); % Initialize to zero
if (nargin == 3 && vesicle.nv > 1)
  for k = 1:vesicle.nv % vesicle of targets
    K = [(1:k-1) (k+1:vesicle.nv)];
    % Loop over all vesicles except k
    for j=1:vesicle.N
      dis2 = (vesicle.X(j,k) - vesicle.X(1:vesicle.N,K)).^2 + ...
          (vesicle.X(j+vesicle.N,k) - ...
            vesicle.X(vesicle.N+1:2*vesicle.N,K)).^2;
      diffxy = [vesicle.X(j,k) - vesicle.X(1:vesicle.N,K) ; ...
          vesicle.X(j+vesicle.N,k) - ...
            vesicle.X(vesicle.N+1:2*vesicle.N,K)];
      % distance squared and difference of source and target location

      coeff = 0.5*log(dis2);
      % first part of single-layer potential for Stokes

      val = coeff.*den(1:vesicle.N,K);
      stokesSLP(j,k) = -sum(val(:));
      val = coeff.*den(vesicle.N+1:2*vesicle.N,K);
      stokesSLP(j+vesicle.N,k) = -sum(val(:));
      % logarithm terms in the single-layer potential

      coeff = (diffxy(1:vesicle.N,:).*den(1:vesicle.N,K) + ...
          diffxy(vesicle.N+1:2*vesicle.N,:).*...
          den(vesicle.N+1:2*vesicle.N,K))./dis2;
      % second part of single-layer potential for Stokes

      val = coeff.*diffxy(1:vesicle.N,:);
      stokesSLP(j,k) = stokesSLP(j,k) + sum(val(:));
      val = coeff.*diffxy(vesicle.N+1:2*vesicle.N,:);
      stokesSLP(j+vesicle.N,k) = stokesSLP(j+vesicle.N,k) + ...
          sum(val(:));
      % r \otimes r term of the single-layer potential

    end % j
  end % k
  % Evaluate single-layer potential at vesicles but oneself
  stokesSLP = 1/(4*pi)*stokesSLP;
  % 1/4/pi is the coefficient in front of the single-layer potential
end % nargin == 4

end % exactStokesSL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokesDLP,stokesDLPtar] = ...
    exactStokesDL(o,vesicle,f,Xtar,K1)
% [stokesDLP,stokesDLPtar] = exactStokesDL(vesicle,f,Xtar,K1) computes
% the double-layer potential due to f around all vesicles except
% itself.  Also can pass a set of target points Xtar and a collection
% of vesicles K1 and the double-layer potential due to vesicles in K1
% will be evaluated at Xtar.  Everything but Xtar is in the 2*N x nv
% format Xtar is in the 2*Ntar x ncol format

normal = [vesicle.xt(vesicle.N+1:2*vesicle.N,:);...
         -vesicle.xt(1:vesicle.N,:)]; 
% Normal vector

if nargin == 5
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  stokesDLPtar = zeros(2*Ntar,ncol);
else
  K1 = [];
  stokesDLPtar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 5, the user does not need the velocity at arbitrary
  % points
end

den = f.*[vesicle.sa;vesicle.sa]*2*pi/vesicle.N;

for k = 1:ncol % loop over columns of target points
  for j = 1:Ntar % loop over rows of target points
    diffxy = [Xtar(j,k) - vesicle.X(1:vesicle.N,K1) ; ...
        Xtar(j+Ntar,k) - vesicle.X(vesicle.N+1:2*vesicle.N,K1)];
    dis2 = diffxy(1:vesicle.N,:).^2 + ...
        diffxy(vesicle.N+1:2*vesicle.N,:).^2;
    % difference of source and target location and distance squared

    rdotnTIMESrdotf = ...
      (diffxy(1:vesicle.N,:).*normal(1:vesicle.N,K1) + ...
      diffxy(vesicle.N+1:2*vesicle.N,:).*...
      normal(vesicle.N+1:2*vesicle.N,K1))./dis2.^2.* ...
      (diffxy(1:vesicle.N,:).*den(1:vesicle.N,K1) + ...
      diffxy(vesicle.N+1:2*vesicle.N,:).*...
      den(vesicle.N+1:2*vesicle.N,K1));
    % \frac{(r \dot n)(r \dot density)}{\rho^{4}} term
    rdotnTIMESrdotf = rdotnTIMESrdotf * diag(1-vesicle.viscCont(K1));


    stokesDLPtar(j,k) = stokesDLPtar(j,k) + ...
        sum(sum(rdotnTIMESrdotf.*diffxy(1:vesicle.N,:)));
    stokesDLPtar(j+Ntar,k) = stokesDLPtar(j+Ntar,k) + ...
        sum(sum(rdotnTIMESrdotf.*diffxy(vesicle.N+1:2*vesicle.N,:)));
    % r \otimes r term of the single-layer potential
  end % j
end % k
%stokesDLPtar = (1-vesicle.viscCont(K1)) * stokesDLPtar/pi;
stokesDLPtar = stokesDLPtar/pi;
% double-layer potential due to vesicles indexed over K1 
% evaluated at arbitrary points

stokesDLP = zeros(2*vesicle.N,vesicle.nv);
if (nargin == 3 && vesicle.nv > 1)
  oc = curve;
  for k = 1:vesicle.nv
    K = [(1:k-1) (k+1:vesicle.nv)];
    [x,y] = oc.getXY(vesicle.X(:,K));
    [nx,ny] = oc.getXY(normal(:,K));
    [denx,deny] = oc.getXY(den(:,K));
    for j=1:vesicle.N
      diffxy = [vesicle.X(j,k) - x ; vesicle.X(j+vesicle.N,k) - y];
      dis2 = diffxy(1:vesicle.N,:).^2 + ...
          diffxy(vesicle.N+1:2*vesicle.N,:).^2;
      % difference of source and target location and distance squared

      rdotfTIMESrdotn = ...
        (diffxy(1:vesicle.N,:).*nx + ...
        diffxy(vesicle.N+1:2*vesicle.N,:).*ny)./dis2.^2 .* ...
        (diffxy(1:vesicle.N,:).*denx + ...
        diffxy(vesicle.N+1:2*vesicle.N,:).*deny);
      % \frac{(r \dot n)(r \dot density)}{\rho^{4}} term

%      diffxy = vesicle.cachedVars.diffxy{j,k};
%      rdotfTIMESrdotn = vesicle.cachedVars.rdotnOVERdis4{j,k}.* ...
%          (diffxy(1:vesicle.N,:).*denx + ...
%          diffxy(vesicle.N+1:2*vesicle.N,:).*deny);

      stokesDLP(j,k) = stokesDLP(j,k) + ...
          sum(sum(rdotfTIMESrdotn.*diffxy(1:vesicle.N,:)));
      stokesDLP(j+vesicle.N,k) = stokesDLP(j+vesicle.N,k) + ...
          sum(sum(rdotfTIMESrdotn.*diffxy(vesicle.N+1:2*vesicle.N,:)));
      % double-layer potential for Stokes
    end
  end

  stokesDLP = stokesDLP*(eye(vesicle.nv) - diag(vesicle.viscCont))/pi;
  % 1/pi is the coefficient in front of the double-layer potential
end
% double-layer potential due to all vesicles except oneself



end % exactStokesDL



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokesSLP,stokesSLPtar] = ...
    exactStokesSLfmm(o,vesicle,f,Xtar,K)
% [stokesSLP,stokeSLPtar] = exactStokesSLfmm(vesicle,f,Xtar,K) uses 
% the FMM to compute the single-layer potential due to all vesicles
% except itself vesicle is a class of object capsules and f is the 
% density function NOT scaled by arclength term.  Xtar is a set of 
% points where the single-layer potential due to all vesicles in index 
% set K needs to be evaulated
global fmms

fmms = fmms + 1;
% count the total number of calls to fmm

N = vesicle.N; % number of points per vesicle
nv = vesicle.nv; % number of vesicles
X = vesicle.X; % vesicle positions
oc = curve;
[x,y] = oc.getXY(X); % seperate x and y coordinates

den = f.*[vesicle.sa;vesicle.sa]*2*pi/N;

if (nargin == 5)
  stokesSLP = [];
else
  [f1,f2] = oc.getXY(den);
  % need to multiply by arclength term.  Seperate it into
  % x and y coordinate

  [u,v] = stokesSLPfmm(f1(:),f2(:),x(:),y(:));
  stokesSLP = zeros(2*N,nv); % initialize
  for k = 1:nv
    is = (k-1)*N+1;
    ie = k*N;
    stokesSLP(1:N,k) = u(is:ie);
    stokesSLP(N+1:2*N,k) = v(is:ie);
  end
  % Wrap the output of the FMM into the usual 
  % [[x1;y1] [x2;y2] ...] format


  for k = 1:nv
    [u,v] = stokesSLPfmm(f1(:,k),f2(:,k),x(:,k),y(:,k));
    stokesSLP(:,k) = stokesSLP(:,k) - [u;v];
  end
  % Subtract potential due to each vesicle on its own.  Nothing
  % to change here for potential at Xtar
end

if nargin == 3
  stokesSLPtar = [];
else
  [x,y] = oc.getXY(X(:,K)); 
  % seperate x and y coordinates at vesicles indexed by K
  [Ntar,ncol] = size(Xtar);
  Ntar = Ntar/2;
  x2 = Xtar(1:Ntar,:);
  x = [x(:);x2(:)];
  y2 = Xtar(Ntar+1:2*Ntar,:);
  y = [y(:);y2(:)];
  % Stack the x and y coordinates of the target points
  [f1,f2] = oc.getXY(den(:,K));
  % seperate x and y coordinates at vesicles indexed by K
  f1 = [f1(:);zeros(Ntar*ncol,1)];
  f2 = [f2(:);zeros(Ntar*ncol,1)];
  % pad density function with zeros so that Xtar doesn't
  % affect the single-layer potential

  [u,v] = stokesSLPfmm(f1,f2,x,y);
  stokesSLPtar = zeros(2*Ntar,ncol); % initialize
  for k = 1:ncol
    is = N*numel(K) + (k-1)*Ntar+1;
    ie = is + Ntar - 1;
    stokesSLPtar(1:Ntar,k) = u(is:ie);
    stokesSLPtar(Ntar+1:2*Ntar,k) = v(is:ie);
  end
  % Wrap the output of the FMM in the usual format
  % for the target points
end

end % exactStokesSLfmm




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokesDLP,stokesDLPtar] = ...
    exactStokesDLfmm(o,vesicle,f,Xtar,K)
% [stokesDLP,stokeDLPtar] = exactStokesDLfmm(vesicle,f,Xtar,K) uses 
% the FMM to compute the double-layer potential due to all vesicles
% except itself vesicle is a class of object capsules and f is the 
% density function NOT scaled by arclength term.  Xtar is a set of 
% points where the single-layer potential due to all vesicles in index 
% set K needs to be evaulated
global fmms

fmms = fmms + 1;
% count the total number of calls to fmm

N = vesicle.N; % number of points per vesicle
nv = vesicle.nv; % number of vesicles
X = vesicle.X; % vesicle positions
oc = curve;
[x,y] = oc.getXY(X); % seperate x and y coordinates
[nx,ny] = oc.getXY(vesicle.normal);
% seperate the x and y coordinates of the normal vector

den = f.*[vesicle.sa;vesicle.sa]*2*pi/N;

if (nargin == 5)
  stokesDLP = [];
else
  [f1,f2] = oc.getXY(den);
  % need to multiply by arclength term.  Seperate it into
  % x and y coordinate

  [u,v] = stokesDLPfmm(f1(:),f2(:),x(:),y(:),nx(:),ny(:));

  stokesDLP = zeros(2*N,nv); % initialize
  for k = 1:nv
    is = (k-1)*N+1;
    ie = k*N;
    stokesDLP(1:N,k) = u(is:ie);
    stokesDLP(N+1:2*N,k) = v(is:ie);
  end
  % Wrap the output of the FMM into the usual 
  % [[x1;y1] [x2;y2] ...] format


  for k = 1:nv
    [u,v] = stokesDLPfmm(f1(:,k),f2(:,k),x(:,k),y(:,k),...
        nx(:,k),ny(:,k));
    stokesDLP(:,k) = stokesDLP(:,k) - [u;v];
  end
  % Subtract potential due to each vesicle on its own.  Nothing
  % to change here for potential at Xtar
end

if nargin == 3
  stokesDLPtar = [];
else
  [x,y] = oc.getXY(X(:,K)); 
  % seperate x and y coordinates at vesicles indexed by K
  [nx,ny] = oc.getXY(vesicle.normal(:,K));
  [Ntar,ncol] = size(Xtar);
  Ntar = Ntar/2;
  x2 = Xtar(1:Ntar,:);
  x = [x(:);x2(:)];
  y2 = Xtar(Ntar+1:2*Ntar,:);
  y = [y(:);y2(:)];
  % Stack the x and y coordinates of the target points
  [f1,f2] = oc.getXY(den(:,K));
  % seperate x and y coordinates at vesicles indexed by K
  f1 = [f1(:);zeros(Ntar*ncol,1)];
  f2 = [f2(:);zeros(Ntar*ncol,1)];
  % pad density function with zeros so that Xtar doesn't
  % affect the double-layer potential
  nx = [nx(:);zeros(Ntar*ncol,1)];
  ny = [ny(:);zeros(Ntar*ncol,1)];
  % pad the normal vctor with zeros so that Xtar doesn't
  % affect the double-layer potential

  [u,v] = stokesDLPfmm(f1,f2,x,y,nx,ny);
  stokesDLPtar = zeros(2*Ntar,ncol); % initialize
  for k = 1:ncol
    is = N*numel(K) + (k-1)*Ntar+1;
    ie = is + Ntar - 1;
    stokesDLPtar(1:Ntar,k) = u(is:ie);
    stokesDLPtar(Ntar+1:2*Ntar,k) = v(is:ie);
  end
  % Wrap the output of the FMM in the usual format
  % for the target points
end


end % exactStokesDLfmm




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [laplaceDLP,laplaceDLPtar] = ...
    exactLaplaceDLfmm(o,vesicle,f,Xtar,K)
% [laplaceDLP,laplaceDLPtar] = exactLaplaceDLfmm(vesicle,f) uses 
% the FMM to compute the double-layer laplace potential due to all 
% vesicles except itself vesicle is a class of object capsules and 
% f is the density function NOT scaled by arclength term.  Xtar is 
% a set of points where the single-layer potential due to all 
% vesicles in index set K

N = vesicle.N; % number of points per vesicle
nv = vesicle.nv; % number of vesicles
X = vesicle.X; % vesicle positions
normal = vesicle.normal; % vesicle's outward normal
oc = curve;
[x,y] = oc.getXY(X); % seperate x and y coordinates
[nx,ny] = oc.getXY(normal); % seperate x and y coordinates of normal

[f,~] = oc.getXY(f);
den = f.*vesicle.sa*2*pi/N;

if (nargin == 5)
  laplaceDLP = [];
else
  [~,rfield,~] = fmm_laplace(nx(:).*den(:),x(:),y(:),1);
  [~,~,cfield] = fmm_laplace(ny(:).*den(:),x(:),y(:),1);
  % double-layer potential due to all vesicles
  potTot = rfield + cfield;

  laplaceDLP = zeros(N,nv); % initialize
  for k = 1:nv
    is = (k-1)*N+1;
    ie = k*N;
    laplaceDLP(1:N,k) = potTot(is:ie);
  end
  % Wrap the output of the FMM into the usual 
  % [[x1;y1] [x2;y2] ...] format

  for k = 1:nv
    [~,rfield,~] = fmm_laplace(nx(:,k).*den(:,k),x(:,k),y(:,k),1);
    [~,~,cfield] = fmm_laplace(ny(:,k).*den(:,k),x(:,k),y(:,k),1);
    laplaceDLP(:,k) = laplaceDLP(:,k) - rfield - cfield;
  end
  % Subtract potential due to each vesicle on its own.  Nothing
  % to change here for potential at Xtar
end
laplaceDLP = -[laplaceDLP;zeros(vesicle.N,vesicle.nv)]/(2*pi);
% pad with zeros so that it is compatible with Stokes near-singular
% integration

if nargin == 3
  laplaceDLPtar = [];
else
  x = x(:,K);
  y = y(:,K);
  nx = nx(:,K);
  ny = ny(:,K);
  den = den(:,K);
  % only care about points on vesicles indexed by K
  [Ntar,ncol] = size(Xtar);
  Ntar = Ntar/2;
  x2 = Xtar(1:Ntar,:);
  x = [x(:);x2(:)];
  y2 = Xtar(Ntar+1:2*Ntar,:);
  y = [y(:);y2(:)];
  % Stack the x and y coordinates of the target points
  den = [den(:);zeros(Ntar*ncol,1)];
  nx = [nx(:);zeros(Ntar*ncol,1)];
  ny = [ny(:);zeros(Ntar*ncol,1)];
  % pad density function and normal with zeros so that Xtar doesn't
  % affect the single-layer potential
  [~,rfield,~] = fmm_laplace(nx.*den,x(:),y(:),1);
  [~,~,cfield] = fmm_laplace(ny.*den,x(:),y(:),1);
  % evalaute single-layer potential using fmm

  laplaceDLPtar = zeros(2*Ntar,ncol); % initialize
  for k = 1:ncol
    is = N*numel(K) + (k-1)*Ntar+1;
    ie = is + Ntar - 1;
    laplaceDLPtar(1:Ntar,k) = rfield(is:ie)+cfield(is:ie);
  end
  % Wrap the output of the FMM in the usual format
  % for the target points
  laplaceDLPtar = -laplaceDLPtar/(2*pi);
end


end % exactLaplaceDLfmm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pot = fmmStokes(o,x,y,f1,f2)
% pot = fmm_stokes(x,y,f1,f2) computes the single-layer potential
% at points (x,y) due to the density function (f1,f2).  Only diagonal
% term is ommited.

[pot1,rfield,cfield] = fmm_laplace(f1,x,y,2);
pot = [-pot1 + x.*rfield;x.*cfield];
% Logarithm term and part of the r \otimes r term

[pot1,rfield,cfield] = fmm_laplace(f2,x,y,2);
pot = pot + [y.*rfield;-pot1 + y.*cfield];
% Logarithm term and part of the r \otimes r term

[~,rfield,cfield] = fmm_laplace(x.*f1 + y.*f2,x,y,1);
pot = pot - [rfield;cfield];
% Other part of r \otimes r term

pot = 0.25*pot/pi;
% single-layer potential has 1/4/pi term in front

end % fmmStokes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [laplaceSLP,laplaceSLPtar] = ...
    exactLaplaceSL(o,vesicle,f,Xtar,K1)
% pot = exactLaplaceSL(vesicle,f,Xtar,K1) 
% computes the single-layer laplace potential due to f around all
% vesicles except itself.  Also can pass a set of target points Xtar
% and a collection of vesicles K1 and the single-layer potential due 
% to vesicles in K1 will be evaluated at Xtar.
% Everything but Xtar is in the 2*N x nv format
% Xtar is in the 2*Ntar x ncol format

oc = curve;
[x,y] = oc.getXY(vesicle.X);
% vesicle position

if nargin == 5
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  laplaceSLPtar = zeros(2*Ntar,ncol);
else
  K1 = [];
  laplaceSLPtar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 5, user does not need the layer potential at arbitrary
  % points
end

den = f.*[vesicle.sa;vesicle.sa]*2*pi/vesicle.N;
% multiply by arclength term

for k2 = 1:ncol % loop over columns of target points
  for j = 1:Ntar % loop over rows of target points
    dis2 = (x(:,K1) - Xtar(j,k2)).^2 + (y(:,K1) - Xtar(j+Ntar,k2)).^2;
    % distance squared

    coeff = 0.5*log(dis2);
    % this is the kernel of the single-layer potential for
    % laplace's equation

    val = coeff.*den(1:vesicle.N,K1);
    laplaceSLPtar(j,k2) = sum(val(:));
    val = coeff.*den(vesicle.N+1:2*vesicle.N,K1);
    laplaceSLPtar(j+Ntar,k2) = sum(val(:));
  end % j

end % k2
% Evaluate single-layer potential at arbitrary target points
laplaceSLPtar = -1/(4*pi)*laplaceSLPtar;
% -1/4/pi is the coefficient in front of the single-layer potential as
% defined in poten.laplaceSLmatrix

laplaceSLP = zeros(vesicle.N,vesicle.nv); % Initialize to zero
% if we only have one vesicle, vesicles of course can not collide
% Don't need to run this loop in this case
if (nargin == 4 && vesicle.nv > 1)
  for k1 = 1:vesicle.nv % vesicle of targets
    K = [(1:k1-1) (k1+1:vesicle.nv)];
    % Loop over all vesicles except k1

    for j=1:vesicle.N
      dis2 = (x(:,K) - x(j,k1)).^2 + (y(:,K) - y(j,k1)).^2;
      % distance squared

      coeff = 0.5*log(dis2);
      % this is the kernel of the single-layer potential for
      % laplace's equation
      val = coeff.*den(1:vesicle.N,K);
      laplaceSLP(j,k1) = sum(val(:));
    end % j
  end % k1
  % Evaluate single-layer potential at vesicles but oneself
  laplaceSLP = -1/(4*pi)*laplaceSLP;
  % -1/4/pi is the coefficient in front of the single-layer potential as
  % defined in poten.laplaceSLmatrix
end % nargin == 4

laplaceSLP = [laplaceSLP;laplaceSLP];

end % exactLaplaceSL


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [laplaceDLP,laplaceDLPtar] = ...
    exactLaplaceDL(o,vesicle,f,Xtar,K1)
% pot = exactLaplaceDL(vesicle,f,Xtar,K1) 
% computes the double-layer laplace potential due to f around all
% vesicles except itself.  Also can pass a set of target points Xtar
% and a collection of vesicles K1 and the double-layer potential due 
% to vesicles in K1 will be evaluated at Xtar.
% Everything but Xtar is in the 2*N x nv format
% Xtar is in the 2*Ntar x ncol format

nv = vesicle.nv; % number of vesicles
X = vesicle.X; % Vesicle positions
normal = vesicle.normal; % Outward normal
sa = vesicle.sa; % Jacobian

if nargin == 5
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  laplaceDLPtar = zeros(2*Ntar,ncol);
else
  K1 = [];
  laplaceDLPtar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 5, user does not need the layer potential at arbitrary
  % points
end

den = f.*[sa;sa]*2*pi/vesicle.N;
% multiply by arclength term

for k2 = 1:ncol % loop over columns of target points
  for j = 1:Ntar % loop over rows of target points
    dis2 = (X(1:vesicle.N,K1) - Xtar(j,k2)).^2 + ... 
        (X(vesicle.N+1:2*vesicle.N,K1) - Xtar(j+Ntar,k2)).^2;
    diffxy = [X(1:vesicle.N,K1) - Xtar(j,k2) ; ...
        X(vesicle.N+1:2*vesicle.N,K1) - Xtar(j+Ntar,k2)];
    % distance squared and difference of source and target location

    coeff = (diffxy(1:vesicle.N,:).*normal(1:vesicle.N,K1) + ...
      diffxy(vesicle.N+1:2*vesicle.N,:).*...
      normal(vesicle.N+1:2*vesicle.N,K1))./dis2;
    % this is the kernel of the double-layer potential for
    % laplace's equation
    val = coeff.*den(1:vesicle.N,K1);
    laplaceDLPtar(j,k2) = sum(val(:));
    val = coeff.*den(vesicle.N+1:2*vesicle.N,K1);
    laplaceDLPtar(j+Ntar,k2) = sum(val(:));
  end % j

end % k2
% Evaluate double-layer potential at arbitrary target points
laplaceDLPtar = 1/(2*pi)*laplaceDLPtar;
% 1/2/pi is the coefficient in front of the double-layer potential

laplaceDLP = zeros(vesicle.N,vesicle.nv); % Initialize to zero
% if we only have one vesicle, vesicles of course can not collide
% Don't need to run this loop in this case
if (nargin == 3 && vesicle.nv > 1)
  for k1 = 1:vesicle.nv % vesicle of targets
    K = [(1:k1-1) (k1+1:vesicle.nv)];
    % Loop over all vesicles except k1

    for j=1:vesicle.N
      dis2 = (X(1:vesicle.N,K) - X(j,k1)).^2 + ...
          (X(vesicle.N+1:2*vesicle.N,K) - X(j+vesicle.N,k1)).^2;
      diffxy = [X(1:vesicle.N,K) - X(j,k1); ...
          X(vesicle.N+1:2*vesicle.N,K) - X(j+vesicle.N,k1)];
      % distance squared and difference of source and target location

      coeff = (diffxy(1:vesicle.N,:).*normal(1:vesicle.N,K) + ...
        diffxy(vesicle.N+1:2*vesicle.N,:).*...
        normal(vesicle.N+1:2*vesicle.N,K))./dis2;
      % this is the kernel of the double-layer potential for
      % laplace's equation
      val = coeff.*den(1:vesicle.N,K);
      laplaceDLP(j,k1) = sum(val(:));
    end % j
  end % k1
  % Evaluate double-layer potential at vesicles but oneself
  laplaceDLP = 1/(2*pi)*laplaceDLP;
  % 1/2/pi is the coefficient in front of the double-layer potential
end % nargin == 3

laplaceDLP = [laplaceDLP;laplaceDLP];

end % exactLaplaceDL



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pressSLP,pressSLPtar] = exactPressSL(o,vesicle,f,Xtar,K1)
% [pressSLP,pressSLPtar] = exactPressSL(vesicle,f,Xtar,K1)
% computes the pressure due to all vesicles contained in vesicle
% and indexed over K1.  Evaluates it at Xtar
% Everything but Xtar is in the 2*N x nv format
% Xtar is in the 2*Ntar x ncol format

N = vesicle.N; % number of points per vesicle
nv = vesicle.nv; % number of vesicles
X = vesicle.X; % Vesicle positions
sa = vesicle.sa; % Jacobian

if nargin == 5
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  pressSLPtar = zeros(Ntar,ncol);
else
  K1 = [];
  pressSLPtar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 5, user does not need the layer potential at arbitrary
  % points
end

den = f.*[sa;sa]*2*pi/vesicle.N;

for k2 = 1:ncol % loop over columns of target points
  for j = 1:Ntar % loop over rows of target points
    dis2 = (Xtar(j,k2) - X(1:vesicle.N,K1)).^2 + ... 
        (Xtar(j+Ntar,k2) - X(vesicle.N+1:2*vesicle.N,K1)).^2;
    diffxy = [Xtar(j,k2) - X(1:vesicle.N,K1) ; ...
        Xtar(j+Ntar,k2) - X(vesicle.N+1:2*vesicle.N,K1)];
    % distance squared and difference of source and target location

    val = (diffxy(1:vesicle.N,:).*den(1:vesicle.N,K1) + ...
        diffxy(vesicle.N+1:2*vesicle.N,:).* ...
        den(vesicle.N+1:2*vesicle.N,K1))./dis2;
    % \frac{(r \dot f){\rho^{2}} term

    pressSLPtar(j,k2) = sum(val(:)); 
  end % j
end % k2
% pressure coming from the single-layer potential for Stokes flow

pressSLP = zeros(N,nv);
% TODO: NOT SURE IF WE WILL EVER NEED THIS BUT SHOULD PUT IT
% IN NONETHELESS

pressSLPtar = [pressSLPtar;pressSLPtar]*1/2/pi;
% near-singular integration needs vector-valued functions
% also need to multiply by 1/(2*pi) as per the pressure of the
% single-layer potential

end % exactPressSL





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pressSLP,pressSLPtar] = exactPressSLfmm(o,vesicle,f,Xtar,K1)
% [pressSLP,pressSLPtar] = exactPressSLfmm(vesicle,f,Xtar,K1)
% computes the pressure due to all vesicles contained in vesicle
% and indexed over K1 using the FMM.  Evaluates it at Xtar
% Everything but Xtar is in the 2*N x nv format
% Xtar is in the 2*Ntar x ncol format

N = vesicle.N; % number of points per vesicle
nv = vesicle.nv; % number of vesicles
X = vesicle.X; % Vesicle positions
sa = vesicle.sa; % Jacobian

if nargin == 5
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  pressSLPtar = zeros(Ntar,ncol);
else
  K1 = [];
  pressSLPtar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 5, user does not need the layer potential at arbitrary
  % points
end

den = f.*[sa;sa]*2*pi/vesicle.N;

x1 = X(1:vesicle.N,K1);
y1 = X(vesicle.N+1:2*vesicle.N,K1);
% source points
x2 = Xtar(1:Ntar,:);
y2 = Xtar(Ntar+1:2*Ntar,:);
% target points
den1 = den(1:vesicle.N,K1);
den2 = den(vesicle.N+1:2*vesicle.N,K1);
% source charges

x = [x1(:); x2(:)];
y = [y1(:); y2(:)];
% stack the sources and targets
den1 = [den1(:);zeros(Ntar*ncol,1)];
den2 = [den2(:);zeros(Ntar*ncol,1)];
% charge of zero at the target locations
[~,rfield,~] = fmm_laplace(den1,x,y,1);
[~,~,cfield] = fmm_laplace(den2,x,y,1);

pressSLPtar = zeros(Ntar,ncol);
for k = 1:ncol
  is = N*numel(K1) + (k-1)*Ntar + 1;
  ie = is + Ntar - 1;
  pressSLPtar(1:Ntar,k) = rfield(is:ie) + cfield(is:ie);
end

pressSLP = zeros(N,nv);
% TODO: NOT SURE IF WE WILL EVER NEED THIS BUT SHOULD PUT IT
% IN NONETHELESS

pressSLPtar = [pressSLPtar;pressSLPtar]*1/2/pi;
% near-singular integration needs vector-valued functions
% also need to multiply by 1/(2*pi) as per the pressure of the
% single-layer potential

end % exactPressSLfmm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pressDLP,pressDLPtar] = exactPressDL(o,vesicle,f,Xtar,K1)
% [pressDLP,pressDLPtar] = exactPressDL(vesicle,f,Xtar,K1)
% computes the pressure due to all vesicles contained in vesicle
% and indexed over K1.  Evaluates it at Xtar
% Everything but Xtar is in the 2*N x nv format
% Xtar is in the 2*Ntar x ncol format

N = vesicle.N; % number of points per vesicle
nv = vesicle.nv; % number of vesicles
X = vesicle.X; % Vesicle positions
normal = vesicle.normal; % Normal vector
sa = vesicle.sa; % Jacobian

if nargin == 5
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  pressDLPtar = zeros(Ntar,ncol);
else
  K1 = [];
  pressDLPtar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 5, user does not need the layer potential at arbitrary
  % points
end

den = f.*[sa;sa]*2*pi/vesicle.N;

for k2 = 1:ncol % loop over columns of target points
  for j = 1:Ntar % loop over rows of target points
    dis2 = (Xtar(j,k2) - X(1:N,K1)).^2 + ... 
        (Xtar(j+Ntar,k2) - X(N+1:2*N,K1)).^2;
    diffxy = [Xtar(j,k2) - X(1:N,K1); ...
        Xtar(j+Ntar,k2) - X(N+1:2*N,K1)];
    % distance squared and difference of source and target location
    rdotn = diffxy(1:N,:).*normal(1:N,K1) + ...
        diffxy(N+1:2*N,:).*normal(N+1:2*N,K1); 
  
    val = (normal(1:N,K1) - 2*rdotn./dis2.*...
        diffxy(1:N,:))./dis2 .* den(1:N,K1);
    val = val + (normal(N+1:2*N,K1) - 2*rdotn./dis2.*...
        diffxy(N+1:2*N,:))./dis2 .* den(N+1:2*N,K1);

    pressDLPtar(j,k2) = sum(val(:)); 
  end % j
end % k2
%% pressure coming from the double-layer potential for Stokes flow

pressDLP = zeros(N,nv);
% TODO: NOT SURE IF WE WILL EVER NEED THIS BUT SHOULD PUT IT
% IN NONETHELESS

pressDLPtar = -[pressDLPtar;pressDLPtar]*1/pi;
% near-singular integration needs vector-valued functions
% also need to multiply by 1/(2*pi) as per the pressure of the
% single-layer potential

end % exactPressDL


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stressSLP,stressSLPtar] = exactStressSL1(o,vesicle,f,Xtar,K1)
% [stressSLP,stressSLPtar] = exactStressSL1(vesicle,f,Xtar,K1)
% computes the stress due to all vesicles contained in vesicle
% and indexed over K1.  Only computes the stress applied to the
% direction [1;0].  Evaluates it at Xtar
% Everything but Xtar is in the 2*N x nv format
% Xtar is in the 2*Ntar x ncol format

N = vesicle.N; % number of points per vesicle
nv = vesicle.nv; % number of vesicles
X = vesicle.X; % Vesicle positions
sa = vesicle.sa; % Jacobian

if nargin == 5
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  stressSLPtar = zeros(2*Ntar,ncol);
else
  K1 = [];
  stressSLPtar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 5, user does not need the layer potential at arbitrary
  % points
end

den = f.*[sa;sa]*2*pi/vesicle.N;

for k2 = 1:ncol % loop over columns of target points
  for j = 1:Ntar % loop over rows of target points
    dis4 = ((Xtar(j,k2) - X(1:vesicle.N,K1)).^2 + ... 
        (Xtar(j+Ntar,k2) - X(vesicle.N+1:2*vesicle.N,K1)).^2).^2;
    diffxy = [Xtar(j,k2) - X(1:vesicle.N,K1) ; ...
        Xtar(j+Ntar,k2) - X(vesicle.N+1:2*vesicle.N,K1)];
    % distance squared and difference of source and target location
    rdotf = diffxy(1:vesicle.N,:).*den(1:vesicle.N,K1) + ...
        diffxy(vesicle.N+1:2*vesicle.N,:).*...
        den(vesicle.N+1:2*vesicle.N,K1);

    val = rdotf./dis4.*diffxy(1:vesicle.N,:).*...
        diffxy(1:vesicle.N,:);
    % \frac{r \dot f}{\rho^{4}}*r1*r1 term
    stressSLPtar(j,k2) = sum(val(:)); 

    val = rdotf./dis4.*diffxy(1:vesicle.N,:).*...
        diffxy(vesicle.N+1:2*vesicle.N,:);
    % \frac{r \dot f}{\rho^{4}}*r1*r2 term
    stressSLPtar(j+Ntar,k2) = sum(val(:)); 
  end % j
end % k2
% stress coming from the single-layer potential for Stokes flow

stressSLP = zeros(2*N,nv);
% TODO: NOT SURE IF WE WILL EVER NEED THIS BUT SHOULD PUT IT
% IN NONETHELESS

stressSLPtar = stressSLPtar/pi;
% 1/pi is the constant in front of the stress of the single-layer 
% potential

end % exactStressSL1



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stressSLP,stressSLPtar] = exactStressSL2(o,vesicle,f,Xtar,K1)
% [stressSLP,stressSLPtar] = exactStressSL2(vesicle,f,Xtar,K1)
% computes the stress due to all vesicles contained in vesicle
% and indexed over K1.  Only computes the stress applied to the
% direction [0;1].  Evaluates it at Xtar
% Everything but Xtar is in the 2*N x nv format
% Xtar is in the 2*Ntar x ncol format

N = vesicle.N; % number of points per vesicle
nv = vesicle.nv; % number of vesicles
X = vesicle.X; % Vesicle positions
sa = vesicle.sa; % Jacobian

if nargin == 5
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  stressSLPtar = zeros(2*Ntar,ncol);
else
  K1 = [];
  stressSLPtar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 5, user does not need the layer potential at arbitrary
  % points
end

den = f.*[sa;sa]*2*pi/vesicle.N;

for k2 = 1:ncol % loop over columns of target points
  for j = 1:Ntar % loop over rows of target points
    dis4 = ((Xtar(j,k2) - X(1:vesicle.N,K1)).^2 + ... 
        (Xtar(j+Ntar,k2) - X(vesicle.N+1:2*vesicle.N,K1)).^2).^2;
    diffxy = [Xtar(j,k2) - X(1:vesicle.N,K1) ; ...
        Xtar(j+Ntar,k2) - X(vesicle.N+1:2*vesicle.N,K1)];
    % distance squared and difference of source and target location
    rdotf = diffxy(1:vesicle.N,:).*den(1:vesicle.N,K1) + ...
        diffxy(vesicle.N+1:2*vesicle.N,:).*...
        den(vesicle.N+1:2*vesicle.N,K1);

    val = rdotf./dis4.*diffxy(1:vesicle.N,:).*...
        diffxy(vesicle.N+1:2*vesicle.N,:);
    % \frac{r \dot f}{\rho^{4}}*r1*r1 term
    stressSLPtar(j,k2) = sum(val(:)); 

    val = rdotf./dis4.*diffxy(vesicle.N+1:2*vesicle.N,:).*...
        diffxy(vesicle.N+1:2*vesicle.N,:);
    % \frac{r \dot f}{\rho^{4}}*r1*r2 term
    stressSLPtar(j+Ntar,k2) = sum(val(:)); 
  end % j
end % k2
% stress coming from the single-layer potential for Stokes flow

stressSLP = zeros(2*N,nv);
% TODO: NOT SURE IF WE WILL EVER NEED THIS BUT SHOULD PUT IT
% IN NONETHELESS

stressSLPtar = stressSLPtar/pi;
% 1/pi is the constant in front of the stress of the single-layer 
% potential

end % exactStressSL2



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stressDLP,stressDLPtar] = exactStressDL1(o,vesicle,f,Xtar,K1)
% [stressDLP,stressDLPtar] = exactStressDL1(vesicle,f,Xtar,K1)
% computes the stress due to the double-layer potential of all vesicles 
% contained in vesicle and indexed over K1.  Only computes the stress 
% applied to the direction [1;0].  Evaluates it at Xtar
% Everything but Xtar is in the 2*N x nv format
% Xtar is in the 2*Ntar x ncol format

oc = curve;
N = vesicle.N; % number of points per vesicle
nv = vesicle.nv; % number of vesicles
X = vesicle.X; % Vesicle positions
sa = vesicle.sa; % Jacobian
normal = vesicle.normal; % normal vector
[nx,ny] = oc.getXY(normal);

if nargin == 5
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  stressDLPtar = zeros(2*Ntar,ncol);
else
  K1 = [];
  stressDLPtar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 5, user does not need the layer potential at arbitrary
  % points
end

[fx,fy] = oc.getXY(f);
% first and second components of the density function

for k2 = 1:ncol % loop over columns of target points
  for j = 1:Ntar % loop over rows of target points
    dis2 = (Xtar(j,k2) - X(1:vesicle.N,K1)).^2 + ... 
        (Xtar(j+Ntar,k2) - X(vesicle.N+1:2*vesicle.N,K1)).^2;
    diffxy = [Xtar(j,k2) - X(1:vesicle.N,K1) ; ...
        Xtar(j+Ntar,k2) - X(vesicle.N+1:2*vesicle.N,K1)];
    % distance squared and difference of source and target location
    [rx,ry] = oc.getXY(diffxy);

    rdotf = rx.*fx(:,K1) + ry.*fy(:,K1);
    % dot product of r and f
    fdotn = fx(:,K1).*nx(:,K1) + fy(:,K1).*ny(:,K1);
    % dot product of f and n
    rdotn = rx.*nx(:,K1) + ry.*ny(:,K1);
    % dot product of r and n

    val = (fdotn./dis2 - ...
        8./dis2.^3.*rdotn.*rdotf.*rx.*rx + ...
        rdotn./dis2.^2.*(2*rx.*fx(:,K1)) + ...
        rdotf./dis2.^2.*(2*rx.*nx(:,K1))).*sa(:,K1);
    % first component of the stress of the double-layer potential
    % applied to [1;0]

    stressDLPtar(j,k2) = sum(val(:)); 
    % scale by arclength
    stressDLPtar(j,k2) = stressDLPtar(j,k2)*2*pi/N;
    % d\theta term

    val = (-8./dis2.^3.*rdotn.*rdotf.*rx.*ry + ...
        rdotn./dis2.^2.*(rx.*fy(:,K1) + ry.*fx(:,K1)) + ...
        rdotf./dis2.^2.*(rx.*ny(:,K1) + ry.*nx(:,K1))).*sa(:,K1);
    % second component of the stress of the double-layer potential
    % applied to [1;0]

    stressDLPtar(j+Ntar,k2) = sum(val(:)); 
    % scale by arclength
    stressDLPtar(j+Ntar,k2) = stressDLPtar(j+Ntar,k2)*2*pi/N;
    % d\theta term

  end % j
end % k2
% stress coming from the single-layer potential for Stokes flow

stressDLP = zeros(2*N,nv);
% TODO: NOT SURE IF WE WILL EVER NEED THIS BUT SHOULD PUT IT
% IN NONETHELESS

stressDLPtar = stressDLPtar/pi;
% 1/pi is the constant in front of the stress of the double-layer 
% potential

end % exactStressDL1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stressDLP,stressDLPtar] = exactStressDL2(o,vesicle,f,Xtar,K1)
% [stressDLP,stressDLPtar] = exactStressDL2(vesicle,f,Xtar,K1)
% computes the stress due to the double-layer potential of all vesicles 
% contained in vesicle and indexed over K1.  Only computes the stress 
% applied to the direction [0;1].  Evaluates it at Xtar
% Everything but Xtar is in the 2*N x nv format
% Xtar is in the 2*Ntar x ncol format

oc = curve;
N = vesicle.N; % number of points per vesicle
nv = vesicle.nv; % number of vesicles
X = vesicle.X; % Vesicle positions
sa = vesicle.sa; % Jacobian
normal = vesicle.normal; % normal vector
[nx,ny] = oc.getXY(normal);

if nargin == 5
  Ntar = size(Xtar,1)/2;
  ncol = size(Xtar,2);
  stressDLPtar = zeros(2*Ntar,ncol);
else
  K1 = [];
  stressDLPtar = [];
  ncol = 0;
  Ntar = 0;
  % if nargin ~= 5, user does not need the layer potential at arbitrary
  % points
end

[fx,fy] = oc.getXY(f);
% first and second components of the density function

for k2 = 1:ncol % loop over columns of target points
  for j = 1:Ntar % loop over rows of target points
    dis2 = (Xtar(j,k2) - X(1:vesicle.N,K1)).^2 + ... 
        (Xtar(j+Ntar,k2) - X(vesicle.N+1:2*vesicle.N,K1)).^2;
    diffxy = [Xtar(j,k2) - X(1:vesicle.N,K1) ; ...
        Xtar(j+Ntar,k2) - X(vesicle.N+1:2*vesicle.N,K1)];
    % distance squared and difference of source and target location
    [rx,ry] = oc.getXY(diffxy);

    rdotf = rx.*fx(:,K1) + ry.*fy(:,K1);
    % dot product of r and f
    fdotn = fx(:,K1).*nx(:,K1) + fy(:,K1).*ny(:,K1);
    % dot product of f and n
    rdotn = rx.*nx(:,K1) + ry.*ny(:,K1);
    % dot product of r and n

    val = (-8./dis2.^3.*rdotn.*rdotf.*ry.*rx + ...
        rdotn./dis2.^2.*(rx.*fy(:,K1) + ry.*fx(:,K1)) + ...
        rdotf./dis2.^2.*(rx.*ny(:,K1) + ry.*nx(:,K1))).*sa(:,K1);
    % second component of the stress of the double-layer potential
    % applied to [0;1]

    stressDLPtar(j,k2) = sum(val(:)); 
    stressDLPtar(j,k2) = stressDLPtar(j,k2)*2*pi/N;
    % d\theta term

    val = (fdotn./dis2 - ...
        8./dis2.^3.*rdotn.*rdotf.*ry.*ry + ...
        rdotn./dis2.^2.*(2*ry.*fy(:,K1)) + ...
        rdotf./dis2.^2.*(2*ry.*ny(:,K1))).*sa(:,K1);
    % first component of the stress of the double-layer potential
    % applied to [0;1]

    stressDLPtar(j+Ntar,k2) = sum(val(:)); 
    stressDLPtar(j+Ntar,k2) = stressDLPtar(j+Ntar,k2)*2*pi/N;
    % d\theta term
  end % j
end % k2
% stress coming from the single-layer potential for Stokes flow

stressDLP = zeros(2*N,nv);
% TODO: NOT SURE IF WE WILL EVER NEED THIS BUT SHOULD PUT IT
% IN NONETHELESS

stressDLPtar = stressDLPtar/pi;
% 1/pi is the constant in front of the stress of the double-layer 
% potential


end % exactStressDL2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [press,stress1,stress2] = pressAndStress(o, ...
    X,sigma,u,kappa,viscCont,...
    walls,pressTar,eta,RS,confined,...
    near,fmm,om)
% [press,stress1,stress2] = op.PressAndStress(...
%    X,sigma,u,kappa,viscCont,...
%    walls,pressTar,eta,RS,confined,...
%    near,fmm,om)
% computes the total pressure and stress due vesicles stored
% in X and the solid walls stored in walls.X.  X is the position
% of the vesicles, sigma is the tension, u is the velocity, kappa
% is the bending coefficient, viscCont is the viscosity contrast, 
% walls is a object for solid walls, pressTar is the target locations
% where the pressure and stress are to be evaluated, eta and RS are
% the density function, rotlets, and stokeslets, confined and near and
% fmm are flags for different algorithms, and om is a monitor object

vesicle = capsules(X,sigma,u,kappa,viscCont);
f = vesicle.tracJump(vesicle.X,vesicle.sig);
% compute traction jump of vesicle

press = vesicle.pressure(f,[],pressTar,near,...
    fmm,'SLP');
% compute pressure due to vesicles

[stress1,stress2] = vesicle.stressTensor(f,...
    [],pressTar,near,fmm,'SLP');
% compute stress due to vesicles at same location as pressure

if confined
  press = press + walls.pressure(...
      eta,RS,pressTar,near,fmm,'DLP');
  % compute pressure due to solid walls and the Stokeslet
  % The Rotlet has a vanishing pressure
  [stress1T,stress2T] = walls.stressTensor(eta,...
      RS,pressTar,near,fmm,'DLP');
  % compute stress due to solid walls.  This includes the
  % stress due to the rotlet and stokeslet terms
  stress1 = stress1 + stress1T;
  stress2 = stress2 + stress2T;
  % update stress
end
% compute pressure and stress due to solid walls if they exist

if om.saveData
  om.writePressure(press);
  om.writeStress(stress1(1:end/2),'11');
  om.writeStress(stress1(1+end/2:end),'12');
  om.writeStress(stress2(1:end/2),'21');
  om.writeStress(stress2(1+end/2:end),'22');
  % write pressure and different components of stress
end

end % pressAndStress

end % methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods (Static)

function LP = lagrangeInterp
% interpMap = lagrangeInterp builds the Lagrange interpolation
% matrix that takes seven function values equally distributed
% in [0,1] and returns the seven polynomial coefficients

interpMat = zeros(7);
LP(1,1) = 6.48e1;
LP(1,2) = -3.888e2;
LP(1,3) = 9.72e2;
LP(1,4) = -1.296e3;
LP(1,5) = 9.72e2;
LP(1,6) = -3.888e2;
LP(1,7) = 6.48e1;

LP(2,1) = -2.268e2;
LP(2,2) = 1.296e3;
LP(2,3) = -3.078e3;
LP(2,4) = 3.888e3;
LP(2,5) = -2.754e3;
LP(2,6) = 1.0368e3;
LP(2,7) = -1.62e2;

LP(3,1) = 3.15e2;
LP(3,2) = -1.674e3;
LP(3,3) = 3.699e3;
LP(3,4) = -4.356e3;
LP(3,5) = 2.889e3;
LP(3,6) = -1.026e3;
LP(3,7) = 1.53e2;

LP(4,1) = -2.205e2;
LP(4,2) = 1.044e3;
LP(4,3) = -2.0745e3;
LP(4,4) = 2.232e3;
LP(4,5) = -1.3815e3;
LP(4,6) = 4.68e2;
LP(4,7) = -6.75e1;

LP(5,1) = 8.12e1;
LP(5,2) = -3.132e2;
LP(5,3) = 5.265e2;
LP(5,4) = -5.08e2;
LP(5,5) = 2.97e2;
LP(5,6) = -9.72e1;
LP(5,7) = 1.37e1;

LP(6,1) = -1.47e1;
LP(6,2) = 3.6e1;
LP(6,3) = -4.5e1;
LP(6,4) = 4.0e1;
LP(6,5) = -2.25e1;
LP(6,6) = 7.2e0;
LP(6,7) = -1e0;

LP(7,1) = 1e0;
% rest of the coefficients are zero


end % lagrangeInterp


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xp = nodesLog

xp = zeros(25,2);
xp(1,1) = 2.379647284118974e-2;
xp(1,2) = 8.795942675593887e-2;
xp(2,1) = 2.935370741501914e-1;
xp(2,2) = 4.989017152913699e-1;
xp(3,1) = 1.023715124251890e+0;
xp(3,2) = 9.131388579526912e-1;
xp(4,1) = 6.531815708567918e-3;
xp(4,2) = 2.462194198995203e-2;
xp(5,1) = 9.086744584657729e-2;
xp(5,2) = 1.701315866854178e-1;
xp(6,1) = 3.967966533375878e-1;
xp(6,2) = 4.609256358650077e-1;
xp(7,1) = 1.027856640525646e+0;
xp(7,2) = 7.947291148621895e-1;
xp(8,1) = 1.945288592909266e+0;
xp(8,2) = 1.008710414337933e+0;
xp(9,1) = 2.980147933889640e+0;
xp(9,2) = 1.036093649726216e+0;
xp(10,1) = 3.998861349951123e+0;
xp(10,2) = 1.004787656533285e+0;
xp(11,1) = 8.371529832014113e-4;
xp(11,2) = 3.190919086626234e-3;
xp(12,1) = 1.239382725542637e-2;
xp(12,2) = 2.423621380426338e-2;
xp(13,1) = 6.009290785739468e-2;
xp(13,2) = 7.740135521653088e-2;
xp(14,1) = 1.805991249601928e-1;
xp(14,2) = 1.704889420286369e-1;
xp(15,1) = 4.142832599028031e-1;
xp(15,2) = 3.029123478511309e-1;
xp(16,1) = 7.964747731112430e-1;
xp(16,2) = 4.652220834914617e-1;
xp(17,1) = 1.348993882467059e+0;
xp(17,2) = 6.401489637096768e-1;
xp(18,1) = 2.073471660264395e+0;
xp(18,2) = 8.051212946181061e-1;
xp(19,1) = 2.947904939031494e+0;
xp(19,2) = 9.362411945698647e-1;
xp(20,1) = 3.928129252248612e+0;
xp(20,2) = 1.014359775369075e+0;
xp(21,1) = 4.957203086563112e+0;
xp(21,2) = 1.035167721053657e+0;
xp(22,1) = 5.986360113977494e+0;
xp(22,2) = 1.020308624984610e+0;
xp(23,1) = 6.997957704791519e+0;
xp(23,2) = 1.004798397441514e+0;
xp(24,1) = 7.999888757524622e+0;
xp(24,2) = 1.000395017352309e+0;
xp(25,1) = 8.999998754306120e+0;
xp(25,2) = 1.000007149422537e+0;

end % nodesLog


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xp = nodesRegular

xp = zeros(35,2);
xp(1,1) = 2.000000000000000e-1;
xp(1,2) = 5.208333333333333e-1;
xp(2,1) = 1.000000000000000e+0;
xp(2,2) = 9.791666666666667e-1;
xp(3,1) = 2.250991042610971e-1;
xp(3,2) = 5.549724327164181e-1;
xp(4,1) = 1.014269060987992e+0;
xp(4,2) = 9.451317411845474e-1;
xp(5,1) = 2.000000000000000e+0;
xp(5,2) = 9.998958260990347e-1;
xp(6,1) = 2.087647422032129e-1;
xp(6,2) = 5.207988277246498e-1;
xp(7,1) = 9.786087373714483e-1;
xp(7,2) = 9.535038018555888e-1;
xp(8,1) = 1.989541386579751e+0;
xp(8,2) = 1.024871626402471e+0;
xp(9,1) = 3.000000000000000e+0;
xp(9,2) = 1.000825744017291e+0;
xp(10,1) = 7.023955461621939e-2;
xp(10,2) = 1.922315977843698e-1;
xp(11,1) = 4.312297857227970e-1;
xp(11,2) = 5.348399530514687e-1;
xp(12,1) = 1.117752734518115e+0;
xp(12,2) = 8.170209442488760e-1;
xp(13,1) = 2.017343724572518e+0;
xp(13,2) = 9.592111521445966e-1;
xp(14,1) = 3.000837842847590e+0;
xp(14,2) = 9.967143408044999e-1;
xp(15,1) = 4.000000000000000e+0;
xp(15,2) = 9.999820119661890e-1;
xp(16,1) = 9.919337841451029e-2;
xp(16,2) = 2.528198928766921e-1;
xp(17,1) = 5.076592669645529e-1;
xp(17,2) = 5.550158230159487e-1;
xp(18,1) = 1.184972925827278e+0;
xp(18,2) = 7.852321453615224e-1;
xp(19,1) = 2.047493467134072e+0;
xp(19,2) = 9.245915673876715e-1;
xp(20,1) = 3.007168911869310e+0;
xp(20,2) = 9.839350200445296e-1;
xp(21,1) = 4.000474996776184e+0;
xp(21,2) = 9.984463448413151e-1;
xp(22,1) = 5.000007879022339e+0;
xp(22,2) = 9.999592378464547e-1;
xp(23,1) = 6.000000000000000e+0;
xp(23,2) = 9.999999686258662e-1;
xp(24,1) = 6.001064731474805e-2;
xp(24,2) = 1.538932104518340e-1;
xp(25,1) = 3.149685016229433e-1;
xp(25,2) = 3.551058128559424e-1;
xp(26,1) = 7.664508240518316e-1;
xp(26,2) = 5.449200036280008e-1;
xp(27,1) = 1.396685781342510e+0;
xp(27,2) = 7.104078497715549e-1;
xp(28,1) = 2.175195903206602e+0;
xp(28,2) = 8.398780940253654e-1;
xp(29,1) = 3.062320575880355e+0;
xp(29,2) = 9.272767950890611e-1;
xp(30,1) = 4.016440988792476e+0;
xp(30,2) = 9.750605697371132e-1;
xp(31,1) = 5.002872064275734e+0;
xp(31,2) = 9.942629650823470e-1;
xp(32,1) = 6.000285453310164e+0;
xp(32,2) = 9.992421778421898e-1;
xp(33,1) = 7.000012964962529e+0;
xp(33,2) = 9.999534370786161e-1;
xp(34,1) = 8.000000175554469e+0;
xp(34,2) = 9.999990854912925e-1;
xp(35,1) = 9.000000000000000e+0;
xp(35,2) = 9.999999989466828e-1;

end % nodesRegular



end % methods (Static)








end % classdef
