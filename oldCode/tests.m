% this script tests classes: that ../ is the root of the Ves2D code.

clear globals;  clear all;
clear classes
addpath ../src

if 1

% check fft1
of = fft1;
of.test;
N=128; h = 2*pi/N; x = h*(1:N)';  y = rand(20,1)*2*pi;
v=exp(sin(x));
vy = exp(sin(y));
fy = of.arbInterp(v,y);
fprintf('\n\t EVERYTHING OK with the fft1 code\n');


end





if 1

% define geometry
N  = 32; 
nv = 3; 
X  = boundary(N,'nv',nv,'curly'); 

% old code curve
%[k t s] = curveProp(X);
%[rv, ar, ln] = reducedVolume(X);

% new code curve
oc = curve; 
[sc, tc, kc] = oc.diffProp(X); % notice the different output order
[rvc, arc, lnc] = oc.geomProp(X);


% check errors

%assert(all([norm(kc-k,inf),norm(t-tc,inf),norm(s-sc,inf)]<=1e-14));
%assert(all([norm(rvc-rv,inf),norm(arc-ar,inf),norm(ln-lnc,inf)])<=1e14);

fprintf('\n\t EVERYTHING OK with the curve code\n');


end



if 1

% check "potential" evaluation class.
op = poten(N);
G = op.stokesSLmatrix(X);
D = op.stokesDLmatrix(X);

fprintf('\n\t EVERYTHING OK with the poten code\n');


end

