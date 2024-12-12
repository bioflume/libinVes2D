clear all; clc



prams.N = 64;                  % points per vesicle
prams.T = 10;                  % time horizon
prams.m = 1000;                  
% number of time steps.  Will be changed by adpativity
prams.Nbd = 256;               % number of points on solid wall
prams.Nbdcoarse = prams.Nbd; 
prams.nvbd = 2;                % number of components to solid walls
prams.viscCont = 1;
prams.gmresTol = 1e-12;        % GMRES tolerance
prams.gmresMaxIter = 1000;
prams.kappa = 1e-1;            % bending coefficient
prams.errorTol = 2e-1;
prams.minSep = 1.0;
% Maximum relative error in area and length before the simulation
% is stopped
prams.rtolArea = 1e-1;
% allowable total error in area
prams.rtolLength = 1e-1;
% allowable total error in length
prams.gCont = 0;

options.farField = 'couette';  % Constricted domain
% background velocity or solid walls (which comes with a default
% boundary condition)
options.order = 1;          % time stepping order
options.inextens = 'method1'; 
% Inextensibility condition can be treated as either
% method1 - div x^{N+1} = div x^{N}
% method2 - div(SLP*tractionjump) = 0;
options.vesves = 'implicit';
% Discretization of vesicle-vesicle intearctions.
% Either 'explicit' or 'implicit'
options.near = true;        % near-singular integration
options.nearStrat = 'interp';
options.fmm = true;
options.fmmDLP = true;
% use FMM to compute single-layer potential
options.confined = true;   % confined or unbounded geometry
% Used for geometry and boundary conditins of solid walls
options.axis = [-21 21 -21 21];
options.usePlot = false;     % View vesicles during run
% Axis for plots if usePlot = true
options.logFile = 'output/couette196VesRigImplcit1000m.log';
% Name of log file for saving messages
options.dataFile = 'output/couette196VesRigImplcit1000m.bin';
% Name of binary data file for storing vesicle information
options.timeAdap = false;
% use time adaptivity
options.profile = false;     % Profile code
options.collision = false;   % Collision detection
options.orderGL = 2;
% order of internal time steps
options.nsdc = 0;
% number of SDC corrections
options.saveData = true;  
options.resolveCol = false;
options.gravity = false;
options.withRigid = true;

rng('default');

nvrd1 = 5;
nvrd2 = 7;
nvrd3 = 9;

idx1 = [1;10;21;31;49];
idx2 = [3;8;27;34;45;49;64];
idx3 = [5;17;24;38;42;55;61;69;76];


nv1 = 50-nvrd1;nv2 = 68-nvrd2;nv3 = 78-nvrd3;
prams.nv = nv1+nv2+nv3;
prams.Nrd = 64;
prams.nvrd = nvrd1+nvrd2+nvrd3;
[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
Xwalls = oc.initConfig(prams.Nbd,...
    options.farField,'nv',prams.nvbd,'center',[[0;0] [0;0]]);
[Xwalls,~,~] = oc.redistributeParameterize(Xwalls,Xwalls*0,Xwalls(1:prams.Nbd,:)*0);
% Build solid walls
%[x,y] = oc.getXY(Xwalls);
%plot(x,y,'-x')
%hold on;


cen1 = [11.85;0];
cen = [];
angle = [];
cenrd = [];
anglerd = [];
for i = 1:(nv1+nvrd1)
  alpha = i*2*pi/(nv1+nvrd1);
  R = [cos(alpha) -sin(alpha);sin(alpha) cos(alpha)];
  rig = false;
  for j = 1:nvrd1
    if i==idx1(j)
      rig = true;
    end
  end
  if ~rig
    angle = [angle;alpha+(rand-0.5)*pi/5];
    cen = [cen R*cen1];
  else
    anglerd = [anglerd;alpha+(rand-0.5)*pi/5];
    cenrd = [cenrd R*cen1];
  end
end

cen2 = [15;0];
for i = 1:(nv2+nvrd2)
  alpha = i*2*pi/(nv2+nvrd2);
  R = [cos(alpha) -sin(alpha);sin(alpha) cos(alpha)];
  rig = false;
  for j = 1:nvrd2
    if i==idx2(j)
      rig = true;
    end
  end
  if ~rig
    angle = [angle;alpha+(rand-0.5)*pi/5];
    cen = [cen R*cen2];
  else
    anglerd = [anglerd;alpha+(rand-0.5)*pi/5];
    cenrd = [cenrd R*cen2];
  end
end


cen3 = [18.15;0];
for i = 1:(nv3+nvrd3)
  alpha = i*2*pi/(nv3+nvrd3);
  R = [cos(alpha) -sin(alpha);sin(alpha) cos(alpha)];
  rig = false;
  for j = 1:nvrd3
    if i==idx3(j)
      rig = true;
    end
  end
  if ~rig
    angle = [angle;alpha+(rand-0.5)*pi/5];
    cen = [cen R*cen3];
  else
    anglerd = [anglerd;alpha+(rand-0.5)*pi/5];
    cenrd = [cenrd R*cen3];
  end
end


X = oc.initConfig(prams.N,'nv',prams.nv,'center',cen,'angle',angle,'reducedArea',0.65,'scale',0.43);
% set up initial vesicle configuration
[X,~,~] = oc.redistributeParameterize(X,X*0,X(1:prams.N,:)*0);


Xrig = oc.initConfig(prams.Nrd,'nv',prams.nvrd,'center',cenrd,'angle',anglerd,'reducedArea',0.65,'scale',0.43);
% set up initial vesicle configuration
[Xrig,~,~] = oc.redistributeParameterize(Xrig,Xrig*0,Xrig(1:prams.N,:)*0);
Xrig = reshape(Xrig,prams.Nrd,[]);
Xrig = Xrig(end:-1:1,:);
Xrig = reshape(Xrig,2*prams.Nrd,prams.nvrd);

%{
[x,y] = oc.getXY(X);
rigid = capsules(X,[],[],[],[],[]);
rigid.normal
quiver(x,y,rigid.normal(1:prams.Nrd,:),rigid.normal(prams.Nrd+1:end,:))
hold on
plot(x,y,'-xr');
axis equal
hold on

[x,y] = oc.getXY(Xrig);
rigid = capsules(Xrig,[],[],[],[],[]);
rigid.normal
plot(x,y,'-xb');
quiver(x,y,rigid.normal(1:prams.Nrd,:),rigid.normal(prams.Nrd+1:end,:))
axis equal
%}

Xfinal = Ves2D(X,Xwalls,Xrig,prams,options);
% Run vesicle code


