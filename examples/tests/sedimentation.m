function [Xfinal] = sedimentation(N,nv,m,Nbd,viscCont,vesves,orderGL,nsdc,saveData,resolveCol,gCont,minSep)
clc
warning off;

prams.N = N;                  % points per vesicle
prams.nv = nv;                  % number of vesicles
prams.T = 10;                  % time horizon
prams.m = m;                % number of time steps
prams.Nbd = Nbd;               % number of points on solid wall
prams.Nbdcoarse = prams.Nbd;
prams.nvbd = 1;                % number of components to solid walls
prams.kappa = 1e-1;            % bending coefficient
prams.viscCont = viscCont;            % viscosity contrast
prams.gmresTol = 1e-12;        % GMRES tolerance
prams.gmresMaxIter = 1000;
prams.errorTol = 2e-1;
prams.minSep = minSep;
prams.rtolArea = 1e-3;
% maxiumum allowable error in area per time step
prams.rtolLength = 1e-3;
% maxiumum allowable error in length per time step

options.farField = 'box';      % Constricted domain
% background velocity or solid walls (which comes with a default
% boundary condition)
options.order = 1;          % time stepping order
options.inextens = 'method1'; 
% Inextensibility condition can be treated as either
% method1 - div x^{N+1} = div x^{N}
% method2 - div(SLP*tractionjump) = 0;
options.vesves = vesves;
% Discretization of vesicle-vesicle and vesicle-boundary 
% intearctions.  Either 'explicit' or 'implicit'
options.near = true;        % near-singular integration
%options.nearStrat = 'cauchy';
options.nearStrat = 'interp';
options.fmm = false;        
options.fmmDLP = false;
% use FMM to compute single-layer potential
options.confined = true;   % confined or unbounded geometry
% Used for geometry and boundary conditins of solid walls
options.usePlot = true;     % View vesicles during run
options.track = false;      % Tracker points on vesicle
options.quiver = false;      % draw velocity vectors on vesicle
options.axis = [-2 2 -3 3]; 
% Axis for plots if usePlot = true
options.verbose = true;
% Decides how much info is written to console
options.profile = false;     % Profile code
options.collision = false;   % Collision detection
options.timeAdap = false;
options.orderGL = orderGL;
options.nsdc = nsdc;
options.saveData = saveData;  
options.resolveCol = resolveCol;
options.gravity = true;
prams.gCont = gCont;
% Save vesicle information and create a log file
filename = ['flow' options.farField 'N' num2str(prams.N) 'nv' num2str(prams.nv) 'Nbd' num2str(prams.Nbd) 'ts' num2str(prams.T/prams.m) 'visc' num2str(prams.viscCont) 'order' ...
num2str(options.order) options.vesves options.nearStrat 'GLorder' num2str(options.orderGL) 'nsdc' num2str(options.nsdc) 'resolveCol' num2str(resolveCol) 'minSep' num2str(prams.minSep) 'gCont' num2str(gCont) 'utest'];
% Name of binary data file for storing vesicle information
options.logFile = ['output/' filename '.log'];
options.dataFile = ['output/' filename '.bin'];

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
if nv == 3
  ang = [pi/2;pi/3;0]; 
  cen = [-0.8,0.7,0;-1.3,-1.2,0.3];
elseif nv == 4
  ang = [pi/2;pi/3;0;pi/3]; 
  cen = [-0.8,0.7,0,0.3;-1.6,-1.5,-0.1,1.5];
elseif nv == 5
  ang = [pi/2;pi/3;0;pi/3;pi/2*0.9]; 
  cen = [-0.8,0.7,0,-0.62,0.76;-1.6,-1.5,-0.1,1.25,1.45];
end

oc = curve;

X = oc.initConfig(prams.N,'nv',prams.nv,'angle',ang,'scale',0.35,...
    'center',cen,'reducedArea',0.90);
[X,~,~] = oc.redistributeParameterize(X,X*0,X(1:prams.N,:)*0);

Xwalls = oc.initConfig(prams.Nbd,options.farField);
[Xwalls,~,~] = oc.redistributeParameterize(Xwalls,Xwalls*0,Xwalls(1:prams.Nbd,:)*0);

if(strcmp(vesves,'implicit'))
  Xfinal = Ves2D(X,Xwalls,[],prams,options);
else
  Xfinal = Ves2DCol(X,Xwalls,[],prams,options);
end
end
