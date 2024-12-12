clear all; clc

disp('One elliptical vesicles in a constricted tube.');
disp('First-order semi-explicit time stepping.');
disp('Explicit vesicle-vesicle interactions.');
disp('Explicit vesicle-boundary interactions.');

prams.N = 128;                 % points per vesicle
prams.nv = 1;                  % number of vesicles
prams.T = 200*1e-2;                  % time horizon
prams.m = 200;                % number of time steps
prams.Nbd = 128;               % number of points on solid wall
prams.Nbdcoarse = prams.Nbd;
prams.nvbd = 1;                % number of components to solid walls
prams.kappa = 1e-1;            % bending coefficient
prams.viscCont = 10;            % viscosity contrast
prams.gmresTol = 1e-10;        % GMRES tolerance
prams.errorTol = 1e-1;
% Maximum relative error in area and length before the simulation
% is stopped
prams.rtolArea = 1e-3;
% maxiumum allowable error in area per time step
prams.rtolLength = 1e-3;
% maxiumum allowable error in length per time step


options.farField = 'choke';      % Constricted domain
% background velocity or solid walls (which comes with a default
% boundary condition)
options.order = 1;          % time stepping order
options.inextens = 'method1'; 
% Inextensibility condition can be treated as either
% method1 - div x^{N+1} = div x^{N}
% method2 - div(SLP*tractionjump) = 0;
options.vesves = 'explicit';
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
options.axis = [-10.5 10.5 -3.5 3.5]; 
% Axis for plots if usePlot = true
options.saveData = false;    
% Save vesicle information and create a log file
options.logFile = 'output/choke1Ves.log';
% Name of log file
options.dataFile = 'output/choke1VesData.bin';
% Name of binary data file for storing vesicle information
options.verbose = true;
% Decides how much info is written to console
options.profile = false;     % Profile code
options.collision = false;   % Collision detection
options.timeAdap = false;
options.expectedOrder = 1;
options.orderGL = 4;
options.nsdc = 2;

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
%X = oc.initConfig(prams.N,'nv',prams.nv,'angle',pi/2,'scale',0.5,...
%    'center',[-9;0],'reducedArea',0.65);
% Initial vesicles has reduced area 0.65 and is a vertical ellipse
% so that it must distort to pass through the contriction

%X = oc.initConfig(prams.N,'nv',prams.nv,...
%    'angle',pi/2*ones(prams.nv,1),'scale',0.48,...
%    'center',[[-3.0;1.5] [-6;0]],'reducedArea',0.9);
X = oc.initConfig(prams.N,'nv',prams.nv,'angle',pi/2,'scale',0.48,...
    'center',[-6.0;1.5],'reducedArea',0.9);

Xwalls = oc.initConfig(prams.Nbd,options.farField);

Xfinal = Ves2D(X,Xwalls,prams,options);
% Run vesicle code


