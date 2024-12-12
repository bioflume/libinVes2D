clear all; clc

fprintf('Multiple vesicles in a couette apparatus.\n');
fprintf('First-order adaptive semi-implicit time stepping.\n');
fprintf('Explicit vesicle-boundary interactions.\n');

prams.N = 92;                  % points per vesicle
prams.T = 50;                  % time horizon
prams.m = 200;                  
% number of time steps.  Will be changed by adpativity
prams.Nbd = 128;               % number of points on solid wall
prams.Nbdcoarse = prams.Nbd; 
prams.nvbd = 2;                % number of components to solid walls
prams.kappa = 1e-1;            % bending coefficient
prams.errorTol = 1e-1;
% Maximum relative error in area and length before the simulation
% is stopped
prams.rtolArea = 1e-1;
% allowable total error in area
prams.rtolLength = 1e-1;
% allowable total error in length
prams.maxdt = 1e0;
% maximum time step size


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
options.fmm = false;
options.fmmDLP = false;
% use FMM to compute single-layer potential
options.confined = true;   % confined or unbounded geometry
% Used for geometry and boundary conditins of solid walls
options.axis = [-21 21 -21 21];
% Axis for plots if usePlot = true
options.logFile = 'output/ex6_couette.log';
% Name of log file for saving messages
options.dataFile = 'output/ex6_couetteData.bin';
% Name of binary data file for storing vesicle information
options.timeAdap = true;
% use time adaptivity
options.orderGL = 4;
% order of internal time steps
options.nsdc = 1;
% number of SDC corrections

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
Xwalls = oc.initConfig(prams.Nbd,...
    options.farField,'nv',prams.nvbd,'center',[[0;0] [-5;0]]);
% Build solid walls

prams.nv = 1;
X = oc.initConfig(prams.N,'nv',prams.nv,'center',[12;0]);
% set up initial vesicle configuration

Xfinal = Ves2D(X,Xwalls,prams,options);
% Run vesicle code


