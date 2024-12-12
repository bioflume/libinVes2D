clear all; clc

fprintf('Multiple vesicles in a double-couette apparatus.\n');
fprintf('First-order semi-implicit time stepping.\n');
fprintf('Explicit vesicle-boundary interactions.\n');

prams.N = 64;                      % points per vesicle
prams.nv = 2;
prams.T = 20;                     % time horizon
prams.m = 2000;                   % number of time steps
prams.Nbd = 64;                    % number of points on solid wall
prams.Nbdcoarse = prams.Nbd;
prams.nvbd = 3;                    % number of components to solid walls
prams.kappa = 1e-1;                % bending coefficient
prams.errorTol = 1e-1;
% Maximum relative error in area and length before the simulation
% is stopped

options.farField = 'doubleCouette'; % Constricted domain
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
options.logFile = 'output/ex7_doubleCouette.log';
% Name of log file for saving messages
options.dataFile = 'output/ex7_doubleCouetteData.bin';
% Name of binary data file for storing vesicle information
options.timeAdap = false;
% use time adaptivity

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
Xwalls = oc.initConfig(prams.Nbd,...
    options.farField,'nv',prams.nvbd,'center',[[0;0] [0;-9] [0;9]]);
% Build solid walls

[X,prams.nv] = oc.initConfig(prams.N,'volFrac',0.10,...
    Xwalls,options.fmm);
prams.viscCont = ones(prams.nv,1);
% Randomly place a template vesicle until the desired volume
% fraction is achieved.  Uses collision detection to decide if they
% a randomly placed vesicle intersects a solid wall or another
% vesicle

Xfinal = Ves2D(X,Xwalls,prams,options);
% Run vesicle code


