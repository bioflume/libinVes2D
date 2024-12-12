function [Xfinal] = chokeTraffic6nv(N,m,Nbd,viscCont,vesves,orderGL,nsdc)
%clear all;
clc

prams.N = N;                  % points per vesicle
prams.nv = 6;                  % number of vesicles
prams.T = 26;                  % time horizon
prams.m = m;                % number of time steps
prams.Nbd = Nbd;               % number of points on solid wall
prams.Nbdcoarse = prams.Nbd;
prams.nvbd = 1;                % number of components to solid walls
prams.kappa = 1e-1;            % bending coefficient
prams.viscCont = viscCont;            % viscosity contrast
prams.gmresTol = 1e-12;        % GMRES tolerance
prams.gmresMaxIter = 100;
prams.errorTol = 2e-1;
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
options.axis = [-10.5 10.5 -3.5 3.5]; 
% Axis for plots if usePlot = true
options.verbose = true;
% Decides how much info is written to console
options.profile = false;     % Profile code
options.collision = false;   % Collision detection
options.timeAdap = false;
options.orderGL = orderGL;
options.nsdc = nsdc;
options.saveData = true;    
% Save vesicle information and create a log file
filename = ['flow' options.farField 'N' num2str(prams.N) 'nv' num2str(prams.nv) 'Nbd' num2str(prams.Nbd) 'ts' num2str(prams.T/prams.m) 'visc' num2str(prams.viscCont) 'order' ...
num2str(options.order) options.vesves options.nearStrat 'GLorder' num2str(options.orderGL) 'nsdc' num2str(options.nsdc)];
% Name of binary data file for storing vesicle information
options.logFile = ['output/' filename '.log'];
options.dataFile = ['output/' filename '.bin'];

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path
cen = [-3,-4.4,-2.6,-4.0,-2.5,-4.5;1.8,1.78,-1.6,-1.76,0.1,0];
ang = [pi/2;pi/2*0.9;pi*0.8/2;-pi*0.92/2;pi*0.9;pi]; 

oc = curve;

% Initial vesicles has reduced area 0.65 and is a vertical ellipse
% so that it must distort to pass through the contriction
X = oc.initConfig(prams.N,'nv',prams.nv,'angle',ang,'scale',0.35,...
    'center',cen,'reducedArea',0.90);

Xwalls = oc.initConfig(prams.Nbd,options.farField);
disp([num2str(prams.nv) ' elliptical vesicles in a constricted tube.']);
disp('First-order semi-implicit time stepping.');
disp([options.vesves ' vesicle-vesicle interactions.']);
disp([options.vesves ' vesicle-boundary interactions.']);
Xfinal = Ves2D(X,Xwalls,prams,options);
% Run vesicle code

end
