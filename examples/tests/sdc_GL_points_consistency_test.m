clear all; clc

disp('One elliptical vesicles in a constricted tube, consistency test.');

prams.N = 128;                 % points per vesicle
prams.nv = 1;                  % number of vesicles
prams.T = 1e-2;                  % time horizon
prams.m = 1;                % number of time steps
prams.Nbd = 128;               % number of points on solid wall
prams.Nbdcoarse = prams.Nbd;
prams.nvbd = 1;                % number of components to solid walls
prams.kappa = 1e-1;            % bending coefficient
prams.viscCont = 2;            % viscosity contrast
prams.gmresTol = 1e-12;        % GMRES tolerance
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
options.verbose = false;
% Decides how much info is written to console
options.profile = false;     % Profile code
options.collision = false;   % Collision detection
options.timeAdap = false;
options.expectedOrder = 1;
options.orderGL = 5;
options.nsdc = 0;

addpath('../examples');
[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
X = oc.initConfig(prams.N,'nv',prams.nv,'angle',pi/2,'scale',0.48,...
    'center',[-6.0;1.5],'reducedArea',0.9);

Xwalls = oc.initConfig(prams.Nbd,options.farField);


GLpts = [2,3,4,5,6,7];
for i = 1:numel(GLpts)
  prams.m = 1;
  options.orderGL = GLpts(i);
  fprintf([ num2str(GLpts(i)) ' GL pts.\n']);
  Xfinal = Ves2D(X,Xwalls,prams,options);
  %XfinalCol = Ves2DCol(X,Xwalls,prams,options);
  XfinalGLtest = Ves2DGLtest(X,Xwalls,prams,options);
  
 
  %fprintf([ 'norm of position diff between orginal code and factorized code: ' num2str(norm(Xfinal-XfinalCol)) '.\n']);
  fprintf([ 'norm of position diff between orginal code with GL steps and orginal code with sdc GL time steps: ' num2str(norm(Xfinal - XfinalGLtest)) '.\n']);
  
end


