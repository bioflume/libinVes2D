function [Xfinal] = taylorrgd(N,m,viscCont,vesves,orderGL,nsdc,saveData,resolveCol,minSep)
clc
warning off;

prams.N = 0;                  % points per vesicle
prams.nv = 0;                  % number of vesicles
prams.Nrd = N;
prams.nvrd = 25;
prams.T = 50;                  % time horizon
prams.m = m;                % number of time steps
prams.Nbd = 0;               % number of points on solid wall
prams.Nbdcoarse = prams.Nbd;
prams.nvbd = 0;                % number of components to solid walls
prams.kappa = 1e-1;            % bending coefficient
prams.viscCont = viscCont;            % viscosity contrast
prams.gmresTol = 1e-12;        % GMRES tolerance
prams.gmresMaxIter = 500;
prams.errorTol = 2e-1;
prams.minSep = minSep;
prams.rtolArea = 1e-3;
% maxiumum allowable error in area per time step
prams.rtolLength = 1e-3;
% maxiumum allowable error in length per time step

options.farField = 'taylorGreen';      % Constricted domain
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
options.confined = false;   % confined or unbounded geometry
% Used for geometry and boundary conditins of solid walls
options.usePlot = true;     % View vesicles during run
options.track = false;      % Tracker points on vesicle
options.quiver = false;      % draw velocity vectors on vesicle
options.axis = [-15 15 -15 15]; 
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
options.withRigid = 1;

filename = ['flow' options.farField 'N' num2str(prams.N) 'nv' num2str(prams.nv) 'ts' num2str(prams.T/prams.m) 'visc' num2str(prams.viscCont) 'order' ...
num2str(options.order) options.vesves options.nearStrat 'GLorder' num2str(options.orderGL) 'nsdc' num2str(options.nsdc) 'resolveCol' num2str(resolveCol) 'minSep' num2str(prams.minSep)];
% Name of binary data file for storing vesicle information
options.logFile = ['output/' filename '.log'];
options.dataFile = ['output/' filename '.bin'];

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path
cen = [-6 -3 0 3 6 -6.1 -3.2 0.2 3.2 6.3 -6 -3 0 3 6 -6.1 -3 0 3 6 -6 -3 0 3 6 ;...
       -2 -2.1 -2.2 -2.5 -2.1 2.3 2.2 2.1 2.3 2.3 0.0 0.0 -0.0 -0.0 0.0 4.1 4.2 4.3 4.4 4.5 -4.1 -4.2 -4.1 -4.0 -3.9];
ang = [0.2;0.3;0.5;0.3;0.25;0.34;0.56;0.78;0.32;0.34;0.2;0.3;0.5;0.3;0.25;0.34;0.56;0.78;0.32;0.34;1;3;2;3;3]; 
oc = curve;

% Initial vesicles has reduced area 0.65 and is a vertical ellipse
% so that it must distort to pass through the contriction
Xrig = oc.initConfig(prams.Nrd,'nv',prams.nvrd,'angle',ang,'scale',0.3,...
    'center',cen,'reducedArea',0.65);
Xrig = reshape(Xrig,prams.Nrd,[]);
Xrig = Xrig(end:-1:1,:);
Xrig = reshape(Xrig,2*prams.Nrd,prams.nvrd);
%[Xrig,~,~] = oc.redistributeParameterize(Xrig,Xrig*0,Xrig(1:prams.Nrd,:)*0);


if(strcmp(vesves,'implicit'))
  Xfinal = Ves2D([],[],Xrig,prams,options);
else
  Xfinal = Ves2DCol([],[],Xrig,prams,options);
end
end
