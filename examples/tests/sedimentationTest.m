%function [Xfinal] = sedimentation(N,nv,m,Nbd,viscCont,vesves,orderGL,nsdc,saveData,resolveCol,gCont,minSep)
warning off;

prams.N = N;                  % points per vesicle
prams.nv = nv;                  % number of vesicles
prams.T = 30;                  % time horizon
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
options.fmm = true;        
options.fmmDLP = true;
% use FMM to compute single-layer potential
options.confined = true;   % confined or unbounded geometry
% Used for geometry and boundary conditins of solid walls
options.usePlot = false;     % View vesicles during run
options.track = false;      % Tracker points on vesicle
options.quiver = false;      % draw velocity vectors on vesicle
options.axis = [-2 2 -3.3 3.3]; 
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
num2str(options.order) options.vesves options.nearStrat 'GLorder' num2str(options.orderGL) 'nsdc' num2str(options.nsdc) 'resolveCol' num2str(resolveCol) 'minSep' num2str(prams.minSep) 'gCont' num2str(gCont)];
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
elseif nv == 100
  xcen =[2.3298,2.2031,2.2470,2.3912,2.2193,2.2779,2.2738,2.2979,2.3805,2.2480,2.2480,2.3026,2.3706,1.3159,1.2291,1.0821,1.1520,1.0737,1.2009,0.9885,1.0077,1.0827,0.9901,1.0308,0.9331,1.0589,-1.4976,-0.0740,0.0062,0.0136,0.0514,0.0099,-0.0297,-0.0498,0.0859,-0.0491,0.1299,-1.0081,-1.0988,-1.1552,-1.2081,-1.1675,-1.1448,-1.0707,-1.1626,-1.1237,-1.0100,-1.0353,-2.3446,-2.3658,-2.2514,-2.3132,-2.2689,-2.2169,-2.2399,-2.2059,-2.2085,-2.1805,-2.3746,-2.3000,1.2000,-2.1000,-0.0000,-1.0500,-0.2000,0.5000,-0.1000,-0.2000,1.2000,-2.3000,-1.3000,-1.2300,-0.1000,-0.8500,1.2000,2.0000,-1.6000,-1.4500,-1.0250,0.3000,2.2000,2.2000,2.2900,2.0900,1.2000,0,0,1.0100,1.0100,1.5000,0.3000,-1.1000,-2.4600,-2.2500,-2.2500,-1.1100,-0.7100,-1.7000,-0.9000,-2.2700];

  ycen =[5.1238,4.0069,2.8248,1.7188,0.5246,-0.6158,-1.7938,-2.9838,-4.1329,-5.2378,-6.3562,-7.4638,-8.5989,-7.7134,-8.9003,-6.0267,-4.1581,-3.0579,-1.9752,-0.3371,0.7847,1.8589,3.7599,4.6818,5.7188,7.6903,7.4281,7.7975,6.6315,4.5854,2.8659,1.0517,-0.4630,-2.0897,-3.9931,-5.8496,-7.1089,-9.0703,-8.0612,-6.4163,-4.7015,-3.6999,-2.4550,-1.3097,0.1447,1.2495,2.7108,4.2937,5.1965,4.1284,2.7839,1.5416,0.3907,-0.6498,-2.2280,-3.7426,-4.9850,-7.3382,-8.4058,-9.5000,-5.1000,-6.0000,-4.9500,-5.5000,-3.1000,-1.2500,2.1000,3.6500,2.7000,6.5000,5.2000,6.0000,5.7000,6.9750,6.8000,6.2000,2.1000,3.6000,-7.2500,-8.0000,6.9700,7.8400,8.7400,9.7400,8.4500,8.6800,9.6300,9.2000,9.9500,10.6400,10.6000,8.2100,7.5000,8.3100,9.0400,9.0000,9.8500,9.8300,10.8000,10.4000];


  ang =[6.1683,1.8411,3.7898,5.1849,2.4000,3.0894,2.2695,3.6927,3.6708,4.1926,2.3456,4.0954,5.4793,5.0257,2.4257,2.7134,6.5978,2.1013,6.4397,2.2409,7.8297,2.5282,2.0975,5.9013,5.3532,3.5262,3.6890,5.1448,6.4664,2.0474,6.3066,7.3337,6.7910,5.4415,3.7698,6.6871,3.1888,5.0091,5.9632,6.2911,3.7095,2.3185,5.8415,6.0279,4.6481,6.3807,1.7872,5.9365,1.8609,6.0070,4.0352,5.8354,1.7952,6.5484,2.4623,7.5849,7.6334,3.3207,7.3097,6.6899,2.6704,3.4558,0.7854,3.1416,3.2987,3.6652,3.1416,3.2987,2.3562,2.3562,3.1416,2.8274,2.5133,0.7854,4.0841,2.3562,-0.5236,3.1416,3.1416,3.1416,-0.7854,-0.7854,-1.0472,-0.8378,0,0,1.4137,0.1000,0.1200,0.1200,0.1200,0.1200,-1.0472,-0.2000,-0.2000,-0.2800,1.5708,1.4137,0,0.9425];

  cen = [xcen;ycen];
  ang = ang';
end

oc = curve;

%X = oc.initConfig(prams.N,'nv',prams.nv,'angle',ang,'scale',0.35,...
%    'center',cen,'reducedArea',0.90);

X = oc.initConfig(prams.N,'nv',prams.nv,'angle',ang,'scale',0.21,...
    'center',cen,'reducedArea',0.90);
[X,~,~] = oc.redistributeParameterize(X,X*0,X(1:prams.N,:)*0);

Xwalls = oc.initConfig(prams.Nbd,options.farField);
[Xwalls,~,~] = oc.redistributeParameterize(Xwalls,Xwalls*0,Xwalls(1:prams.Nbd,:)*0);

if(strcmp(vesves,'implicit'))
  Xfinal = Ves2D(X,Xwalls,[],prams,options);
else
  Xfinal = Ves2DCol(X,Xwalls,[],prams,options);
end
