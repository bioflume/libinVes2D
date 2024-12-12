%N =16;m =3000;viscCont = 10;vesves = 'explicit';orderGL=3;nsdc=2;saveData=1;resolveCol=1;minSep=0.25;Nbd = 1024;
warning off;

prams.N = N;                  % points per vesicle
prams.nv = 80;                  % number of vesicles
prams.T = 30;                  % time horizon
prams.m = m;                % number of time steps
prams.Nbd = Nbd;               % number of points on solid wall
prams.Nbdcoarse = prams.Nbd;
prams.nvbd = 1;                % number of components to solid walls
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

options.farField = 'chokeerf';      % Constricted domain
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
options.saveData = saveData;    
options.resolveCol = resolveCol;
options.restart = false;
% options.restartFile = '';

filename = ['flow' options.farField 'N' num2str(prams.N) 'nv' num2str(prams.nv) 'Nbd' num2str(prams.Nbd) 'ts' num2str(prams.T/prams.m) 'visc' num2str(prams.viscCont) 'order' ...
num2str(options.order) options.vesves options.nearStrat 'GLorder' num2str(options.orderGL) 'nsdc' num2str(options.nsdc) 'resolveCol' num2str(resolveCol) 'minSep' num2str(prams.minSep) 'fmm' num2str(options.fmm) 'near' num2str(options.near)];
% Name of binary data file for storing vesicle information
options.logFile = ['output/' filename '.log'];
options.dataFile = ['output/' filename '.bin'];

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

cenx =[-4.8762,-5.9931,-7.1752,-8.2812,-9.4754,-10.6158,-11.7938,-12.9838,-14.1329,-15.2378,-16.3562,-17.4638,-18.5989,-17.7134,-18.9003,-16.0267,-14.1581,-13.0579,-11.9752,-10.3371,-9.2153,-8.1411,-6.2401,-5.3182,-4.2812,-2.3097,-2.5719,-2.2025,-3.3685,-5.4146,-7.1341,-8.9483,-10.4630,-12.0897,-13.9931,-15.8496,-17.1089,-19.0703,-18.0612,-16.4163,-14.7015,-13.6999,-12.4550,-11.3097,-9.8553,-8.7505,-7.2892,-5.7063,-4.8035,-5.8716,-7.2161,-8.4584,-9.6093,-10.6498,-12.2280,-13.7426,-14.9850,-17.3382,-18.4058,-19.5000,-15.1000,-16.0000,-14.9500,-15.5000,-13.1000,-11.2500,-7.9000,-6.3500,-7.3000,-3.5000,-4.8000,-4.0000,-4.3000,-3.0250,-3.2000,-3.8000,-7.9000,-6.4000,-17.2500,-18.0000];

ceny =[-2.3298,-2.2031,-2.2470,-2.3912,-2.2193,-2.2779,-2.2738,-2.2979,-2.3805,-2.2480,-2.2480,-2.3026,-2.3706,-1.3159,-1.2291,-1.0821,-1.1520,-1.0737,-1.2009,-0.9885,-1.0077,-1.0827,-0.9901,-1.0308,-0.9331,-1.0589,1.4976,0.0740,-0.0062,-0.0136,-0.0514,-0.0099,0.0297,0.0498,-0.0859,0.0491,-0.1299,1.0081,1.0988,1.1552,1.2081,1.1675,1.1448,1.0707,1.1626,1.1237,1.0100,1.0353,2.3446,2.3658,2.2514,2.3132,2.2689,2.2169,2.2399,2.2059,2.2085,2.1805,2.3746,2.3000,-1.2000,2.1000,0,1.0500,0.2000,-0.5000,0.1000,0.2000,-1.2000,2.3000,1.3000,1.2300,0.1000,0.8500,-1.2000,-2.0000,1.6000,1.4500,1.0250,-0.3000];

ang =[4.5975,0.2703,2.2190,3.6141,0.8292,1.5186,0.6987,2.1220,2.1000,2.6218,0.7748,2.5246,3.9085,3.4549,0.8549,1.1426,5.0270,0.5305,4.8689,0.6701,6.2589,0.9574,0.5267,4.3305,3.7824,1.9554,2.1182,3.5740,4.8957,0.4766,4.7358,5.7629,5.2203,3.8707,2.1990,5.1163,1.6180,3.4383,4.3924,4.7204,2.1387,0.7477,4.2707,4.4571,3.0773,4.8099,0.2164,4.3657,0.2901,4.4362,2.4644,4.2646,0.2244,4.9776,0.8915,6.0141,6.0626,1.7499,5.7389,5.1191,1.0996,1.8850,-0.7854,1.5708,1.7279,2.0944,1.5708,1.7279,0.7854,0.7854,1.5708,1.2566,0.9425,-0.7854,2.5133,0.7854,-2.0944,1.5708,1.5708,1.5708];

cen = [cenx;ceny];
ang = ang';

oc = curve;
% Initial vesicles has reduced area 0.65 and is a vertical ellipse
% so that it must distort to pass through the contriction
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
