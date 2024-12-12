function [order] = relaxation_sdc_conv_order_test(method, nsdc)


disp('one vesicle relaxation convergence test with sdc corrections.');

prams.N = 96;                 % points per vesicle
prams.nv = 1;                  % number of vesicles
prams.T = 0.1;                  % time horizon
prams.m = 10;                % number of time steps
prams.kappa = 1e-1;            % bending coefficient
prams.viscCont = 1;            % viscosity contrast
prams.gmresTol = 1e-14;        % GMRES tolerance
prams.errorTol = 1e-1;
% Maximum relative error in area and length before the simulation
% is stopped
prams.rtolArea = 1e-3;
% maxiumum allowable error in area per time step
prams.rtolLength = 1e-3;
% maxiumum allowable error in length per time step


options.farField = 'relaxation';      % Constricted domain
% background velocity or solid walls (which comes with a default
% boundary condition)
options.order = 1;          % time stepping order
options.inextens = 'method1'; 
% Inextensibility condition can be treated as either
% method1 - div x^{N+1} = div x^{N}
% method2 - div(SLP*tractionjump) = 0;
options.vesves = method;
% Discretization of vesicle-vesicle and vesicle-boundary 
% intearctions.  Either 'explicit' or 'implicit'
options.near = true;        % near-singular integration
%options.nearStrat = 'cauchy';
options.nearStrat = 'interp';
options.fmm = false;        
options.fmmDLP = false;
% use FMM to compute single-layer potential
% Used for geometry and boundary conditins of solid walls
options.usePlot = true;     % View vesicles during run
options.track = false;      % Tracker points on vesicle
options.quiver = false;      % draw velocity vectors on vesicle
options.axis = [-3 3 -4 4]; 
% Axis for plots if usePlot = true
options.saveData = false;    
% Save vesicle information and create a log file
options.logFile = 'relaxation_conv.log';
% Name of log file
options.dataFile = 'relaxation_conv.bin';
% Name of binary data file for storing vesicle information
options.verbose = false;
% Decides how much info is written to console
options.profile = false;     % Profile code
options.collision = false;   % Collision detection
options.timeAdap = false;
%options.expectedOrder = 1;
options.orderGL = 5;
options.nsdc = 0;

addpath('../examples');
[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
X = oc.initConfig(prams.N,'curly');
[ra0,a0,l0] = oc.geomProp(X);
sampleSteps = 4;

options.nsdc = nsdc;
prams.m = 10;
dt = [];
lengtherror = [];
for m = 1:sampleSteps
  dt = [dt;prams.T/prams.m];
  Xfinal = Ves2D(X,[],prams,options);
  [ra1,a1,l1] = oc.geomProp(Xfinal);
  lengtherror = [lengtherror;max(abs(l1-l0)/l0)];
  fprintf(['dt: ' num2str(prams.T/prams.m) ',\t area error: ' num2str(abs(a0-a1)/abs(a0)) '    ---\tlength error: ' num2str(abs(l0-l1)/abs(l0)) '\n']);
  prams.m = prams.m*2;
end
a = fitlm(log(dt),log(lengtherror));
order =  a.Coefficients.Estimate(2);
fprintf(['Estimated convergence order is: ' num2str(order) '\n']);
end
