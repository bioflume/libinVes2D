clear all; clc
warning off;
%fprintf('Single vesi cle in a shear flow.\n');
%fprintf('First-order semi-implicit time stepping.\n');

% Physics parameters
prams.N =64;               % points per vesicle
prams.nv = 1;               % number of vesicles
prams.Nrd =64;
prams.nvrd = 1;
prams.T = 20;               % time horizon
prams.m = 500;              % number of time steps
prams.kappa = 1e-1;         % bending coefficient
prams.gmresTol = 1e-12;     %
prams.gmresMaxIter = 100;
prams.viscCont = 100;        % viscosity contrast
prams.errorTol = 2e-1;
prams.minSep = 1;
options.farField = 'extensional'; % background velocity
options.order = 1;          % time stepping order
options.vesves = 'explicit';
options.inextens = 'method1';
options.withRigid = true;
 
options.near= true;
options.nearStrat = 'interp';
options.fmm = false;        
options.fmmDLP = false;

options.saveData = false; 
options.logFile = 'output/shear2Ves.log';
% Name of log file for saving messages
options.dataFile = 'output/shear2Ves.bin';
% Name of binary data file for storing vesicle information
options.orderGL = 2;
options.nsdc = 0;
options.resolveCol = true;
[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path
oc = curve;
angle = pi/2*ones(prams.nv,1);
scale = 0.5;
X = oc.initConfig(prams.N,'nv',prams.nv,...
    'scale',scale,...
    'reducedArea',0.97,...
    'angle',angle,...
    'center',[ -2; 0]);
% Initial configuration of reduce area 0.65 and aligned


Xrig = oc.initConfig(prams.Nrd,'nv',prams.nvrd,...
    'scale',scale,...
    'reducedArea',0.97,...
    'angle',pi/2*ones(prams.nvrd,1),...
    'center',[2 ; 0]);
Xrig = reshape(Xrig,prams.Nrd,[]);
Xrig = Xrig(end:-1:1,:);
Xrig = reshape(Xrig,2*prams.Nrd,prams.nvrd);

[Xf1] = Ves2DCol(X,[],Xrig,prams,options);