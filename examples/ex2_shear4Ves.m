clear all; clc

fprintf('Four vesicles in a shear flow.\n');
fprintf('Second-order semi-implicit time stepping.\n');
fprintf('Implicit vesicle-vesicle interactions.\n');

% Physics parameters
prams.N = 64;               % points per vesicle
%prams.N = 64;               % points per vesicle
prams.nv = 4;               % number of vesicles
prams.T = 1e-5;               % time horizon
prams.m = 1;              % number of time steps
prams.kappa = 1e-1;         % bending coefficient
prams.gmresTol = 1e-10;
prams.gmresMaxIter = 100;
prams.viscCont = 4*ones(1,prams.nv);

options.farField = 'shear'; % background velocity
options.inextens = 'method1';
% method1 - div x^{N+1} = div x^{N}
% method2 - div(SLP*tractionjump) = 0;
options.order = 1;          % time stepping order
options.vesves = 'implicit';
% Discretization of vesicle-vesicle interactions.
% Either 'explicit' or 'implicit'
options.near = true;        % near-singular integration
options.fmm = false;
options.logFile = 'output/ex2_shear4Ves.log';
% Name of log file for saving messages
options.dataFile = 'output/ex2_shear4VesData.bin';
% Name of binary data file for storing vesicle information
options.axis = [-6 6 -5 5];
% Axis for the plot
options.usePlot = 1;
options.saveData = 1;
options.orderGL = 3;
options.nsdc = 1;
options.adhesion = false;

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
angle = pi/2*ones(prams.nv,1);
scale = 0.4;
X = oc.initConfig(prams.N,'nv',prams.nv,...
    'scale',scale,...
    'reducedArea',0.65,...
    'angle',angle,...
    'center',[[-3 -1 1 3];[0 0 0 0]]);
% Initial configuration


%tic;
options.nearStrat = 'interp';
Xfinal1 = Ves2D(X,[],prams,options);
%oldtime=toc;

%tic;
%options.nearStrat = 'cauchy';
%Xfinal2 = Ves2D(X,[],prams,options);
%newtime=toc;

%fprintf('Timing old code: %g\n',oldtime);
%fprintf('Timing new code: %g\n',newtime);



