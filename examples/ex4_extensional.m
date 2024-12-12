clear all; clc

fprintf('Two vesicles in an extensional flow.\n');
fprintf('Second-order semi-implicit time stepping.\n');
fprintf('Implicit vesicle-vesicle interactions.\n');

% Physics parameters
prams.N = 96;                     % points per vesicle
prams.nv = 2;                     % number of vesicles
prams.T = 24;                     % time horizon
prams.m = 4800;                   % number of time steps
prams.kappa = 1e-1;               % bending coefficient
prams.gmresTol = 1e-10;
%prams.gmresMaxIter = 400;

options.farField = 'extensional'; % background velocity
options.order = 1;                % time stepping order
options.vesves = 'implicit';
% Discretization of vesicle-vesicle interactions.
% Either 'explicit' or 'implicit'
options.near = true;              % near-singular integration
%options.nearStrat = 'cauchy';
options.nearStrat = 'interp';
options.fmm = false;
options.axis = [-5 5 -4 4];       % Axis for plots if usePlot=true
options.logFile = 'output/ex4_extensional.log';
% Name of log file for saving messages
options.dataFile = 'output/ex4_extensionalData.bin';
% Name of binary data file for storing vesicle information
options.orderGL = 4;
options.nsdc = 2;

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
X = oc.initConfig(prams.N,'nv',prams.nv,...
    'reducedArea',0.65,...
    'angle',pi/2*ones(prams.nv,1),...
    'center',[[-1.2;0] [1.2;0]],...
    'scale',0.5);
% initial configuration

Xfinal = Ves2D(X,[],prams,options);
% run vesicle simulation

