clear all; clc

fprintf('Single vesicle in a shear flow.\n');
fprintf('First-order semi-implicit time stepping.\n');

% Physics parameters
prams.N = 32;               % points per vesicle
prams.nv = 1;               % number of vesicles
prams.T = 10;               % time horizon
prams.m = 200;              % number of time steps
prams.kappa = 1e-1;         % bending coefficient
prams.viscCont = 1;         % viscosity contrast

options.farField = 'shear'; % background velocity
% Save vesicle information and create a log file
options.logFile = 'output/ex1_shear1Ves.log';
% Name of log file for saving messages
options.dataFile = 'output/ex1_shear1VesData.bin';
% Name of binary data file for storing vesicle information
options.correctShape = true;

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
X = oc.initConfig(prams.N,'reducedArea',0.65,'angle',pi/2);
% Initial configuration of reduce area 0.65 and aligned

Ves2D(X,[],prams,options);
% Run vesicle code

