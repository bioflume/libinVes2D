clear all; clc

fprintf('Simple curly vesicle in a relaxation flow.\n');
fprintf('Second-order semi-implicit time stepping.\n');

% Physics parameters
prams.N = 32;                    % points per vesicle
prams.nv = 1;                    % number of vesicles
prams.T = 1;                     % time horizon
prams.m = 100;                   % number of time steps
prams.kappa = 1e-1;              % bending coefficient
prams.viscCont = 1;              % viscosity contrast

options.farField = 'relaxation'; % background velocity
options.order = 1;               % time stepping order
options.axis = [-3 3 -4 4];      % Axis for plot
options.logFile = 'output/ex3_relaxation.log';
% Name of log file for saving messages
options.dataFile = 'output/ex3_relaxationData.bin';
% Name of binary data file for storing vesicle information
%options.profile = true;

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
X = oc.initConfig(prams.N,'curly');
% Initial configuration
%theta = (0:prams.N-1)'*2*pi/prams.N;
%X = [cos(theta);3*sin(theta)];

Xfinal = Ves2D(X,[],prams,options);
% Run vesicle code

