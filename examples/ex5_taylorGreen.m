clear all; clc

fprintf('Nine vesicles in a Taylor-Green flow.\n');
fprintf('Second-order semi-implicit time stepping.\n');
fprintf('Implicit vesicle-vesicle interactions.\n');

% Physics parameters
prams.N = 64;                     % points per vesicle
prams.nv = 9;                     % number of vesicles
prams.T = 10;                     % time horizon
prams.m = 200;                    % number of time steps
prams.kappa = 1e-1;               % bending coefficient
prams.viscCont = 1;               % viscosity contrast

options.farField = 'taylorGreen'; % background velocity
options.inextens = 'method1';
% method of enforcing inextensibility.
% Can be 'method1' or 'method2'
options.order = 1;                % time stepping order
options.vesves = 'implicit';
% Discretization of vesicle-vesicle intearctions.
% Either 'explicit' or 'implicit'
options.near = true;        % near-singular integration
options.fmm = false;        
% use FMM to compute single-layer potential
options.axis = [0 pi 0 pi] + [-0.1 0.1 -0.1 0.1]; 
% axis for plot
options.logFile = 'output/ex5_taylorGreen.log';
% Name of log file for saving messages
options.dataFile = 'output/ex5_taylorGreenData.bin';
% Name of binary data file for storing vesicle information

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path

oc = curve;
sx = [.02 -.01 -.014 .12 .04 -.08 .08 .02 -.06]; 
sy = [-.02 -.03 -.011 .04 .02 .02 .03 .06 -.02]; 
cenx = kron(ones(1,3),[pi/4 pi/2 3*pi/4]);
ceny = kron([pi/4 pi/2 3*pi/4],ones(1,3));
cenx = cenx + sx;
ceny = ceny + sy;
% x and y coordinates of centers
% Have added pertubation to centers so that flow is 
% more interesting
angle = ones(prams.nv,1);
X = oc.initConfig(prams.N,'nv',prams.nv,...
    'reducedArea',0.65,...
    'angle',angle,...
    'center',[cenx;ceny],...
    'scale',0.15);
% Initial configuration

Ves2D(X,[],prams,options);
% Run vesicle code

