function [options,prams] = initVes2D(options,prams)
% Set a path pointing to src directory and set options and
% prams to default values if they are left blank

P = path; ii = find(pwd == filesep); ii = ii(end);
subPath = pwd; subPath = [subPath(1:ii) 'src'];
if isempty(strfind(P, subPath)),addpath(subPath);end
subPath=[subPath '/lsc2d_code'];
if isempty(strfind(P, subPath)),addpath(subPath);end

PramList = {'N','nv','T','m','Nbd','nvbd','kappa','viscCont',...
    'gmresTol','gmresMaxIter','errorTol',...
    'rtolArea','rtolLength',...
    'betaUp','betaDown','alpha','adStrength','adRange'};
defaultPram.N = 64;
defaultPram.nv = 1;
defaultPram.Nbd = 0;
defaultPram.nvbd = 0;
defaultPram.T = 1;
defaultPram.m = 100;
defaultPram.kappa = 1e-1;
defaultPram.viscCont = 1;
defaultPram.gmresTol = 1e-12;
defaultPram.gmresMaxIter = 50;
defaultPram.errorTol = 1e-1;
defaultPram.rtolArea = 1e-2;
defaultPram.rtolLength = 1e-2;
defaultPram.betaUp = 1.5;
defaultPram.betaDown = 0.6;
defaultPram.alpha = sqrt(0.9);
defaultPram.adRange = 8e-1;
defaultPram.adStrength = 4e-1;

for k = 1:length(PramList)
  if ~isfield(prams,PramList{k})
    eval(['prams.' PramList{k} '=defaultPram.' PramList{k} ';'])
    % Set any unassigned parameters to a default value
  end
end

OptList = {'order','inextens','farField','vesves','near',...
    'nearStrat','fmm','fmmDLP','confined','usePlot','track',...
    'quiver','axis','saveData','logFile','dataFile','verbose',...
    'profile','collision','timeAdap','bending','tracers',...
    'SDCcorrect','pressure','orderGL','nsdc','adhesion',...
    'correctShape'};

defaultOpt.order = 1;
defaultOpt.inextens = 'method1';
defaultOpt.farField = 'shear';
defaultOpt.vesves = 'implicit';
defaultOpt.near = true;
defaultOpt.nearStrat = 'interp';
defaultOpt.fmm = false;
defaultOpt.fmmDLP = false;
defaultOpt.confined = false;
defaultOpt.usePlot = true;
defaultOpt.track = false;
defaultOpt.quiver = false;
defaultOpt.axis = [-5 5 -5 5];
defaultOpt.saveData = true;
defaultOpt.logFile = 'output/example.log';
defaultOpt.dataFile = 'output/exampleData.bin';
defaultOpt.verbose = true;
defaultOpt.profile = false;
defaultOpt.collision = false;
defaultOpt.timeAdap = false;
defaultOpt.bending = false;
defaultOpt.tracers = false;
defaultOpt.SDCcorrect = false;
defaultOpt.pressure = false;
defaultOpt.orderGL = 2;
defaultOpt.nsdc = 0;
defaultOpt.adhesion = false;
defaultOpt.correctShape = false;


for k = 1:length(OptList)
  if ~isfield(options,OptList{k})
    eval(['options.' OptList{k} '=defaultOpt.' OptList{k} ';'])
    % Set any unassigned options to a default value
  end
end

if ~options.confined
  prams.Nbd = 0;
  prams.nvbd = 0;
end
% If the geometry is unbounded, make sure to set the number
% of points and number of components of the solid walls
% to 0.  Otherwise, later components will crash

if numel(prams.viscCont) ~=prams.nv
  prams.viscCont = prams.viscCont*ones(1,prams.nv);
end

if options.nsdc > 0
  if options.order > 1
    fprintf('***************************************************\n')
    fprintf('Can only do sdc updates with first-order\n');
    fprintf('Setting sdc corrections to zero\n');
    fprintf('PUSH ANY KEY TO CONTINUE\n');
    pause
    fprintf('***************************************************\n')
    options.nsdc = 0;
  end

  if (any(prams.viscCont ~= 1) && strcmp(options.vesves,'explicit') && options.nsdc > 0)
    fprintf('***************************************************\n')
    fprintf('Not sure if this combination works\n');
    fprintf('See how rhs1 is built at SDC corrections\n');
    fprintf('PUSH ANY KEY TO CONTINUE\n');
    pause
    fprintf('***************************************************\n')
%    options.nsdc = 0;
  end

end





