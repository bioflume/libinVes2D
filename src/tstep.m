classdef tstep < handle
% This class defines the functions required to advance the geometry
% forward in time.  Handles both implicit and explicit vesicle-vesicle
% interactions, different inextensibility conditions, viscosity
% contrast, solid walls vs. unbounded flows.  This class also
% implements the adaptive time stepping strategy where the errors in
% length, area, and the residual are monitored to adaptively choose a
% time step size.


properties
Xcoeff        % Explicit discretization
rhsCoeff      % Explicit terms from the discretization of 
              % the derivative
beta          % Term multiplying implicit term in discretization of
              % derivative
order         % Time stepping order
expectedOrder % order that is used to pick new time step size
dt            % Time step size
currentTime   % current time needed for adaptive time stepping
finalTime     % time horizon
solver        % method1, method2, method3, or method4.  Changes how
              % inextensiblity condition is treated
Galpert       % Single-layer stokes potential matrix using Alpert
Gbarnett      % layer potential in equation 4.12 of Alex and 
              % Shravan's paper
D             % Double-layer stokes potential matrix
lapDLP        % Double-layer laplace potential matrix
wallRestrict  % Rate that the preconditioner is restricted for
              % the solid wall components
wallDLP       % Double-layer potential due to solid walls
wallN0        % Modificiation of double-layer potential to 
              % remove the rank deficiency on outer boundary
farField      % Far field boundary condition
confined      % whether geometry is bounded or not
bdiagVes      % precomputed inverse of block-diagonal precondtioner
              % only for vesicle-vesicle interactions
bdiagTen
bdiagWall     % precomputed inverse of block-diagonal precondtioner
              % only for wall-wall interactions
bdiagWallCoarse
bdiagRigid
vesves        % Discretization of vesicle-vesicle interactions.
              % Can be implicit or explicit
fmm           % with or without the FMM
fmmDLP        % use FMM for the double-layer potential
near          % with or without near-singular integration
nearStrat     % strategy for doing near-singular integration
profile       % display profile times throughout simulation
bending       % bending is defined so that the vesicle tends to
              % a prescribed configuration
gmresTol      % GMRES tolerance
gmresMaxIter  % maximum number of gmres iterations
collDete      % decides if we do collision detection
timeAdap      % Using adaptive time stepping
orderGL       % Order of Gauss-Lobatto points
GLpts         % Gauss-Lobatto points 
SDCcorrect    % timeStep changes depending on if we are doing
              % an sdc correction or forming a provisional
              % solution
rtolArea      % maximum allowable error in area
rtolLength    % maximum allowable error in length
nsdc          % number of sdc iterations to take
NearV2V       % near-singular integration for vesicle to vesicle
NearW2V       % near-singular integration for wall to vesicle
NearV2W       % near-singular integration for vesicle to wall 
NearW2W       % near-singular integration for wall to wall 
NearV2R
NearW2R
NearR2W
NearR2R
NearR2V
op            % class for evaluating potential so that we don't have
              % to keep building quadrature matricies
betaUp        % maximum amount time step is scaled up
betaDown      % maximum amount time step is scaled down
alpha         % safety factor scaling for new time step size
adhesion      % use adhesion in the model
periodic      % fudges in periodicity into the geometry
adStrength    % strength of adhesion
adRange       % range of adhesion
antiAlias     % use anti-aliasing
resolveCol    % resolve collision
minSep        % minimal seperation distance
gravity
gCont
withRigid
rigidDLP
D1up
upMx
downMx
uprate
end % properties

methods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = tstep(options,prams)
% o=tstep(options,prams: constructor.  Initialize class.  Take all
% elements of options and prams needed by the time stepper

o.order = options.order; % Time stepping order
o.expectedOrder = options.expectedOrder; % Expected time stepping order
o.dt = prams.T/prams.m; % Time step size
[o.Xcoeff,o.rhsCoeff,o.beta] = o.getCoeff(o.order);
% Get coefficients for doing time integration
o.currentTime = 0;
% Method always starts at time 0
o.finalTime = prams.T;
% Need the time horizon for adaptive time stepping
o.solver = options.inextens;
% Discretization of inextensibility
o.vesves = options.vesves;
% Vesicle-vesicle interactions
o.fmm = options.fmm;
o.fmmDLP = options.fmmDLP;
% fast-multipole method
o.near = options.near;
% near-singular integration
o.nearStrat = options.nearStrat;
% strategy for doing near-singular integration.
% either 'cauchy' or 'interp'
o.profile = options.profile;
% disp profile times
o.bending = options.bending;
% user-defined prefered configuration of vesicle
o.gmresTol = prams.gmresTol;
% GMRES tolerance
o.gmresMaxIter = prams.gmresMaxIter;
% maximum number of GMRES iterations
o.farField = @(X) o.bgFlow(X,options.farField,options.farFieldSpeed);
% Far field boundary condition built as a function handle
o.confined = options.confined;
% Confined or unbounded geometry
o.collDete = options.collision;
% Collision detection
o.timeAdap = options.timeAdap;
% Time adaptivity flag
o.orderGL = options.orderGL;
% Gauss-Lobatto order
o.GLpts = o.gaussLobatto(o.orderGL);
% load Gauss-Lobatto points
o.rtolArea = prams.rtolArea/prams.T;
o.rtolLength = prams.rtolLength/prams.T;
o.nsdc = options.nsdc;;
% number of sdc iterations to take
if prams.nv
  o.op = poten(prams.N,options.antiAlias);
  
  o.uprate = 256/prams.N;%4;
  o.D1up = fft1.D1(prams.N*o.uprate);
  o.upMx = o.upSampleMatrix(prams.N,o.uprate);
  o.downMx = o.downSampleMatrix(prams.N,o.uprate);
elseif prams.nvrd
  o.op = poten(prams.Nrd,options.antiAlias);
end
% Build class with quadrature points and weights
% as well as lagrange interpolation matrix
if o.order > 1
  o.SDCcorrect = false;
end

o.betaUp = prams.betaUp;
o.betaDown = prams.betaDown;
o.alpha = prams.alpha;
% safety factors for adaptive time stepping

o.adhesion = options.adhesion;
o.periodic = options.periodic;
o.wallRestrict = 1;

o.adStrength = prams.adStrength;
% strength of adhesion
o.adRange = prams.adRange;
% range of adhesion

o.antiAlias = options.antiAlias;
% anti-aliasing flag

o.resolveCol = options.resolveCol;

o.minSep = prams.minSep;

o.gravity = options.gravity;
o.gCont = prams.gCont;
o.withRigid = options.withRigid;
end % tstep: constructor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xcoeff,rhsCoeff,beta] = getCoeff(o,order)
% [Xcoeff,rhsCoeff,beta] = getCoeff(order) generates the coefficients
% required to discretize the derivative.  First-order time  derivatives
% are approximated by beta*x^{N+1} + rhscoeff.*[x^{N} x^{N-1} ...]
% Explicit terms (operators) are discretized at Xcoeff.*[x^{N} x^{N-1}
% ...] All rules are from Ascher, Ruuth, Wetton 1995.

if (order == 4) % fourth-order
  beta = 25/12;
  Xcoeff = [-1; 4; -6; 4]; 
  rhsCoeff = [-1/4; 4/3; -3; 4];
elseif (order == 3) % third-order
  beta = 11/6;
  Xcoeff = [1; -3; 3];
  rhsCoeff = [1/3; -3/2; 3];
elseif (order == 2) % second-order
  beta = 1.5;
  Xcoeff = [-1; 2];
  rhsCoeff = [-0.5; 2];
else % first-order
  beta = 1;
  Xcoeff = 1;
  rhsCoeff = 1;
end


end % getCoeff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [walls,wallsCoarse] = initialConfined(o,prams,Xwalls)
% walls = initialConfined(prams,Xwalls) builds an object of class
% capsules for the solid walls

op = poten(prams.Nbd);
uwalls = o.farField(Xwalls);
% velocity on solid walls coming from no-slip boundary condition
walls = capsules(Xwalls,[],uwalls,zeros(prams.nvbd,1),...
    zeros(prams.nvbd,1),o.antiAlias);

XwallsCoarse = curve.restrict(Xwalls,prams.Nbd,prams.Nbdcoarse);
wallsCoarse = capsules(XwallsCoarse,[],[],zeros(prams.nvbd,1),...
    zeros(prams.nvbd,1),o.antiAlias);

o.wallDLP = op.stokesDLmatrix(wallsCoarse);
o.wallN0 = op.stokesN0matrix(wallsCoarse);
%o.wallDLP = op.stokesDLmatrix(walls);
%o.wallN0 = op.stokesN0matrix(walls);
o.bdiagWall = o.wallsPrecond(wallsCoarse);
% block diagonal preconditioner for solid walls build the double-layer
% potential and the matrix N0 that removes the rank one null space from
% the double-layer potential for the solid walls

end % initialConfined


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xstore,sigStore,uStore,etaStore,RSstore,Xtra] = ...
  firstSteps(o,options,prams,Xinit,sigInit,uInit,Xrig,etaRig,...
  walls,wallsCoarse,om,Xtra,pressTar)
% [Xstore,sigStore,uStore,etaStore,RSstore] = ...
%   firstSteps(options,prams,Xinit,sigInit,uInit,...
%   walls,wallsCoarse,om,pressTar)
% refines the first time step [0,dt] and uses a first-order, then
% second-order, then third-order, etc time stepping to find the vesicle
% position, tension, velocity, and density function eta defined on the
% solid walls and the rotlets and stokeslets at t = dt.  Returns
% Xstore, sigStore, uStore etaStore, so that they are immediatly ready
% to use for higher-order time stepping.  This routine also computes
% the tension and density functions of the initial condition.  This is
% needed in SDC since we require these functions at all Gauss-Lobatto
% points in order to compute the residual

N = size(Xinit,1)/2; % points per vesicle
nv = size(Xinit,2); % number of vesicles
if o.confined 
  Xwalls = walls.X; % discretization points of solid walls
  Nbd = size(Xwalls,1)/2; % points per wall
  nvbd = size(Xwalls,2); % number of wall components
else
  Xwalls = [];
  Nbd = 0;
  nvbd = 0;
end

prams.T = prams.T/prams.m*(o.order-1);
% time horizon is enough time steps so that we can
% continue with higher-order time integrator
if o.order ~=1
  mR = min(ceil(prams.m/32),100)*o.order^2;
  mR = mR * (o.order - 1);
else
  mR = 1;
end
% For second-order time stepping, this keeps the local error from time
% step 0 to dt on the order of dt^3
prams.m = mR*(o.order-1);
% number of time steps to take in range [0,dt]

Xstore = zeros(2*N,nv,prams.m+1);
sigStore = zeros(N,nv,prams.m+1);
uStore = zeros(2*N,nv,prams.m+1);
etaStore = zeros(2*Nbd,nvbd,prams.m+1);
RSstore = zeros(3,nvbd,prams.m+1);

Xstore(:,:,1) = Xinit;
vesicle = capsules(Xinit,zeros(N,nv),[],...
    prams.kappa,prams.viscCont,o.antiAlias);


[sigStore(:,:,1),etaStore(:,:,1),RSstore(:,:,1),uStore(:,:,1),iter] = ...
    vesicle.computeSigAndEta(o,walls);
%[sigStore(:,:,1),etaStore(:,:,1),RSstore(:,:,1),u,iter] = ...
%    vesicle.computeSigAndEta(o,walls);
% need intial tension, density function, rotlets, and stokeslets so
% that we can do SDC updates, couling with advection-diffusion system,
% etc

om.initializeFiles(Xinit,sigStore(:,:,1),etaStore(:,:,1),...
    RSstore(:,:,1),Xrig,etaRig,Xwalls,Xtra,pressTar);
% delete previous versions of files and write some initial
% options/parameters to files and the console

message = ['GMRES took ',num2str(iter),...
    ' iterations to find intial tension and density function'];
om.writeMessage(message,'%s\n');
om.writeMessage(' ','%s\n');

for n = 1:prams.m
  time = n*prams.T/prams.m;

  options.order = min(n,o.order);
  tt = tstep(options,prams);
  % take the highester order possible
  tt.SDCcorrect = false;
  % no SDC corrections
  tt.wallDLP = o.wallDLP;
  tt.wallN0 = o.wallN0;
  tt.bdiagWall = o.bdiagWall;
  % build inverse of the wall-to-wall intearctions includes diagonal and
  % off-diagonal terms.  That is, it could be used to solve the stokes
  % equations in a multiply-connected domain with no vesicles

  updatePreco = true;
  [X,sigma,u,eta,RS,iter,iflag] = tt.timeStep(...
      Xstore(:,:,n-tt.order+1:n),...
      sigStore(:,:,n-tt.order+1:n),...
      uStore(:,:,n-tt.order+1:n),...
      etaStore(:,:,n-tt.order+1:n),...
      RSstore(:,:,n-tt.order+1:n),...
      [],[],[],[],[],[],[],...
      prams.kappa,prams.viscCont,walls,wallsCoarse,...
      updatePreco);
  % take one time step

  if numel(Xtra) > 1
    vel = o.tracersVel(X,sigma,u,...
        prams.kappa,prams.viscCont,walls,eta,RS,Xtra);
    Xtra = Xtra + tt.dt*vel;
  end

  accept = true;
  dtScale = 0;
  res = 0;
  % Required for adaptive time stepping which have not been
  % implemented yet for high-order time stepping

  Xstore(:,:,n+1) = X;
  sigStore(:,:,n+1) = sigma;
  uStore(:,:,n+1) = u;
  etaStore(:,:,n+1) = eta;
  RSstore(:,:,n+1) = RS;
  % save and update position, tension, velocity, and
  % density function

  terminate = om.outputInfo(X,sigma,u,eta,RS,...
      Xwalls,Xtra,time,iter,dtScale,res,iflag);
  % save data, write to log file, write to console as
  % requested
end 
% end of using small time steps to get the simulation far enough
% to use the desired order time integrator

if o.order > 1
  if options.pressure
    op = poten(N);
    [press,stress1,stress2] = op.pressAndStress(...
        X,sigma,u,prams.kappa,prams.viscCont,...
        walls,pressTar,eta,RS,options.confined,...
        options.fmm);
    % compute the pressure and stress due to the vesicles
    % and the solid walls
    if options.saveData
      om.writePressure(press);
      om.writeStress(stress1(1:end/2),'11');
      om.writeStress(stress1(1+end/2:end),'12');
      om.writeStress(stress2(1:end/2),'21');
      om.writeStress(stress2(1+end/2:end),'22');
    end
  end % ~options.pressure
  % write pressure and stress to output file at time dt

  Xstore = Xstore(:,:,1:mR:end);
  sigStore = sigStore(:,:,1:mR:end);
  uStore = uStore(:,:,1:mR:end);
  etaStore = etaStore(:,:,1:mR:end);
  RSstore = RSstore(:,:,1:mR:end);
end
% Only want every mR^{th} solution as these correspond to
% time steps that are multiples of dt
% Only need to do this if the time stepping order is not 1

end % firstSteps


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,sigma,u,eta,RS,Xrig,etaRig,iter,accept,dtScale,normRes,iflag] = ...
    timeStepGL(o,Xstore,sigStore,uStore,...
    etaStore,RSstore,XrigStore,etaRigStore,kappa,viscCont,walls,wallsCoarse,om,time,accept)
% [X,sigma,u,eta,RS,iter,dt,accept,dtScale,normRes,iflag] = ...
%    timeStepGL(Xstore,sigStore,uStore,...
%    etaStore,RSstore,kappa,viscCont,walls,wallsCoarse,a0,l0)
% takes the desired number of time steps at Gauss-Lobatto points to
% find solution at dt given solution at 0.  Calls o.timeStep which is
% what we use for the constant sized time step integration.  Returns a
% flag that says if the solution was accepted or not and the amount the
% time step was scaled.  errors is the norm of the residual and iflag
% is tells if any of the gmres runs failed

N = size(Xstore,1)/2; % number of points per vesicle
nv = size(Xstore,2); % number of vesicles
Nbd = size(etaStore,1)/2; % number of points per solid wall
nvbd = size(etaStore,2); % number of solid walls
a0 = om.area; % input area
l0 = om.length; % input length
oc = curve;

Xprov = zeros(2*N,nv,abs(o.orderGL));
sigmaProv = zeros(N,nv,abs(o.orderGL));
uProv = zeros(2*N,nv,abs(o.orderGL));
etaProv = zeros(2*Nbd,nvbd,abs(o.orderGL));
RSprov = zeros(3,nvbd,abs(o.orderGL));
% need to save the position, tension, velocity, density
% functions, rotlets/stokeslets at the Gauss-Lobatto
% time steps
for k = 1:nv
  Xprov(:,k,1) = Xstore(:,k);
  sigmaProv(:,k,1) = sigStore(:,k);
  uProv(:,k,1) = uStore(:,k);
end
% Store the initial conditions in the provisional
% solution which will store the solution at all
% Gauss-Lobatto points

for k = 1:nvbd
  etaProv(:,k,1) = etaStore(:,k);
  RSprov(:,k,1) = RSstore(:,k);
end

Galpert = zeros(2*N,2*N,nv,abs(o.orderGL));
%Gbarnett = zeros(N,N,nv,abs(o.orderGL));
Gbarnett = [];
% need the single-layer potential at all levels
% of the provisional solution
if any(viscCont ~= 1)
  D = zeros(2*N,2*N,nv,abs(o.orderGL));
else
  D = [];
end
% need the double-layer potential at all levels
% of the provisional solution

deltaX = zeros(2*N,nv,abs(o.orderGL));
deltaSigma = zeros(N,nv,abs(o.orderGL));
deltaEta = zeros(2*Nbd,nvbd,abs(o.orderGL));
deltaRS = zeros(3,nvbd,abs(o.orderGL));
% Errors that are solved for in each iteration.  The first column of
% deltaX will always be zero since we assume the inital condition is
% exact

X = Xprov(:,:,1);
sigma = sigmaProv(:,:,1);
u = uProv(:,:,1);

vesicle(1) = capsules(X,sigma,u,kappa,viscCont,o.antiAlias);
% vesicle configuration at first Gauss-Lobatto point

if o.timeAdap
  [~,aInit,lInit] = oc.geomProp(X);
  eaInit = abs(aInit - a0)./a0;
  elInit = abs(lInit - l0)./l0;
end
% compute intial reduced area, area, and length if using adaptive
% time stepping

dt = o.dt;
% need to save the time step size
dtimes = diff(o.GLpts)*dt/2;
% time step sizes use for Gauss-Lobatto points

% START OF FORMING PROVISIONAL SOLUTION
iflag = 0; % initialize gmres flag as everything okay
iter = 0; % initialize number of gmres iterations
updatePreco = true; 
% Form preconditioner at initial configuartion
for k = 1:abs(o.orderGL)-1
  o.dt = dtimes(k);
  o.SDCcorrect = false;

  [X,sigma,u,eta,RS,Xrig,etaRig,subIter,iflagTemp] = o.timeStep(...
      Xprov(:,:,k-o.order+1:k),sigmaProv(:,:,k-o.order+1:k),...
      uProv(:,:,k-o.order+1:k),etaProv(:,:,k-o.order+1:k),...
      RSprov(:,:,k-o.order+1:k),XrigStore,etaRigStore,[],[],[],[],[],[],[],...
      kappa,viscCont,walls,wallsCoarse,updatePreco);
  % form provisional solution at next Gauss-Lobatto point
  iter = iter + subIter;
  % running counter for the total number of gmres iterations
  % per time step.  GMRES is used at each Gauss-lobatto time
  % substep which we need to add up for a fair comparison to
  % the traditional time stepping method

  updatePreco = false;
  % don't want to precompute the preconditioner at the next time
  % step to save CPU time.  Since the configuration doesn't 
  % change by much, the preconditioner shouldn't change
  % by much either and the old preconditioner should
  % work fine

  if iflagTemp ~= 0
    iflag = iflagTemp;
  end
  % if any of the gmres iterations fail, assign the failure to
  % flag to iflag for monitor to output

  if o.nsdc > 0
    Galpert(:,:,:,k) = o.Galpert;
%    Gbarnett(:,:,:,k) = o.Gbarnett;
    if any(viscCont ~= 1)
      D(:,:,:,k) = o.D;
    end
    % want to save single-layer potential for computing the residual.
    % This will not save the last one, but this is computed in
    % computeResidual
    NearV2V{k} = o.NearV2V;
    NearW2V{k} = o.NearW2V;
    NearV2W{k} = o.NearV2W;
    NearW2W{k} = o.NearW2W;

    vesicle(k+1) = capsules(X,sigma,u,kappa,viscCont,o.antiAlias);
    % need to save if doing SDC corrections 
  end

  Xprov(:,:,k+1) = X;
  sigmaProv(:,:,k+1) = sigma;
  uProv(:,:,k+1) = u;
  etaProv(:,:,k+1) = eta;
  RSprov(:,:,k+1) = RS;
  % save the provisional solution at the intermediate 
  % Gauss-Lobatto point

  if o.resolveCol
    
    Nrd = size(XrigStore,1)/2; % number of points per vesicle
    nvrd = size(XrigStore,2); % number of vesicles
    if o.withRigid
      X0 = [Xprov(:,:,k-o.order+1:k) XrigStore(:,:,1)];
    else
      X0 = [Xprov(:,:,k-o.order+1:k)];
    end
    X1 = [X Xrig];
    oc = curve;
    cellSize = 999;
    if nv
      Nmax = N;
      [ra,area,length] = oc.geomProp(X0(:,1:nv));
      edgelength = length/N;
      cellSize = min(cellSize,max(edgelength));
    end
    if o.withRigid
      Nmax = Nrd;
      [ra,area,length] = oc.geomProp(X0(:,nv+1:end));
      edgelengthRig = length/Nrd;
      cellSize = min(cellSize,max(edgelengthRig));
    end
    if o.confined
      [ra,area,length] = oc.geomProp(walls.X);
      walllength = max(length/Nbd);
      wallupSamp = ceil(walllength/cellSize);
    else
      wallupSamp = 1;
    end
    upSampleFactor = 1;
    nexten = 0;
    c_tol = 1e-12;
    minSep = 0;
    maxIter = 1000;

    [Ns,totalnv,Xstart,Xend,totalPts] = o.preColCheck(X0,X1,walls,wallupSamp);
    [vgrad, iv, ids, vols] = getCollision(Ns, totalnv, Xstart, Xend, minSep, maxIter, totalPts, c_tol, nv+nvrd, nvbd, Nmax*upSampleFactor, Nbd*wallupSamp, nexten,cellSize*0.7);
    om.writeMessage(['ivolume: ' num2str(iv)],'%s\n');

    if(iv<0)
      om.writeMessage('collision, exit matlab.','%s\n');
      exit;
    end

  end
end % k


% need to save the near-singular integration structure and
% layer-potential matricies at the final Gauss-Lobatto point for future
% SDC sweeps and computing the residual, if there are any SDC sweeps
if o.nsdc > 0
  if o.near
    if o.confined
      [NearV2V{abs(o.orderGL)},NearV2W{abs(o.orderGL)}] = ...
          vesicle(abs(o.orderGL)).getZone(walls,3);
      % Need vesicle to vesicle and vesicle to wall interactions
      if nvbd == 1 
        [NearW2W{abs(o.orderGL)},NearW2V{abs(o.orderGL)}] = ...
            walls.getZone(vesicle(abs(o.orderGL)),2);
      else
        [NearW2W{abs(o.orderGL)},NearW2V{abs(o.orderGL)}] = ...
            walls.getZone(vesicle(abs(o.orderGL)),3);
      end
      % Only need wall to vesicle interactions.  Wall to wall 
      % interactions should also use near-singular integration since
      % they may be close to one another
    else
      NearV2V{abs(o.orderGL)} = ...
          vesicle(abs(o.orderGL)).getZone([],1);
      % no solid walls, so only need vesicle-vesicle intearactions
      NearV2W{abs(o.orderGL)} = [];
      NearW2V{abs(o.orderGL)} = [];
      NearW2W{abs(o.orderGL)} = [];
    end
  else
    NearV2V{abs(o.orderGL)} = [];
    NearV2W{abs(o.orderGL)} = [];
    NearW2V{abs(o.orderGL)} = [];
    NearW2W{abs(o.orderGL)} = [];
  end
  % Need near-singular integration strucutures for final state
  % for sdc updates and computing the residual
  Galpert(:,:,:,abs(o.orderGL)) = ...
      o.op.stokesSLmatrix(vesicle(o.orderGL));
%  Gbarnett(:,:,:,abs(o.orderGL)) = ...
%      o.op.laplaceSLcomplexMatrix(vesicle(o.orderGL));
  if any(viscCont ~= 1)
    D(:,:,:,abs(o.orderGL)) = ...
        o.op.stokesDLmatrix(vesicle(o.orderGL));
  end
  % need single-layer potential at the final state of the 
  % provisional solution to form the residual
end
% END OF FORMING PROVISIONAL SOLUTION AND THE SINGLE-LAYER
% POTENTIAL AT THESE INTERMEDIATE STATES

o.dt = dt;
% change back to original time step
 
%color = ['r' 'g' 'b' 'k' 'c' 'm' 'y'];
%figure(2); clf; hold on
%for k = 1:abs(o.orderGL)
%  plot(sigmaProv(:,1,k),color(k))
%end
%figure(3); clf; hold on
%for k = 1:abs(o.orderGL)
%  plot(squeeze(etaProv(1:end/2,1,k)),color(k))
%end
%figure(5);clf;hold on;
%for k = 1:abs(o.orderGL)
%  plot(etaProv(end/2+1:end,1,k),color(k))
%end
%pause
%% DEBUG: TO MAKE SURE THAT THE TENSION IS CONTINUOUS FROM BETWEEN
%% THE DIFFERENT GAUSS-LOBATTO POINTS

Xinit = Xprov;
uInit = uProv;
sigmaInit = sigmaProv;
etaInit = etaProv;
RSinit = RSprov;
% save the very first provisional solution so that it can be taken
% as the solution if sdc corrections increase the error but this
% original provisional solution has the smallest error
if o.nsdc > 0
  [vesVel,divVesVel,wallVel,residual] = ...
      o.computeVelocities(vesicle,Galpert,Gbarnett,D,walls,...
        etaProv,RSprov,NearV2V,NearV2W,NearW2V,NearW2W,uProv);
  % form the residual as well as the velocities due to the different
  % components which are necessary to form the right-hand sides when
  % doing sdc corrections

  normRes = max(abs(residual(:,:,end)));
  normRes = max(normRes);
  % save the size of the residual at the final Gauss-Lobatto
  % point.  We may have more of these later if we do a full
  % deferred correction method.  Use the maximum L2 error where
  % the maximum is taken over all the vesicles
else
  normRes = 0;
end
% compute the integrand that is in the residual and the residual 
% of the picard integral formulation of the dynamic equation
% If the residual is not desired, then errors is set to 0 and
% this will never refine the time step size

[~,a,l] = oc.geomProp(Xprov(:,:,end));
ea = abs(a - a0)./abs(a0);
el = abs(l - l0)./abs(l0);
eaVec = ea;
elVec = el;
message = ['sdcCount   ' num2str(0,'%2d') ...
  ': Residual = ' num2str(normRes(end),'%5.2e') ...
  ', eA = ' num2str(max(ea),'%5.2e') ...
  ', eL = ' num2str(max(el),'%5.2e')];
om.writeMessage(message,'%s\n');
% print the provisional solution's residual and errors

% Start doing SDC corrections
for sdcCount = 1:o.nsdc
uDeltaX = 0*deltaX(:,:,1);
  for n = 1:abs(o.orderGL) - 1
    o.dt = dtimes(n);

    updatePreco = false;
    o.SDCcorrect = true;
    o.Galpert = Galpert(:,:,:,n+1);
%    o.Gbarnett = Gbarnett(:,:,:,n+1);
    if any(viscCont ~= 1)
      o.D = D(:,:,:,n+1);
    end
    o.NearV2V = NearV2V{n+1};
    o.NearV2W = NearV2W{n+1};
    o.NearW2V = NearW2V{n+1};
    o.NearW2W = NearW2W{n+1};
    [X,sigma,u,eta,RS,tmp1,tmp2,subIter,iflagTemp] = o.timeStep(...
        Xprov(:,:,n+1),sigmaProv(:,:,n+1),...
        uDeltaX,etaProv(:,:,n+1),RSprov(:,:,n+1),[],[],...
        deltaX(:,:,n),deltaSigma(:,:,n),...
        deltaEta(:,:,n),deltaRS(:,:,n),...
        residual(:,:,n+1) - residual(:,:,n),...
        vesVel(:,:,n+1),wallVel(:,:,n+1),...
        kappa,viscCont,walls,wallsCoarse,updatePreco,...
        vesicle(1).sa,vesicle(1).IK,Xprov(:,:,n));
    % Form the sdc update
    iter = iter + subIter;

    if iflagTemp ~= 0
      iflag = iflagTemp;
    end
    o.SDCcorrect = false;
    % turn correct off since it is not used when forming the
    % provisional solution
    deltaX(:,:,n+1) = X;
    deltaSigma(:,:,n+1) = sigma;
    deltaEta(:,:,n+1) = eta;
    deltaRS(:,:,n+1) = RS;
    % approximations of the error
    uDeltaX = (deltaX(:,:,n+1) - deltaX(:,:,n))/dtimes(n);
    uProv(:,:,n+1) = uProv(:,:,n+1) + uDeltaX;
      
    if o.resolveCol  

      X0 = Xprov(:,:,n+1);
      X1 = Xprov(:,:,n+1)+X;
      oc = curve;
      cellSize = 999;
      if nv
	Nmax = N;
	[ra,area,length] = oc.geomProp(X0(:,1:nv));
	edgelength = length/N;
	cellSize = min(cellSize,max(edgelength));
      end
      if o.confined
	[ra,area,length] = oc.geomProp(walls.X);
	walllength = max(length/Nbd);
	wallupSamp = ceil(walllength/cellSize);
      else
	wallupSamp = 1;
      end
      upSampleFactor = 1;
      nexten = 0;
      c_tol = 1e-12;
      minSep = 0;
      maxIter = 1000;

      [Ns,totalnv,Xstart,Xend,totalPts] = o.preColCheck(X0,X1,walls,wallupSamp);
      [vgrad, iv, ids, vols] = getCollision(Ns, totalnv, Xstart, Xend, minSep, maxIter, totalPts, c_tol, nv, nvbd, Nmax*upSampleFactor, Nbd*wallupSamp, nexten,cellSize*0.7);
      
      if(iv<0)
	om.writeMessage('collision, exit matlab.','%s\n');
	exit;
      end

    end
  end
  o.dt = dt;
  % go back to original time step
  
  Xprov = Xprov + deltaX;
  sigmaProv = sigmaProv + deltaSigma;
  etaProv = etaProv + deltaEta;
  RSprov = RSprov + deltaRS;
  % update provision solution velocities
  % update provision solution
  
  if sdcCount < o.nsdc
    for n = 1:abs(o.orderGL)
      vesicle(n) = capsules(Xprov(:,:,n),sigmaProv(:,:,n),...
          [],kappa,viscCont,o.antiAlias);
      Galpert(:,:,:,n) = o.op.stokesSLmatrix(vesicle(n));
%      Gbarnett(:,:,:,n) = o.op.laplaceSLcomplexMatrix(vesicle(n));
      if any(viscCont ~= 1)
        D(:,:,:,n) = o.op.stokesDLmatrix(vesicle(n));
      else
        D = [];
      end
    end
    % update the vesicle objects and the single-layer potentials

    [vesVel,divVesVel,wallVel,residual] = ...
      o.computeVelocities(vesicle,Galpert,Gbarnett,D,walls,...
          etaProv,RSprov,NearV2V,NearV2W,NearW2V,NearW2W,uProv);
  end
  % if we are only recording the error in area and length, don't need
  % to compute the final residual

  normRes = max(abs(residual(:,:,end)));
  normRes = max(normRes);

  [~,a,l] = oc.geomProp(Xprov(:,:,end));
  ea = abs(a - a0)./abs(a0);
  el = abs(l - l0)./abs(l0);
  eaVec = [eaVec ea];
  elVec = [elVec el];
  message = ['sdcCount   ' num2str(sdcCount,'%2d') ...
    ': Residual = ' num2str(normRes(end),'%5.2e') ...
    ', eA = ' num2str(max(ea),'%5.2e') ...
    ', eL = ' num2str(max(el),'%5.2e')];
  om.writeMessage(message,'%s\n');
end
% End of doing SDC corrections
% update solution with an SDC iteration

if o.timeAdap
  [accept,dtScale] = o.newTimeStepSize(a,l,...
      aInit,lInit,accept,om);
  % if doing adaptive time stepping, get new time step
  % size and time step scaling
else
  accept = true;
  dtScale = 1;
  % if not doing adaptive time stepping, keep the same
  % time step size and always accept solution
end

if accept
  % take the new solution
  X = Xprov(:,:,end);
  sigma = sigmaProv(:,:,end);
  if o.nsdc > 0
    u = vesVel(:,:,end);
  else
    u = uProv(:,:,end);
  end
  eta = etaProv(:,:,end);
  RS = RSprov(:,:,end);
else
  % revert to the old solution
  X = Xstore;
  sigma = sigStore;
  u = uStore;
  eta = etaStore;
  RS = RSstore;
end

if accept && o.periodic
  X = oc.addAndRemove(X,walls,o.near,o.fmm);
end

end % timeStepGL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [accept,dtScale] = newTimeStepSize(o,...
      aFinal,lFinal,aInit,lInit,accept,om)
% [accept,dtScale] = newTimeStepSize(...
%      aFinal,lFinal,aInit,lInit,accept,om)
% finds a new time step size based on the change in the area and length
% from the first to final Gauss- Lobatto points.  The output is a flag
% indicating acceptance or rejection, and the amount the time step is
% scaled.  To increase the likelihood that the a new time step size is
% accepted, the time step size is never increased if the previous time
% step was rejected

alpha = o.alpha;
% buffer so that we aren't trying to keep error exactly at 1.  This is a
% safeguard so that the next time step is accepted with higher
% probability
betaUp = o.betaUp;
% allowable upscaling change in time step size
betaDown = o.betaDown;
% allowable downscaling change in time step size

errArea = abs(aInit - aFinal);
errLength = abs(lInit - lFinal);
% absolute errors in area and length
tauArea = max(aInit,aFinal)*o.rtolArea*o.dt;
tauLength = max(lInit,lFinal)*o.rtolLength*o.dt;
% Tolerance for errArea and errLength
err = max(max(errArea./tauArea),max(errLength./tauLength));
% Maximum of relative errors in area and length
% Want this quantity to be as close to 1 as possible

actualOrder = o.expectedOrder;
% order of time stepping method
dtOPT = err^(-1/(actualOrder))*o.dt;
% optimal time step size

dtOld = o.dt;
if accept
  o.dt = alpha^(1/actualOrder) * ...
      min(betaUp*o.dt,max(dtOPT,betaDown*o.dt));
else
  o.dt = alpha^(1/(actualOrder)) * ...
      min(o.dt,max(dtOPT,betaDown*o.dt));
  % don't want to scale up this time step if it was previously rejected.
  % This hopefully gets rid of pattern where a solution alternates
  % between being accepted and rejected.
end
% safety factor added to the optimal time step size also, time step size
% is not scaled up or down too fast For safety factor, take 1/p root.
% In our formulation, this makes alpha the desired value for err
% regardless of the order of the method.
dtScale = o.dt/dtOld;
% time step scaling


if err > 1
  accept = false;
  % reject time step because the error is too large
  message = ['Time Step REJECTED with error ' ...
      num2str(err,'%4.2e')];
  om.writeMessage(message,'%s\n')
  message = ['Time step scaled by           ' ... 
      num2str(dtScale,'%4.2e')];
  om.writeMessage(message,'%s\n')
  message = ['New time step size is         ' ...
      num2str(o.dt,'%4.2e')];
  om.writeMessage(message,'%s\n')
  om.writeMessage(' ','%s\n')
else
  accept = true;
  % accept the solution because the error is small
  message = ['Time Step ACCEPTED with error ' ...
      num2str(err,'%4.2e')];
  om.writeMessage(message,'%s\n')
  message = ['Time step scaled by           ' ... 
      num2str(dtScale,'%4.2e')];
  om.writeMessage(message,'%s\n')
  message = ['New time step size is         ' ...
      num2str(o.dt,'%4.2e')];
  om.writeMessage(message,'%s\n')
end

end % newTimeStepSize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,sigma,u,eta,RS,Xrig,etaRig,iter,iflag] = timeStep(o,...
    Xstore,sigStore,uStore,etaStore,RSstore,XrigStore,etaRigStore,...
    deltaX,deltaSig,deltaEta,deltaRS,...
    diffResidual,vesVel,wallVel,kappa,viscCont,walls,wallsCoarse,...
    updatePreco,sa,IK,Xold)
% [X,sigma,u,eta,RS,iter,iflag] = timeStep(o,...
%    Xstore,sigStore,uStore,etaStore,RSstore,...
%    deltaX,deltaSig,deltaEta,deltaRS,...
%    diffResidual,vesVel,wallVel,kappa,viscCont,walls,wallsCoarse,...
%    updatePreco,sa,IK)
% uses explicit or implicit vesicle-vesicle interactions and
% discretizes the inextensibility condition in three different way
% (method1, method2, or method 3).  Must pass in the vesicle positions,
% tension, and velocity from enough previous time steps (depends on
% o.order).  Returns a new positions, tension, velocity, density
% function defined on the solid walls and the number of required GMRES
% iterations if o.SDCcorrect=true, then it uses deltaX, deltaSig, etc
% to compute the right-hand sides needed for sdc updates updatePreco is
% a flog that decides if the block-diagonal preconditioner should be
% updated or not
% NOTE THAT A LOT OF THE FEATURES ARE NOT IMPLEMENTED WITH EXPLICIT
% VESICLE-VESICLE INTERACTIONS.

N = size(Xstore,1)/2; % Number of points per vesicle
nv = size(Xstore,2); % Number of vesicles
if o.confined
  Xwalls = walls.X; % discretization points of solid walls
else
  Xwalls = [];
end
Nbd = size(Xwalls,1)/2; % Number of points on the solid walls
nvbd = size(Xwalls,2); % number of solid wall components
alpha = (1 + viscCont)/2; 
Nrd = size(XrigStore,1)/2;
nvrd = size(XrigStore,2);
% constant that appears in front of time derivative in
% vesicle dynamical equations

Xm = zeros(2*N,nv);
sigmaM = zeros(N,nv);
uM = zeros(2*N,nv);
Xo = zeros(2*N,nv);
etaM = zeros(2*Nbd,nvbd);
RSm = zeros(3,nvbd);
XrigM = zeros(2*Nrd,nvrd);
Xrig = zeros(2*Nrd,nvrd);
etaRigM = zeros(2*Nrd,nvrd);
etaRig = zeros(2*Nrd,nvrd);
for k = 1:o.order
  Xm = Xm + Xstore(:,:,k)*o.Xcoeff(k);
  sigmaM = sigmaM + sigStore(:,:,k)*o.Xcoeff(k);
  uM = uM + uStore(:,:,k)*o.Xcoeff(k);
  if o.confined
    etaM = etaM + etaStore(:,:,k)*o.Xcoeff(k);
    RSm = RSm + RSstore(:,:,k)*o.Xcoeff(k);
  end
  if o.withRigid
    XrigM = XrigM + XrigStore(:,:,k)*o.Xcoeff(k);
    etaRigM = etaRigM + etaRigStore(:,:,k)*o.Xcoeff(k);
  end
  Xo = Xo + Xstore(:,:,k)*o.rhsCoeff(k);
end
% Form linear combinations of previous time steps needed
% for Ascher, Ruuth, and Wetton IMEX methods

vesicle = capsules(Xm,sigmaM,uM,kappa,viscCont,o.antiAlias);
% build an object vesicle that contains tangent vector, jacobian, etc.
rigid = capsules(XrigM,[],[],zeros(nvrd,1),zeros(nvrd,1),o.antiAlias);

op = o.op;
if ~o.SDCcorrect
  o.Galpert = op.stokesSLmatrix(vesicle);
%  o.Gbarnett = op.laplaceSLcomplexMatrix(vesicle);
end
% Build single layer potential matrix and put it in current object
% If we are doing an sdc update, this is already precomputed and 
% stored from when we formed the provisional solution

if ~o.SDCcorrect
  if any(viscCont ~= 1)
    o.D = op.stokesDLmatrix(vesicle);
  else
    o.D = [];
  end
end
% Compute double-layer potential matrix due to each vesicle
% independent of the others.  Matrix is zero if there is no
% viscosity contrast

if ~o.SDCcorrect
  if o.withRigid
    %op = poten(Nrd);
    o.rigidDLP = op.stokesDLmatrix(rigid);
  end
end

if ~o.SDCcorrect
  if o.near
    if o.confined
      [o.NearV2V,o.NearV2W] = vesicle.getZone(walls,3);
      % Need vesicle to vesicle and vesicle to wall interactions
      if nvbd == 1 
        [o.NearW2W,o.NearW2V] = walls.getZone(vesicle,2);
      else
        [o.NearW2W,o.NearW2V] = walls.getZone(vesicle,3);
      end
      % Only need wall to vesicle interactions.  Wall to wall
      % interactions should also use near-singular integration since
      % they may be close to one another
      if o.withRigid
        [o.NearR2R,o.NearR2W] = rigid.getZone(walls,3);
	      [~,o.NearW2R] = walls.getZone(rigid,2);
	      if nv
	        [~,o.NearV2R] = vesicle.getZone(rigid,2);
	        [~,o.NearR2V] = rigid.getZone(vesicle,2);
	      end
      end
    else
      if nv
        o.NearV2V = vesicle.getZone([],1);
      end
      % no solid walls, so only need vesicle-vesicle intearactions
      o.NearV2W = [];
      o.NearW2V = [];
      o.NearW2W = [];
      if o.withRigid
        o.NearR2R = rigid.getZone([],1);
	      if nv
	        [~,o.NearV2R] = vesicle.getZone(rigid,2);
	        [~,o.NearR2V] = rigid.getZone(vesicle,2);
	      end
      end
    end
  else
    o.NearV2V = [];
    o.NearV2W = [];
    o.NearW2V = [];
    o.NearW2W = [];
    o.NearR2R = [];
    o.NearR2V = [];
    o.NearR2W = [];
    o.NearV2R = [];
    o.NearW2R = [];
  end
  % If using near-singular integration, need structures for deciding who
  % is close, how close it is, who is closest, etc.
end
% Only form near-singular integration structure if not doing an SDC
% update.  Otherwise this was formed and saved when forming the
% provisional solution

if o.collDete
  o.lapDLP = op.laplaceDLmatrix(vesicle);
  % build matrix that evaluates double-layer laplace potential for each
  % vesicle

  oc = curve;
  [icollisionVes,icollisionWall] = ...
    oc.collision(vesicle,walls,o.NearV2V,o.NearV2W,o.fmm,o.near);
  % Check for collisions 
  if icollisionVes
    fprintf('VESICLES HAVE CROSSED\n')
    pause
  end
  if icollisionWall
    fprintf('VESICLES HAVE CROSSED SOLID WALL\n')
    pause
  end
end
% check for collisions

if ~o.SDCcorrect
  rhs1 = Xo;
  rhs2 = zeros(N,nv);
  if o.confined
    rhs3 = walls.u;
  else
    rhs3 = [];
  end
  if o.withRigid
    rhs4 = zeros(2*Nrd,nvrd);
    rhs5 = zeros(2,nvrd);
    rhs6 = zeros(1,nvrd);
  else
    rhs4 = [];
    rhs5 = [];
    rhs6 = [];
  end
else
  if strcmp(o.vesves,'explicit')
    rhs1 = diffResidual;
  else
    rhs1 = deltaX + diffResidual;
  end
  if any(vesicle.viscCont ~= 1)
    z = zeros(2*N*nv,1);
    for k = 1:nv
      z(2*(k-1)*N+1:2*k*N) = rhs1(:,k);
    end
    if strcmp(o.vesves,'explicit')
      z = o.IminusDexp(z,vesicle);
    else
      z = o.IminusD(z,vesicle);
    end
    for k = 1:nv
      rhs1(:,k) = z(2*(k-1)*N+1:2*k*N)/alpha(k);
    end
  end
  if strcmp(o.vesves,'explicit')
    rhs1 = rhs1 + deltaX;
  end
  if strcmp(o.solver,'method1')
    rhs2 = ones(N,nv);
  else
    rhs2 = -vesicle.surfaceDiv(vesVel);
  end
  if o.confined
    rhs3 = -wallVel + walls.u;
  else
    rhs3 = [];
  end
  rhs4 = [];
  rhs5 = [];
  rhs6 = [];
end
% Parts of rhs from previous solution.  The right-hand-side depends on
% whether we are doing an SDC correction or forming the provisional
% solution.

% START TO COMPUTE RIGHT-HAND SIDE DUE TO VESICLE TRACTION JUMP
if o.gravity
  % For now, treat gravity explicit same as backgroud flow
  fgrav = o.gravityForce(Xm,vesicle,o.gCont);
else
  fgrav = zeros(2*N,nv);
end
if (strcmp(o.vesves,'explicit') || o.gravity)
  if o.profile
    tic
  end
  if ~o.SDCcorrect
    if strcmp(o.vesves,'implicit')
      f = fgrav;
    else
      f = vesicle.tracJump(Xm,sigmaM) + fgrav;
    end
  else
    if strcmp(o.vesves,'implicit')
      f = vesicle.tracJump(deltaX,deltaSig)*0;
    else
      f = vesicle.tracJump(deltaX,deltaSig);
    end
  end
  if ~o.fmm
    kernel = @op.exactStokesSL;
  else
    kernel = @op.exactStokesSLfmm;
  end

  if ~o.near
    Fslp = kernel(vesicle,f,[]);
    % Evaulate single-layer potential on all vesicles but itself without
    % near-singular integration
    if o.confined
      [~,FSLPwall] = kernel(vesicle,f,[],walls.X,(1:nv));
      % Evaluate single-layer potential on solid walls due to all
      % vesicles
    else
      FSLPwall = [];
    end

  else
    if o.nearStrat == 'cauchy'
      Fslp = op.nearSingStokesSLP(vesicle,f,vesicle,...
          o.Gbarnett,true,o.fmm);
    else
      SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);
      SLPtrap = SLP;
      kernelDirect = kernel;
      Fslp = op.nearSingInt(vesicle,f,SLP,SLPtrap,...
          o.NearV2V,kernel,kernelDirect,vesicle,true);
    end
    % Use near-singular integration to compute single-layer potential
    % due to all other vesicles.  Need to pass function op.exactStokesSL
    % so that the layer potential can be computed at far points and
    % Lagrange interpolation points

    if o.confined
      if o.nearStrat == 'cauchy'
        FSLPwall = op.nearSingStokesSLP(vesicle,f,walls,...
            o.Gbarnett,false,o.fmm);
      else
        FSLPwall = op.nearSingInt(vesicle,f,SLP,SLPtrap,...
            o.NearV2W,kernel,kernelDirect,walls,false);
        % Evaluate the velocity on the walls due to the vesicles
      end
    else
      FSLPwall = [];
    end
    if o.withRigid
    % add sedimentation for rigid particle with implicit scheme
    else
      FSLPrigid = [];
    end
  end
  if o.profile
    fprintf('Build right-hand side V2V and V2W   %5.1e\n',toc);
  end
else
  Fslp = zeros(2*N,nv);
  FSLPwall = zeros(2*Nbd,nvbd);
  FSLPrigid = zeros(2*Nrd,nvrd);
  % vesicle-vesicle and vesicle-wall interactions are handled
  % implicitly in TimeMatVec
end
for k = 1:nv
  rhs1(:,k) = rhs1(:,k) + o.dt/alpha(k)*Fslp(:,k);
end

rhs3 = rhs3 - FSLPwall;
if (o.gravity && ~o.SDCcorrect)
Gfg = op.exactStokesSLdiag(vesicle,o.Galpert,fgrav);
rhs1 = rhs1 + o.dt*Gfg*diag(1./alpha);
end

rhs4 = rhs4 + FSLPrigid;
% END TO COMPUTE RIGHT-HAND SIDE DUE TO VESICLE TRACTION JUMP

% START TO COMPUTE RIGHT-HAND SIDE DUE TO VISCOSITY CONTRAST
if any(vesicle.viscCont ~= 1)
  if o.near
    jump = 1/2*(1-vesicle.viscCont);
    DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle,o.D,X);
  end
  % Need to add jump to double-layer potential if using near-singular
  % integration so that we can compute boundary values for the
  % near-singular integration algorithm

  if ~o.fmmDLP
    kernel = @op.exactStokesDL;
  else
    kernel = @op.exactStokesDLfmm;
  end
  if strcmp(o.vesves,'implicit')
    if ~o.SDCcorrect
      density = Xo;
    else
      density = deltaX;
    end

    if ~o.near
      Fdlp = kernel(vesicle,density,[]);
      if o.confined
        [~,FDLPwall] = kernel(vesicle,density,[],walls.X,(1:nv));
      else
        FDLPwall = [];
      end
    else
      kernelDirect = kernel;
      Fdlp = op.nearSingInt(vesicle,density,DLP,DLP,...
          o.NearV2V,kernel,kernelDirect,vesicle,true);
      % Use near-singular integration to compute double-layer
      % potential from previous solution
      if o.confined
        FDLPwall = op.nearSingInt(vesicle,density,DLP,DLP,...
          o.NearV2W,kernel,kernelDirect,walls,false);
      else
        FDLPwall = [];
      end
      if o.withRigid
        FDLPrigid = op.nearSingInt(vesicle,density,DLP,DLP,...
          o.NearV2R,kernel,kernelDirect,rigid,false);
      else
        FDLPrigid = [];
      end
    end
  elseif strcmp(o.vesves,'explicit')
    if ~o.near
      Fdlp = -o.dt * kernel(vesicle,uM,[]);
      % Evaulate the velocity due to the viscosity contrast
      % on all vesicles but itself WITHOUT near-singular
      % integration
      if o.confined
        [~,FDLPwall] = kernel(vesicle,uM,[],walls.X,(1:nv));
        FDLPwall = -o.dt*FDLPwall;
        % Evaulate the velocity due to the viscosity contrast on the
        % walls WITHOUT near-singulation integration
      else
        FDLPwall = [];
      end
    else
      kernelDirect = kernel;
      Fdlp = -o.dt * op.nearSingInt(vesicle,uM,DLP,DLP,...
          o.NearV2V,kernel,kernelDirect,vesicle,true);
      % Evaulate the velocity due to the viscosity contrast on all
      % vesicles but itself WITH near-singular integration

      if o.confined
        FDLPwall = -o.dt * op.nearSingInt(vesicle,uM,DLP,DLP,...
            o.NearV2W,kernel,kernelDirect,walls,false);
        % Evaulate the velocity due to the viscosity contrast on the
        % walls WITH near-singulation integration
      else
        FDLPwall = [];
      end
    end
    FDLPrigid = zeros(2*Nrd,nvrd);
  end
else
  Fdlp = zeros(2*N,nv);
  FDLPwall = zeros(2*Nbd,nvbd);
  FDLPrigid = zeros(2*Nrd,nvrd);
  % If no viscosity contrast, there is no velocity induced due to a
  % viscosity contrast
end

if (any(viscCont ~= 1) && o.SDCcorrect && strcmp(o.vesves,'explicit'))
  DXo = op.exactStokesDLdiag(vesicle,o.D,deltaX);
  rhs1 = rhs1 - (Fdlp + DXo) * diag(1./alpha);
end
if (any(viscCont ~= 1) && ~o.SDCcorrect)
  DXo = op.exactStokesDLdiag(vesicle,o.D,Xo);
  rhs1 = rhs1 - (Fdlp + DXo) * diag(1./alpha);
end
% add in viscosity contrast term due to each vesicle independent of the
% others (o.D * Xo) from the previous solution followed by the term due
% to all other vesicles (Fdlp)

rhs3 = rhs3 + FDLPwall/o.dt;
rhs4 = rhs4 - FDLPrigid/o.dt;
% compute the double-layer potential due to all other vesicles from the
% appropriate linear combination of previous time steps.  Depends on
% time stepping order and vesicle-vesicle discretization

if o.adhesion
  adhesion = vesicle.adhesionTerm(o.adStrength,o.adRange);
  if ~o.fmm
    kernel = @op.exactStokesSL;
  else
    kernel = @op.exactStokesSLfmm;
  end

  Fadhesion = op.exactStokesSLdiag(vesicle,o.Galpert,adhesion);
  % diagonal term of adhesion

  if ~o.near
    Fadhesion = Fadhesion + kernel(vesicle,adhesion);
    % Evaulate single-layer potential on all vesicles but itself without
    % near-singular integration
  else
    SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);
    SLPtrap = SLP;
    kernelDirect = kernel;
    if o.nearStrat == 'cauchy'
      Fadhesion = Fadhesion + ...
          op.nearSingStokesSLP(vesicle,adhesion,vesicle,...
          o.Gbarnett,true,o.fmm);
    else
      Fadhesion = Fadhesion + ...
          op.nearSingInt(vesicle,adhesion,SLP,SLPtrap,...
          o.NearV2V,kernel,kernelDirect,vesicle,true);
    end
    % Use near-singular integration to compute single-layer potential
    % due to all other vesicles.  Need to pass function op.exactStokesSL
    % so that the layer potential can be computed at far points and
    % Lagrange interpolation points
  end

  rhs1 = rhs1 + o.dt*Fadhesion*diag(1./alpha);
end

% START TO COMPUTE RIGHT-HAND SIDE DUE TO SOLID WALLS
if o.confined
  if strcmp(o.vesves,'explicit')
    if o.profile
      tic
    end
    if o.near
      jump = -1/2;
      DLP = @(X) jump*X + op.exactStokesDLdiag(walls,o.wallDLP,X);
    end
    % compute the matrix for doing evaluating the double-layer potential
    % on the solid walls required for near-singular integration

    if ~o.fmmDLP
      kernel = @op.exactStokesDL;
    else
      kernel = @op.exactStokesDLfmm;
    end
    if ~o.near
      if ~o.SDCcorrect
        charge = etaM;
      else
        charge = deltaEta;
      end
      [~,Fwall2Ves] = kernel(walls,charge,...
          vesicle.X,1:nvbd);
      % velocity field due to the walls evaluated on the vesicle
    else
      if ~o.SDCcorrect
        charge = etaM;
      else
        charge = deltaEta;
      end
      kernelDirect = kernel;
      Fwall2Ves = op.nearSingInt(walls,charge,DLP,DLP,...
          o.NearW2V,kernel,kernelDirect,vesicle,false);
    end

    for k = 2:nvbd
      if ~o.SDCcorrect
        stokeslet = RSm(1:2,k);
        rotlet = RSm(3,k);
      else
        stokeslet = deltaRS(1:2,k);
        rotlet = deltaRS(3,k);
      end
      Fwall2Ves = Fwall2Ves + ...
        o.RSlets(vesicle.X,walls.center(:,k),stokeslet,rotlet);
    end
    if o.profile
      fprintf('Build right-hand side W2V           %5.1e\n',toc);
    end
    Fwall2Rigid = zeros(2*Nrd,nvrd);
  elseif strcmp(o.vesves,'implicit')
    Fwall2Ves = zeros(2*N,nv);
    Fwall2Rigid = zeros(2*Nrd,nvrd);
  end
else
  Fwall2Ves = zeros(2*N,nv);
  Fwall2Rigid = zeros(2*Nrd,nvrd);
  if ~o.SDCcorrect
    rhs1 = rhs1 + o.dt*o.farField(Xm)*diag(1./alpha);
    if o.withRigid
      rhs4 = rhs4 + o.farField(rigid.X);
    end
    % Add in far-field condition (extensional, shear, etc.)
  else
    rhs1 = rhs1 + o.dt*(o.farField(Xold+deltaX)-o.farField(Xold))*diag(1./alpha)*0;
  end
end
rhs1 = rhs1 + o.dt*Fwall2Ves*diag(1./alpha);
rhs4 = rhs4 + Fwall2Rigid;
% right-hand side of the velocity evaluated on the solid walls
% END TO COMPUTE RIGHT-HAND SIDE DUE TO SOLID WALLS


% START TO COMPUTE THE RIGHT-HAND SIDE FOR THE INEXTENSIBILITY CONDITION
if (strcmp(o.solver,'method1'))
  % If using method1, vesicle-vesicle interactions and the presence of a
  % viscosity contrast is irrelevent
  if ~o.SDCcorrect
    rhs2 = rhs2 + vesicle.surfaceDiv(Xo);
  else
    divf = zeros(N,nv);
    for k = 1:nv
      divf(:,k) = curve.arcDeriv(Xm(1:N,k),1,1./sa(:,k),...
          IK(:,k)).^2 + ...
                  curve.arcDeriv(Xm(N+1:2*N,k),1,1./sa(:,k),...
          IK(:,k)).^2; 
    end
    rhs2 = 1/2*(rhs2 - divf);
  end
else 
  % If using method2, method3, or method4, vesicle-vesicle interaction
  % affects the right-hand side of the inextensibility condition

  if ~o.confined && ~o.SDCcorrect
    if any(viscCont ~= 1)
      if o.profile
        tic
      end
      rhs2 = rhs2 - vesicle.surfaceDiv(...
          o.solveIminusD(o.farField(Xm),vesicle));
      if o.profile
        fprintf('Solve system alpha*I - DLP          %5.1e\n',toc);
      end
    else
      rhs2 = rhs2 - vesicle.surfaceDiv(o.farField(Xm));
    end
  end
  % add in term from farfield

  if (strcmp(o.vesves,'explicit'))
    rhs2 = rhs2 - vesicle.surfaceDiv(Fslp); 
    rhs2 = rhs2 - vesicle.surfaceDiv(Fwall2Ves);
    if o.adhesion
      rhs2 = rhs2 + vesicle.surfaceDiv(Fadhesion);
    end
  end
end
% rhs2 is the right-hand side for the inextensibility condition
% END TO COMPUTE THE RIGHT-HAND SIDE FOR THE INEXTENSIBILITY CONDITION

if (any(vesicle.viscCont ~= 1) && ...
      o.confined)
      %strcmp(o.vesves,'implicit') && o.confined)
  rhs3 = rhs3 * o.dt;
end
% This makes sure that the rhs is all order one rather than have rhs3
% being order 1/o.dt and other two parts (rhs1 and rhs2) being order 1.
% This of course needs to be compensated in the TimeMatVec routine

% START TO COMPUTE RIGHT-HAND SIDE DUE TO NEW BENDING TERM
if o.bending
  f = vesicle.newBending;
  % New part of traction jump if we have variable bending
  rhs1 = rhs + o.dt * op.exactStokesSLdiag(vesicle,o.Galpert,f);
end
% This term is always treated fully explicitly.  The highest derivative
% it contains is second-order and it appears through a curvature
% END TO COMPUTE RIGHT-HAND SIDE DUE TO NEW BENDING TERM.  THIS WILL
% ONLY SUPPORT A SINGLE VESICLE IN AN UNBOUNDED FLOW


rhs = [rhs1; rhs2];
rhs = rhs(:);
rhs = [rhs; rhs3(:)];
% Stack the right-hand sides in an alternating with respect to the
% vesicle fashion
% Add on the no-slip boundary conditions on the solid walls
rhs = [rhs; zeros(3*(nvbd-1),1)];
rhs4 = [rhs4;rhs5;rhs6];
rhs = [rhs(:);rhs4(:)];
% Rotlet and Stokeslet equations
usePreco = true;
%usePreco = false;

% START BUILDING BLOCK-DIAGONAL PRECONDITIONER
if usePreco
  if updatePreco
    % only build the preconditioner if updatePreco == true
    if o.profile
      tic
    end
    if nv
      %[Ben,Ten,Div] = vesicle.computeDerivs;
      [Ben1,Ten1,Div1] = o.computeUpDerivs(vesicle);
      [Ben,Ten,Div] = o.upSampleDerivs(Ben1,Ten1,Div1);
    end
    if o.profile
      fprintf('Build differential operators        %5.1e\n',toc);
    end
    % Compute bending, tension, and surface divergence of current
    % vesicle configuration

    if o.profile
      tic
    end

    if (strcmp(o.solver,'method1') || strcmp(o.solver,'method2'))
      bdiagVes.L = zeros(3*N,3*N,nv);
      bdiagVes.U = zeros(3*N,3*N,nv);
    elseif strcmp(o.solver,'method3')
      bdiagVes.GT = zeros(2*N,N,nv); % SLP of tension
      bdiagVes.DGT = zeros(N,N,nv); % divergence of SLP of tension
      bdiagVes.DGB = zeros(N,2*N,nv); % divergence of SLP of bending
      bdiagVes.schur = zeros(2*N,2*N,nv); 
      % schur complement of the lower right N by N block
    elseif strcmp(o.solver,'method4')
      bdiagVes.GT = zeros(2*N,N,nv); % SLP of tension
      bdiagVes.IpBen = zeros(2*N,2*N,nv); 
      % inverse of identity plus SLP of Bending
      bdiagVes.DGB = zeros(N,2*N,nv); % divergence of SLP of bending
      bdiagVes.schur = zeros(N,N,nv); 
      % schur complement of the upper left 2*N by 2*N block
    end

    for k=1:nv
      if strcmp(o.solver,'method1')
        if any(vesicle.viscCont ~= 1)
          [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu(...
            [o.beta*(eye(2*N) - o.D(:,:,k)/alpha(k)) + ...
                o.dt/alpha(k)*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
            -o.dt/alpha(k)*o.Galpert(:,:,k)*Ten(:,:,k); ...
            o.beta*Div(:,:,k) zeros(N)]);
        else
          [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu(...
            [o.beta*eye(2*N) + ...
                o.dt*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
            -o.dt/alpha(k)*o.Galpert(:,:,k)*Ten(:,:,k); ...
            o.beta*Div(:,:,k) zeros(N)]);
        end
      elseif strcmp(o.solver,'method2')
        if any(vesicle.viscCont ~= 1)
          [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu( ...
            [o.beta*(eye(2*N) - o.D(:,:,k)/alpha(k)) + ...
                o.dt/alpha(k)*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
            -o.dt/alpha(k)*o.Galpert(:,:,k)*Ten(:,:,k); ...
            -vesicle.kappa*Div(:,:,k)*(alpha(k)*eye(2*N) - o.D(:,:,k))\...
                (o.Galpert(:,:,k)*Ben(:,:,k)) ...
            Div(:,:,k)*(alpha(k)*eye(2*N) - o.D(:,:,k))\...
                (o.Galpert(:,:,k)*Ten(:,:,k))]);
        else
          [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu(...
            [o.beta*eye(2*N) + ...
                o.dt*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
            -o.dt*o.Galpert(:,:,k)*Ten(:,:,k); ...
            -vesicle.kappa*Div(:,:,k)*o.Galpert(:,:,k)*Ben(:,:,k) ...
            Div(:,:,k)*o.Galpert(:,:,k)*Ten(:,:,k)]);
        end
      elseif strcmp(o.solver,'method3')
        % schur complement of lower right block of method2
        bdiagVes.GT(:,:,k) = o.Galpert(:,:,k)*Ten(:,:,k);
        bdiagVes.DGT(:,:,k) = (Div(:,:,k)*o.Galpert(:,:,k)*...
            Ten(:,:,k))\eye(N);
        bdiagVes.DGB(:,:,k) = ...
          vesicle.kappa*Div(:,:,k)*o.Galpert(:,:,k)*Ben(:,:,k);
        bdiagVes.schur(:,:,k) = ...
          inv((o.beta*eye(2*N) + vesicle.kappa*o.dt*...
          o.Galpert(:,:,k)*Ben(:,:,k)) - ...
          o.dt*o.Galpert(:,:,k)*Ten(:,:,k)*bdiagVes.DGT(:,:,k)*...
          vesicle.kappa*Div(:,:,k)*o.Galpert(:,:,k)*Ben(:,:,k));

      elseif strcmp(o.solver,'method4')
        % schur complement of upper left block of method2
        bdiagVes.GT(:,:,k) = o.Galpert(:,:,k)*Ten(:,:,k);
        bdiagVes.IpBen(:,:,k) = inv(o.beta*eye(2*N) + ...
            vesicle.kappa*o.dt*o.Galpert(:,:,k)*Ben(:,:,k));
        bdiagVes.DGB(:,:,k) = ...
          vesicle.kappa*Div(:,:,k)*o.Galpert(:,:,k)*Ben(:,:,k);
        bdiagVes.schur(:,:,k) = ...
          inv(Div(:,:,k)*(-vesicle.kappa*o.dt*o.Galpert(:,:,k)*...
          Ben(:,:,k)*bdiagVes.IpBen(:,:,k) + eye(2*N))*...
          bdiagVes.GT(:,:,k));
      end
    end
    o.bdiagVes = bdiagVes;
    % Build block-diagonal preconditioner of self-vesicle 
    % intearctions in matrix form
    if o.withRigid
      bdiagRigid.L = zeros(2*Nrd+3,2*Nrd+3,nvrd);
      bdiagRigid.U = zeros(2*Nrd+3,2*Nrd+3,nvrd);
      for k=1:nvrd
        cmR = rigid.center(:,k);
        Xr = rigid.X(:,k);
        saR = rigid.sa(:,k);
        [bdiagRigid.L(:,:,k),bdiagRigid.U(:,:,k)] = ... 
	lu([0.5*eye(2*Nrd)-o.rigidDLP(:,:,k) ... 
	[ones(Nrd,1);zeros(Nrd,1)] [zeros(Nrd,1);ones(Nrd,1)] ... 
	[Xr(Nrd+1:end)-cmR(2);-(Xr(1:Nrd)-cmR(1))];...
        [saR'*2*pi/Nrd zeros(1,Nrd) 0 0 0];[zeros(1,Nrd) saR'*2*pi/Nrd 0 0 0];...
	[2*pi/Nrd*(Xr(Nrd+1:end)'-cmR(2)).*saR' -2*pi/Nrd*(Xr(1:Nrd)'-cmR(1)).*saR' 0 0 0]]);
      end
      o.bdiagRigid = bdiagRigid;
    end
    if o.profile
      fprintf('Build block-diagonal preconditioner %5.1e\n\n',toc);
    end
  end % updatePreco
end % usePreco

if o.profile
  tGMRES = tic;
end


if 1
  warning off
  % any warning is printed to the terminal and the log file so
  % don't need the native matlab version
  if usePreco
    % when doing corrections, expect that these guys will be small
    [Xn,iflag,R,I,resvec] = gmres(@(X) o.TimeMatVec(X,vesicle,walls,0,rigid),...
        rhs,[],o.gmresTol,o.gmresMaxIter...
        ,@o.preconditionerBD);
  else
    [Xn,iflag,R,I,resvec] = gmres(@(X) o.TimeMatVec(X,vesicle,walls,0),...
        rhs,[],o.gmresTol,o.gmresMaxIter);
  end
  warning on
end
% Use GMRES to solve for new positions, tension, density
% function defined on the solid walls, and rotlets/stokeslets
if o.profile
  fprintf('Time to do one time step is         %5.1e\n\n',toc(tGMRES))
end

iter = I(2);
% Number of GMRES iterations

X = zeros(2*N,nv);
sigma = zeros(N,nv);
eta = zeros(2*Nbd,nvbd);
RS = zeros(3,nvbd);
% allocate space for positions, tension, and density function

for k=1:nv
  X(:,k) = Xn((3*k-3)*N+1:(3*k-1)*N);
  sigma(:,k) = Xn((3*k-1)*N+1:3*k*N);
end
% unstack the positions and tensions

% Unstack the positions and tensions
Xn = Xn(3*nv*N+1:end);
for k = 1:nvbd
  eta(:,k) = Xn((k-1)*2*Nbd+1:2*k*Nbd);
end
% unstack the density function
otlets = Xn(2*nvbd*Nbd+1:end);
for k = 2:nvbd
  istart = (k-2)*3+1;
  iend = 3*(k-1);
  RS(:,k) = otlets(istart:iend);
end
% unstack the rotlets and stokeslets
if o.withRigid
  XnRig = Xn(2*nvbd*Nbd+3*(max([(nvbd-1);0]))+1:end);
  [Xrig,etaRig,uRig,omgRig] = o.extractRHSrigid(XnRig,rigid);
end

u = (o.beta*X - Xo)/o.dt;
% Compute the velocity using the differencing stencil

end % timeStep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = TimeMatVec(o,Xn,vesicle,walls,evalOnce,rigid)
% val = TimeMatVec(Xn,vesicle,vesicle,walls) operates the
% vesicle-vesicle and vesicle-boundary interaction formulas to the
% function Xn which contains both the position, tension, and density
% function
global matvecs

matvecs = matvecs + 1;
% counter for the number of matrix-vector multiplications
% that are required for the entire simulation

op = o.op; % poten class
N = vesicle.N; % Number of points
nv = vesicle.nv; % Number of vesicles
if o.confined
  Nbd = walls.N; % Number of points on walls
  nvbd = walls.nv; % Number of components to walls
else
  Nbd = 0;
  nvbd = 0;
end
if o.withRigid
  Nrd = rigid.N;
  nvrd = rigid.nv;
else
  Nrd = 0;
  nvrd = 0;
end

valPos = zeros(2*N,nv);
% right-hand side that corresponds to position equation
valTen = zeros(N,nv);
% right-hand side that corresponds to inextensibilty equation
valWalls = zeros(2*Nbd,nvbd);
% right-hand side that corresponds to solid wall equation
valLets = zeros(3*(nvbd-1),1);
% right-hand side corresponding to the rotlets and stokeslets
valRigid = zeros(2*Nrd,nvrd);
valRigU = zeros(2,nvrd);
valRigOmg = zeros(1,nvrd);
% right-hand side corresponding to RBM

Xm = zeros(2*N,nv);
sigmaM = zeros(N,nv);
for k=1:nv
  Xm(1:2*N,k) = Xn((3*k-3)*N+1:(3*k-1)*N);
  sigmaM(:,k) = Xn((3*k-1)*N+1:3*k*N);
end
% Unstack the position and tension from the input

eta = Xn(3*nv*N+1:3*nv*N+2*Nbd*nvbd+3*max([(nvbd-1);0]));
%eta = Xn(3*nv*N+1:end);
etaM = zeros(2*Nbd,nvbd);
for k = 1:nvbd
  etaM(:,k) = eta((k-1)*2*Nbd+1:2*k*Nbd);
end
% Unstack the density function from the input
otlets = Xn(3*nv*N+2*nvbd*Nbd+1:3*nv*N+2*nvbd*Nbd+3*(max([(nvbd-1);0])));
% stokeslets and rotlets of each vesicle.  Ordered as
% [stokeslet1(component 1);stokeslet1(component 2);rotlet1;...
%  stokeslet2(component 1);stokeslet2(component 2);rotlet2;...];

etaRig = zeros(2*Nrd,nvrd);
uRig = zeros(2,nvrd);
omgRig = zeros(1,nvrd);
if o.withRigid
  XnRig = Xn(3*nv*N+2*nvbd*Nbd+3*(max([(nvbd-1);0]))+1:end);
  XnRig = reshape(XnRig,[],nvrd);
  for k = 1:nvrd
    etaRig(:,k) = XnRig(1:2*Nrd,k);
    uRig(:,k) = XnRig(2*Nrd+1:2*Nrd+2,k);
    omgRig(:,k) = XnRig(end,k);
  end
end


if(evalOnce&&o.gravity)
  fgrav = o.gravityForce(Xm,vesicle,o.gCont);
else
  fgrav = zeros(2*N,nv);
end
f = vesicle.tracJump(Xm,sigmaM) + fgrav;
% f is the traction jump stored as a 2N x nv matrix

alpha = (1+vesicle.viscCont)/2; 
% constant that multiplies the time derivative in the 
% vesicle position equation

Gf = op.exactStokesSLdiag(vesicle,o.Galpert,f);
% Gf is the single-layer potential applied to the traction jump. 

if any(vesicle.viscCont ~= 1)
  DXm = op.exactStokesDLdiag(vesicle,o.D,Xm);
else
  DXm = zeros(2*N,nv);
end
% DXm is the double-layer potential applied to the position
if nv
% START COMPUTING REQUIRED SINGLE-LAYER POTENTIALS
if strcmp(o.vesves,'implicit')
  if o.profile
    tic
  end
  if ~o.fmm
    kernel = @op.exactStokesSL;
  else
    kernel = @op.exactStokesSLfmm;
  end

  if ~o.near
    Fslp = kernel(vesicle,f,[]);
    % Evaulate single-layer potential due to all other vesicles
    % WITHOUT near singular integration.  FMM is optional
    if o.confined
      [~,FSLPwall] = kernel(vesicle,f,[],walls.X,(1:nv));
      % Evaluate single-layer potential due to all vesicles on
      % the solid walls WITHOUT near-singular integration
    else
      FSLPwall = [];
    end
    if o.withRigid
      [~,FSLPrigid] = kernel(vesicle,f,[],rigid.X,(1:nv));
    else
      FSLPrigid = [];
    end
  else
    if o.nearStrat == 'cauchy'
      Fslp = op.nearSingStokesSLP(vesicle,f,vesicle,...
          o.Gbarnett,true,o.fmm);
    else
      SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);
      SLPtrap = SLP;
      kernelDirect = kernel;
      Fslp = op.nearSingInt(vesicle,f,SLP,SLPtrap,...
          o.NearV2V,kernel,kernelDirect,vesicle,true);
    end
    % Evaulate single-layer potential due to all other vesicles
    % WITH near-singular integration.  FMM is optional

    if o.confined
      if o.nearStrat == 'cauchy'
        FSLPwall = op.nearSingStokesSLP(vesicle,f,walls,...
            o.Gbarnett,false,o.fmm);
      else
        FSLPwall = op.nearSingInt(vesicle,f,SLP,SLPtrap,...
          o.NearV2W,kernel,kernelDirect,walls,false);
      end
      % Evaluate single-layer potential due to all vesicles on
      % the solid walls WITH near-singular integration
    else
      FSLPwall = [];
    end

    if o.withRigid
      if o.nearStrat == 'cauchy'
        FSLPrigid = op.nearSingStokesSLP(vesicle,f,rigid,...
            o.Gbarnett,false,o.fmm);
      else
        FSLPrigid = op.nearSingInt(vesicle,f,SLP,SLPtrap,...
          o.NearV2R,kernel,kernelDirect,rigid,false);
      end
    else
      FSLPrigid = [];
    end

  end
  if o.profile
    fprintf('Apply V2V and V2W interactions      %5.1e\n',toc) 
  end
else
  Fslp = zeros(2*N,nv);
  FSLPwall = zeros(2*Nbd,nvbd); 
  FSLPrigid = [];
  % These terms is handled explicitly if o.vesves is 'explicit'
end
% END COMPUTING REQUIRED SINGLE-LAYER POTENTIALS
% Evaluate single-layer potential due to all vesicles except itself and
% the single-layer potential due to all vesicles evaluated on the solid
% walls.  Sets the layer-potential to zero if using explicit
% interactions

% START COMPUTING REQUIRED DOUBLE-LAYER POTENTIALS FOR VISCOSITY
% CONTRAST
if any(vesicle.viscCont ~= 1)
  if o.near
    jump = 1/2*(1-vesicle.viscCont);
    DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle,o.D,X);
  end

  if strcmp(o.vesves,'implicit')
    if ~o.fmmDLP
      kernel = @op.exactStokesDL;
    else
      kernel = @op.exactStokesDLfmm;
    end
    if ~o.near
      Fdlp = kernel(vesicle,Xm,[]);
      if o.confined
        [~,FDLPwall] = kernel(vesicle,Xm,[],walls.X,(1:nv));
      else
        FDLPwall = [];
      end
      if o.withRigid
      else
        FDLPrigid = [];
      end
    else
      kernelDirect = kernel;
      Fdlp = op.nearSingInt(vesicle,Xm,DLP,DLP,...
          o.NearV2V,kernel,kernelDirect,vesicle,true);
      % Use near-singular integration to compute double-layer
      % potential 
      if o.confined
        FDLPwall = op.nearSingInt(vesicle,Xm,DLP,DLP,...
          o.NearV2W,kernel,kernelDirect,walls,false);
      else
        FDLPwall = [];
      end
      if o.withRigid
        FDLPrigid = op.nearSingInt(vesicle,Xm,DLP,DLP,...
          o.NearV2R,kernel,kernelDirect,rigid,false);
      else
        FDLPrigid = [];
      end
    end
  elseif strcmp(o.vesves,'explicit')
    Fdlp = zeros(2*N,nv);
    FDLPwall = zeros(2*Nbd,nvbd);
    FDLPrigid = zeros(2*Nrd,nvrd);
  end
else
  Fdlp = zeros(2*N,nv);
  FDLPwall = zeros(2*Nbd,nvbd);
  FDLPrigid = zeros(2*Nrd,nvrd);
end
% END COMPUTING REQUIRED DOUBLE-LAYER POTENTIALS FOR VISCOSITY CONTRAST
else
  Fslp = [];
  Fdlp = [];
  FSLPwall = zeros(2*Nbd,nvbd);
  FSLPrigid = zeros(2*Nrd,nvrd);
  FDLPwall = zeros(2*Nbd,nvbd);
  FDLPrigid = zeros(2*Nrd,nvrd);
end

% START OF EVALUATING DOUBLE-LAYER POTENTIALS DUE TO SOLID WALLS
if o.confined
  if strcmp(o.vesves,'implicit')
    if ~o.fmmDLP
      kernel = @op.exactStokesDL;
    else
      kernel = @op.exactStokesDLfmm;
    end
    if ~o.near
      [~,Fwall2Ves] = kernel(walls,etaM,[],vesicle.X,1:nvbd);
    else
      jump = -1/2;
      DLP = @(X) jump*X + op.exactStokesDLdiag(walls,o.wallDLP,X);
      kernelDirect = kernel;
      Fwall2Ves = op.nearSingInt(walls,etaM,DLP,DLP,...
          o.NearW2V,kernel,kernelDirect,vesicle,false);
    end

    if o.withRigid
      if ~o.near
        [~,Fwall2Rigid] = kernel(walls,etaM,[],rigid.X,1:nvbd);
      else
        jump = -1/2;
        DLP = @(X) jump*X + op.exactStokesDLdiag(walls,o.wallDLP,X);
        kernelDirect = kernel;
        Fwall2Rigid = op.nearSingInt(walls,etaM,DLP,DLP,...
            o.NearW2R,kernel,kernelDirect,rigid,false);
      end
    else
      Fwall2Rigid = [];
    end

    if o.profile
      fprintf('Apply W2V interaction               %5.1e\n',toc) 
    end
    % compute the velocity on the vesicles due to the solid walls
  else
    Fwall2Ves = zeros(2*N,nv);
    Fwall2Rigid = zeros(2*Nrd,nvrd);
  end
else
  Fwall2Ves = zeros(2*N,nv);
  Fwall2Rigid = zeros(2*Nrd,nvrd);
end
% END OF EVALUATING DOUBLE-LAYER POTENTIALS DUE TO SOLID WALLS

% START OF EVALUATING WALL TO WALL INTERACTIONS
if (o.confined && nvbd > 1)
  if o.profile
    tic
  end
  % only need to do wall to wall interactions if the domain is multiply
  % connected
  if ~o.fmmDLP
    kernel = @op.exactStokesDL;
    FDLPwall2wall = kernel(walls,etaM,[]);
  else
    kernel = @op.exactStokesDLfmm;
    FDLPwall2wall = zeros(2*walls.N,walls.nv);
    for k = 1:nvbd
      isou = [(1:k-1) (k+1:nvbd)];
      [~,FDLPwall2wall(:,k)] = kernel(walls,etaM,[],walls.X(:,k),isou);
    end
    % Since the double-layer potential is still expensive even with the
    % FMM, this eliminates the need to do one big FMM followed by a
    % bunch of small ones to subtract off the self-interaction term
    % which is calculated using the precomputed matrix
  end

  if o.profile
    fprintf('Apply W2W interaction               %5.1e\n\n',toc)
  end
end
% END OF EVALUATING WALL TO WALL INTERACTIONS

% START OF EVALUATING POTENTIAL DUE TO STOKESLETS AND ROTLETS
if nvbd > 1
  LetsWalls = zeros(2*Nbd,nvbd);
  LetsVes = zeros(2*N,nv);
  LetsRigid = zeros(2*Nrd,nvrd);
  for k = 2:nvbd
    stokeslet = otlets(3*(k-2)+1:3*(k-2)+2);
    rotlet = otlets(3*(k-1));
    if strcmp(o.vesves,'implicit')
      LetsVes = LetsVes + o.RSlets(vesicle.X,walls.center(:,k),...
          stokeslet,rotlet);
      if o.withRigid
        LetsRigid = LetsRigid + o.RSlets(rigid.X,walls.center(:,k),...
            stokeslet,rotlet);
      end
    end
    % compute velocity due to rotlets and stokeslets on the vesicles

    LetsWalls = LetsWalls + o.RSlets(walls.X,walls.center(:,k),...
        stokeslet,rotlet);
    % compute velocity due to rotlets and stokeslets on the solid walls
  end
  valLets = o.letsIntegrals(otlets,etaM,walls);
  % Integral constraints on the density function eta related
  % to the weights of the stokeslets and rotlets
else
  LetsVes = zeros(2*N,nv);
  LetsWalls = zeros(2*Nbd,nvbd);
  FDLPwall2wall = zeros(2*Nbd,nvbd);
  LetsRigid = zeros(2*Nrd,nvrd);
end
% END OF EVALUATING POTENTIAL DUE TO STOKESLETS AND ROTLETS

if o.withRigid
  if strcmp(o.vesves,'implicit')
    if o.profile
      tic
    end

    if ~o.fmmDLP
      kernel = @op.exactStokesDL;
    else
      kernel = @op.exactStokesDLfmm;
    end

    if o.near
      jump = -1/2;
      DLP = @(X) jump*X + op.exactStokesDLdiag(rigid,o.rigidDLP,X);
      %DLP = @(X) jump*X + op.exactStokesRigidDLdiag(rigid,o.rigidDLP,theta,X);
    end

    if ~o.near
      if ~o.SDCcorrect
        charge = etaRig;
      else
        charge = deltaRigEta;
      end
      if nv
        [~,Frigid2Ves] = kernel(rigid,charge,...
            vesicle.X,1:nvrd);
      else
        Frigid2Ves = [];
      end

      if o.confined
        [~,Frigid2Wall] = kernel(rigid,charge,...
            walls.X,1:nvrd);
      else
        Frigid2Wall = [];
      end
      
      FDLPrigid2rigid = kernel(rigid,charge,[]);
    else
      if ~o.SDCcorrect
        charge = etaRig;
      else
        charge = deltaRigEta;
      end
      kernelDirect = kernel;

      if nv
        Frigid2Ves = op.nearSingInt(rigid,charge,DLP,DLP,...
            o.NearR2V,kernel,kernelDirect,vesicle,false);
      else
        Frigid2Ves = [];
      end

      if o.confined
        Frigid2Wall = op.nearSingInt(rigid,charge,DLP,DLP,...
            o.NearR2W,kernel,kernelDirect,walls,false);
      else
        Frigid2Wall = [];
      end

      FDLPrigid2rigid = op.nearSingInt(rigid,charge,DLP,DLP,...
          o.NearR2R,kernel,kernelDirect,rigid,true);
    end

    if o.profile
      fprintf('Build right-hand side R2V           %5.1e\n',toc);
    end
  else
    Frigid2Ves = zeros(2*N,nv);
    Frigid2Wall = zeros(2*Nbd,nvbd);
    FDLPrigid2rigid = zeros(2*Nrd,nvrd);
  end
else
  FDLPrigid2rigid = zeros(2*Nrd,nvrd);
  Frigid2Ves = zeros(2*N,nv);
  Frigid2Wall = zeros(2*Nbd,nvbd);
end


% START OF EVALUATING VELOCITY ON VESICLES
valPos = valPos - o.dt*Gf*diag(1./alpha);
% self-bending and self-tension terms
valPos = valPos - o.beta*DXm*diag(1./alpha);
% self-viscosity contrast term
valPos = valPos - o.dt*Fslp*diag(1./alpha);
% single-layer potential due to all other vesicles
valPos = valPos - o.beta*Fdlp*diag(1./alpha);
% double-layer potential due to all other vesicles
valPos = valPos - o.dt*Fwall2Ves*diag(1./alpha);
% velocity due to solid walls evaluated on vesicles
valPos = valPos - o.dt*LetsVes*diag(1./alpha);
valPos = valPos - o.dt*Frigid2Ves*diag(1./alpha);
% velocity on vesicles due to the rotlets and stokeslets
% END OF EVALUATING VELOCITY ON VESICLES

% START OF EVALUATING VELOCITY ON WALLS
if o.confined
  valWalls = valWalls - 1/2*etaM + ...
    op.exactStokesDLdiag(walls,o.wallDLP,etaM);
  valWalls(:,1) = valWalls(:,1) + ...
    op.exactStokesN0diag(walls,o.wallN0,etaM(:,1));
end
% evaluate velocity on solid walls due to the density function.
% self solid wall interaction

valWalls = valWalls + FSLPwall;
% velocity on walls due to the vesicle traction jump
valWalls = valWalls + o.beta*FDLPwall/o.dt;
% velocity on walls due to the vesicle viscosity jump
valWalls = valWalls + FDLPwall2wall;
% velocity on walls due to all other walls
valWalls = valWalls + LetsWalls;
valWalls = valWalls + Frigid2Wall;
% velocity on walls due to the rotlets and stokeslets
% END OF EVALUATING VELOCITY ON WALLS

% START OF EVALUATING INEXTENSIBILITY CONDITION
if (strcmp(o.solver,'method1'))
  valTen = o.beta * vesicle.surfaceDiv(Xm);
  % compute surface divergence of the current GMRES iterate
  % method1 sets this equal to the surface divergence of
  % the previous time step
else
  if any(vesicle.viscCont ~= 1)
    if o.profile
      tic
    end
    valTen = vesicle.surfaceDiv(...
        o.solveIminusD(Gf+Fslp+Fwall2Ves+LetsVes,vesicle));
    if o.profile
      fprintf('Solve system alpha*I - DLP          %5.1e\n',toc);
    end
  else
    valTen = -1/o.dt*vesicle.surfaceDiv(valPos);
  end
  % method2, method3, and method4 sets the surface divergence of the sum
  % of single-layer potentials due to bending and tension plus the
  % farField to zero.  The only difference between the two methods is
  % how the preconditioner is used.  method2 uses the possibly
  % ill-conditioned full matrix where as method3 and method 4 use the
  % two schur complements.  Eventually will phase out method2 and it
  % will be fully replaced by method3 and method4
end
% Two possible discretizations of the inextensibility condition
% END OF EVALUATING INEXTENSIBILITY CONDITION

valPos = valPos + o.beta*Xm;
% beta times solution coming from time derivative

val = zeros(3*N*nv,1);
% Initialize output from vesicle and inextensibility equations to zero
for k=1:nv
  val((k-1)*3*N+1:3*k*N) = [valPos(:,k);valTen(:,k)];
end
% Stack val as [x-coordinate;ycoordinate;tension] repeated
% nv times for each vesicle

if o.withRigid
  for i = 1:nvrd
    X = rigid.X(:,i);
    sa = rigid.sa(:,i);
    cm = rigid.center(:,i);
    eta = etaRig(:,i);

    stokeslet = zeros(2,1);
    rotlet = 0;
    stokeslet(1) = sum(eta(1:Nrd).*sa)*2*pi/Nrd; % force_X
    stokeslet(2) = sum(eta(Nrd+1:end).*sa)*2*pi/Nrd; % force_y
    rotlet = sum(...
      ((X(Nrd+1:end)-cm(2)).*eta(1:Nrd) - ...
      (X(1:Nrd)-cm(1)).*eta(Nrd+1:end)).*...
      sa)*2*pi/Nrd; % torque
    
    valRigU(:,i) = stokeslet;
    valRigOmg(:,i) = rotlet;
    
    stokeslet = stokeslet/(2*pi);
    rotlet = rotlet/(2*pi);

    letsvel = o.RSlets(X,cm,stokeslet,rotlet);
    
    U = repmat(uRig(:,i)',Nrd,1) + omgRig(i)*[X(Nrd+1:end)-cm(2),-(X(1:Nrd)-cm(1))];
    valRigid(:,i) = 1/2*eta - o.rigidDLP(:,:,i)*eta + U(:);% - letsvel;
  end
end

valRigid = valRigid - FSLPrigid;
% velocity on walls due to the vesicle traction jump
valRigid = valRigid - o.beta*FDLPrigid/o.dt;
% velocity on walls due to the vesicle viscosity jump
valRigid = valRigid - Fwall2Rigid;
% velocity on walls due to all other walls
valRigid = valRigid - LetsRigid;
valRigid = valRigid - FDLPrigid2rigid;
valRigid = [valRigid;valRigU;valRigOmg];

if (any(vesicle.viscCont ~= 1) && ...
      o.confined)
      %strcmp(o.vesves,'implicit') && o.confined)
  % This combination of options causes problems with
  % the scaling of the preconditioner.  Need to
  % get rid of the potentially small value o.dt
  valWalls = valWalls * o.dt;
end

val = [val;valWalls(:);valLets;valRigid(:)];
% Stack velocity along the solid walls in same manner as above
% Stack the stokeslets and rotlet componenets at the end

end % TimeMatVec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vesVel,divVesVel,wallVel,residual] = ...
      computeVelocities(o,vesicle,Galpert,Gbarnett,D,walls,...
          etaProv,RSprov,NearV2V,NearV2W,NearW2V,NearW2W,uProv)
% [vesVel,divVesVel,wallVel,residual] = ...
%      computeVelocities(vesicle,Galpert,Gbarnett,D,walls,...
%       etaProv,RSprov,NearV2V,NearV2W,NearW2V,NearW2W);
% computes the velocity induced by the provisional vesicle position,
% tension, and density function on the vesicles and the walls.  Also
% returns the vesicle divergence of the velocity field and the residual
% of the picard integral formulation of the vesicle velocity.  These
% quantities are needed to form the modificiations of the right-hand
% side when doing SDC updates

N = vesicle.N;
nv = vesicle.nv;
if o.confined
  Nbd = walls.N;
  nvbd = walls.nv;
else
  Nbd = 0;
  nvbd = 0;
end

vesVel = zeros(2*N,nv,abs(o.orderGL));
% velocity on the vesicles due to the provisional solution 
divVesVel = zeros(N,nv,abs(o.orderGL));
% vesicle divergence of the vesicles due to the provisional
% solution
wallVel = zeros(2*Nbd,nvbd,abs(o.orderGL));
% velocity on the solid walls due to the provisional solution
residual = zeros(2*N,nv,abs(o.orderGL));
% residual of the Picard integral coming from the time
% derivative term of the vesicle position

vesves = o.vesves;
o.vesves = 'implicit';
% need to use implicit so that other vesicles are used to compute
% the integrand z
order = o.order;
o.order = 1;
[o.Xcoeff,o.rhsCoeff,o.beta] = o.getCoeff(o.order);
% need to save the time stepping order
dt = o.dt;
% need to save current time step
o.dt = 1;
% avoid introducing numerical error by multiplying and dividing
% by a potentially small number

if ~o.confined
  z = zeros(3*N*nv,1);
else
  z = zeros(3*N*nv + 2*Nbd*nvbd + 3*(nvbd-1),1);
end
op = o.op;
for n = 1:abs(o.orderGL)
  % loop over the configurations
  o.Galpert = Galpert(:,:,:,n);
%  o.Gbarnett = Gbarnett(:,:,:,n);
  if any(vesicle(n).viscCont ~= 1)
    o.D = D(:,:,:,n);
  end
  % Already computed the single-layer potential when forming the
  % provisional solution Store the single- and double-layer potentials
  % in the object o
  o.NearV2V = NearV2V{n};
  o.NearV2W = NearV2W{n};
  o.NearW2V = NearW2V{n};
  o.NearW2W = NearW2W{n};
  % retrieve near-singular integration strucutre at Gauss-Lobatto time
  % step n

  for k = 1:nv
    z(3*(k-1)*N+1:3*k*N) = [vesicle(n).X(:,k);vesicle(n).sig(:,k)];
  end
  % put positions and tension into z the way that TimeMatVec requires

  for k = 1:nvbd
    istart = 3*nv*N+2*(k-1)*Nbd+1;
    iend = istart + 2*Nbd - 1;
    z(istart:iend) = etaProv(:,k,n);
  end
  % put in density function on solid walls into z the way that
  % TimeMatVec requires them

  for k = 2:nvbd
    istart = 3*nv*N+2*nvbd*Nbd+3*(k-2)+1;
    iend = istart + 2;
    z(istart:iend) = RSprov(:,k,n);
  end
  % put the rotlet and stokeslet coefficients into z

  viscCont = vesicle(n).viscCont;
  vesicle(n).viscCont = ones(1,nv);
  % don't want the double-layer contribution when computing the velocity
  % of the vesicle
  z = o.TimeMatVec(z,vesicle(n),walls,1);
  % use TimeMatVec to find velocity on vesicles and solid
  % walls due to the provisional solution
  vesicle(n).viscCont = viscCont;
  % set the viscosity contrast back to its original value

  rhs = zeros(2*N,nv);
  for k = 1:nv
    istart = 3*(k-1)*N+1;
    iend = istart + 2*N - 1;
    rhs(:,k) = -z(istart:iend);
  end
  if ~o.confined
    rhs = rhs + o.farField(vesicle(n).X);
  end

  vesVel(:,:,n) = vesicle(n).X + rhs;
  % form the velocity on the vesicle due to the current provisional
  % solution

  if any(vesicle(n).viscCont ~= 1)
    if o.profile
      tic
    end
    if strcmp(vesves,'explicit')
     
      jump = 1/2*(1-vesicle(n).viscCont);
      DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle(n),o.D,X);
      if ~o.fmmDLP
        kernel = @op.exactStokesDL;
      else
        kernel = @op.exactStokesDLfmm;
      end
      kernelDirect = kernel;
      den = uProv(:,:,n);
      Fdlp = op.nearSingInt(vesicle(n),den,DLP,DLP,...
          o.NearV2V,kernel,kernelDirect,vesicle(n),true);
      vesVel(:,:,n) = vesVel(:,:,n) + Fdlp;
      vesVel(:,:,n) = o.solveIminusDexp(vesVel(:,:,n),vesicle(n));
    
      %vesVel(:,:,n) = o.solveIminusD(vesVel(:,:,n),vesicle(n));
    else
      vesVel(:,:,n) = o.solveIminusD(vesVel(:,:,n),vesicle(n));
    end
    if o.profile
      fprintf('Solve system alpha*I - DLP          %5.1e\n',toc);
    end
  end 
  % need to apply inv(alpha*I - DLP) if there is a viscosity contrast to
  % obtain the velocity of the vesicles.

  divVesVel(:,:,n) = vesicle(n).surfaceDiv(vesVel(:,:,n));
  % surface divergence of the velocity of the vesicle due
  % to the provisional solution

  for k = 1:nvbd
    istart = 3*nv*N + (k-1)*2*Nbd + 1;
    iend = istart + 2*Nbd - 1;
    wallVel(:,k,n) = z(istart:iend);
  end

  if any(vesicle(n).viscCont ~= 1) && o.confined
    jump = 1/2*(1-vesicle(n).viscCont);
    DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle(n),o.D,X);
    if ~o.fmmDLP
      kernel = @op.exactStokesDL;
    else
      kernel = @op.exactStokesDLfmm;
    end

    kernelDirect = kernel;
    den = vesVel(:,:,n);
    FDLPwall = op.nearSingInt(vesicle(n),den,DLP,DLP,...
        o.NearV2W,kernel,kernelDirect,walls,false);
    % Evaulate the velocity due to the viscosity contrast
    % on the walls WITH near-singulation integration
    wallVel(:,:,n) = wallVel(:,:,n) + FDLPwall;
  end
  % velocity on the solid walls due to the provisional solution
end

o.dt = dt;
o.order = order;
[o.Xcoeff,o.rhsCoeff,o.beta] = o.getCoeff(o.order);
o.vesves = vesves;
% change back to original time step and vesicle-vesicle interaction

IvesVel = o.lobattoInt(vesVel);
% integrate the vesicles velocity using quadrature rules 
% that are exact for polynomials defined at the 
% Gauss-Lobatto points
for n = 1:abs(o.orderGL)
  residual(:,:,n) = vesicle(1).X - vesicle(n).X + ...
      o.dt/2 * IvesVel(:,:,n);
end
% compute residual by adding the initial vesicle configuartion and
% subtracting the current vesicle configuartion

end % computeVelocities

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zOut = solveIminusD(o,zIn,vesicle)
% zOut = o.solveIminusD(zIn,vesicle) inverts the linear system \alpha*I
% - DLP where alpha=0.5*(1+viscCont) and DLP is the double-layer
% potential.  This is the opeartor that needs to be inverted to find the
% velocity when there is a viscosity contrast involved.  If there is no
% viscosity contrast, alpha = 1 and DLP = 0 so that zOut = zIn


warning off
%tGMRES = tic;
%[zIn,flag,relres,iter] = gmres(@(X) o.IminusD(X,vesicle),zIn(:),...
%    [],1e-2*o.gmresTol,min(2*vesicle.N*vesicle.nv,o.gmresMaxIter));
[zIn,flag,relres,iter] = gmres(@(X) o.IminusD(X,vesicle),zIn(:),...
    [],1e-2*o.gmresTol,min(2*vesicle.N*vesicle.nv,o.gmresMaxIter),...
    @(z) o.precoIminusD(z,vesicle));
warning on
%fprintf('GMRES time is %4.2e\n',toc(tGMRES))
% solve with block-diagonal preconditioned GMRES.  Integral equation is
% of the form identity + compact.  Need a bit more accuracy in this
% gmres solve as it is an inner iteration within an outer GMRES
% iteration.  This hides the fact that the GMRES solver is not linear

zOut = zeros(2*vesicle.N,vesicle.nv);
for k = 1:vesicle.nv
  zOut(:,k) = zIn((k-1)*2*vesicle.N+1:2*k*vesicle.N);
end
% Sort output in appropriate format

end % solveIminusD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = IminusD(o,X,vesicle)
% val = IminusD(X,vesicle) does a matrix-vector multiply with the matrix
% (alpha*I - D) where alpha = 0.5*(1+viscCont) and D is the double-layer
% potential.  This matrix needs to be inverted when solving for the
% velocity field of a given vesicle configuration and tension

op = o.op;
N = vesicle.N;
nv = vesicle.nv;
alpha = 0.5*(1+vesicle.viscCont);

Xm = zeros(2*N,nv);
for k = 1:nv
  Xm(:,k) = X((k-1)*2*N+1:k*2*N);
end

val = Xm*diag(alpha);
% "jump" term since we are computing alpha * I - DLP

if vesicle.N > 128
  val = val - op.exactStokesDLdiag(vesicle,o.D,Xm);
else
  val = val - op.exactStokesDLdiag(vesicle,[],Xm);
end
% self-interaction term

if o.near
  jump = 0.5*(1-vesicle.viscCont);
  DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle,o.D,X);
end
% compute term for evaluating limiting value along boundary of domain
if ~o.fmmDLP || N*nv <=  2000
  kernel = @op.exactStokesDL;
else
  kernel = @op.exactStokesDLfmm;
end
% kernel for evaulating the double-layer potential

%tic
if ~o.near
  Fdlp = kernel(vesicle,Xm,[]);
else
  kernelDirect = kernel;
  Fdlp = op.nearSingInt(vesicle,Xm,DLP,DLP,...
    o.NearV2V,kernel,kernelDirect,vesicle,true);
end
%fprintf('FMM other interactions %4.2e\n',toc)
% potential due to other vesicles

val = val - Fdlp;
% add potential due to self and other vesicles.  Jump is already
% included
val = val(:);

end % IminusD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = sigDenMatVec(o,sigma,vesicle,walls)
% val = sigDenMatVec(sigma,vesicle,walls,NearV2V,NearV2W) does the
% matvec multiply required to find the tension and density function of a
% given vesicle configuration and farfield or solid wall velocity field

N = vesicle.N; % Number of points
nv = vesicle.nv; % Number of vesicles
if o.confined
  Nbd = walls.N; % Number of points per wall
  nvbd = walls.nv; % Number of walls
else
  Nbd = 0;
  nvbd = 0;
end

u = zeros(2*N,nv);
% need to do some manipulations to the TimeMatVec matrix applied only
% to the first component, which is the velocity
for k = 1:nv
  istart = (k-1)*3*N + 1;
  iend = istart + 2*N - 1;
  u(:,k) = sigma(istart:iend);
end

vesves = o.vesves;
o.vesves = 'implicit';
% want interactions to be implicit so that the most accurate tension and
% density functions are found
inextens = o.solver;
o.solver = 'method1';
dt = o.dt;
o.dt = 1;

op = o.op;
sigma = o.TimeMatVec(sigma,vesicle,walls,0);
% Do a matvec but let the incoming postitions be zero since they are
% handled by the initial condition
% Can overwrite sigma as we don't need it from this point onwards
o.vesves = vesves;
o.solver = inextens;
o.dt = dt;
% change back to old vesicle-vesicle and vesicle-boundary interactions

valVel = zeros(2*N,nv);
valTen = zeros(N,nv);
valDen = zeros(2*Nbd,nvbd);
valRS = zeros(3,(nvbd-1));
for k = 1:nv
  valVel(:,k) = sigma((k-1)*3*N+1:(3*k-1)*N);
  valTen(:,k) = sigma((3*k-1)*N+1:3*k*N);
end
% part that corresponds to the velocity and tension
for k = 1:nvbd
  valDen(:,k) = sigma(3*nv*N+(k-1)*2*Nbd+1:3*nv*N+k*2*Nbd);
end
% part that corresponds to the density function
for k = 2:nvbd
  valRS(:,k-1) = sigma(3*nv*N+2*Nbd*nvbd+1:end);
end
% part that corresponds to the Stokeslets and Rotlets

if ~o.fmm
  kernel = @op.exactStokesSL;
else
  kernel = @op.exactStokesSLfmm;
end
% kernel for single-layer potential.  Only difference is if the FMM is
% used or not

f = vesicle.tracJump(u,zeros(N,nv));
% bending due to the velocity
if ~o.near
  Fslp = kernel(vesicle,f,[]);
  if o.confined
    [~,FSLPwall] = kernel(vesicle,f,[],walls.X,(1:nv));
  else
    FSLPwall = [];
  end
else
  SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);
  SLPtrap = SLP;
  kernelDirect = kernel;
  Fslp = op.nearSingInt(vesicle,f,SLP,SLPtrap,...
      o.NearV2V,kernel,kernelDirect,vesicle,true);

  if o.confined
    FSLPwall = op.nearSingInt(vesicle,f,SLP,SLPtrap,...
      o.NearV2W,kernel,kernelDirect,walls,false);
    % Evaluate single-layer potential due to all vesicles on
    % the solid walls WITH near-singular integration
  else
    FSLPwall = [];
  end
end

alpha = (1+vesicle.viscCont)/2; 
valVel = valVel * diag(alpha);
% multiply top row of matrix by alpha

valVel = valVel + op.exactStokesSLdiag(vesicle,o.Galpert,f) + Fslp;
valDen = valDen - FSLPwall;
% subtract off terms that TimeMatVec introduces but we do not have in
% this linear system

val = [valVel;valTen];
val = [val(:);valDen(:);valRS(:)];
% stack the different components coming from the inextensibility, solid
% walls, and rotlets/stokeslets

end % sigDenMatVec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = letsIntegrals(o,otlets,etaM,walls)
% z = letsIntegrals(stokeslet,rotlet,etaM,walls) integrates the density
% function to enforce constraints on stokeslets and rotlets

Nbd = walls.N;
nvbd = walls.nv;
z = zeros(3*(nvbd-1),1);

for k = 2:nvbd
  stokeslet = otlets(3*(k-2)+1:3*(k-2)+2);
  % two stokeslet terms per inner boundary
  rotlet = otlets(3*(k-1));
  % one rotlet term per inner boundary
  ind = 3*(k-2)+1;
  z(ind) = -2*pi*stokeslet(1) + ...
    sum(etaM(1:Nbd,k).*walls.sa(:,k))*2*pi/Nbd;
  % integral of density function dotted with [1;0]
  % is one stokeslet
  z(ind+1) = -2*pi*stokeslet(2) + ...
    sum(etaM(Nbd+1:2*Nbd,k).*walls.sa(:,k))*2*pi/Nbd;
  % integral of density fuction dotted with [0;1]
  % is the other stokeslet
  z(ind+2) = -2*pi*rotlet + sum(...
    ((walls.X(Nbd+1:2*Nbd,k)).*etaM(1:Nbd,k) - ...
    (walls.X(1:Nbd,k)).*etaM(Nbd+1:2*Nbd,k)).*...
    walls.sa(:,k))*2*pi/Nbd;
  % integral of density function dotted with (-y,x)
  % is the rotlet
end % k

end % letsIntegrals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vel = RSlets(o,X,center,stokeslet,rotlet)
% vel = RSlets(o,X,center,stokeslet,rotlet) evaluates the velocity due
% to the stokeslet and rotlet terms.  Center of the rotlet and
% stokeslet is contained in center

oc = curve;
[x,y] = oc.getXY(X);
% set of points where we are evaluating the velocity
[cx,cy] = oc.getXY(center);
% the center of the rotlet/stokeslet terms

rho2 = (x-cx).^2 + (y-cy).^2;
% distance squared

LogTerm = -0.5*log(rho2)*stokeslet(1);
rorTerm = 1./rho2.*((x-cx).*(x-cx)*stokeslet(1) + ...
    (x-cx).*(y-cy)*stokeslet(2));
RotTerm = (y-cy)./rho2*rotlet;
velx = 1/4/pi*(LogTerm + rorTerm) + RotTerm;
% x component of velocity due to the stokeslet and rotlet

LogTerm = -0.5*log(rho2)*stokeslet(2);
rorTerm = 1./rho2.*((y-cy).*(x-cx)*stokeslet(1) + ...
    (y-cy).*(y-cy)*stokeslet(2));
RotTerm = -(x-cx)./rho2*rotlet;
vely = 1/4/pi*(LogTerm + rorTerm) + RotTerm;
% y component of velocity due to the stokeslet and rotlet

vel = [velx;vely];
% velocity

end % RSlets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF DIFFERENT PRECONDITIONERS INCLUDING BLOCK-DIAGONAL, ONE FOR
% THE SYSTEM THAT SOLVES FOR THE TENSION AND DENSITY GIVEN A POSITION,
% MULTIGRID IDEAS, SCHUR COMPLEMENTS, AND ANALYTIC (BASED ON A
% CIRCLE).  THE BLOCK-DIAGONAL PRECONDITIONER IS THE MOST ROBUST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = preconditionerBD(o,z)
% val = preconditionBD(z) applies the block diagonal preconditioner
% required by preconditioned-GMRES to the vector z

if (strcmp(o.solver,'method1') || strcmp(o.solver,'method2'))
  nv = size(o.bdiagVes.L,3); % number of vesicles
  N = size(o.bdiagVes.L,1)/3; % number of points
  if o.withRigid
    nvrd = size(o.bdiagRigid.L,3);
    Nrd = (size(o.bdiagRigid.L,1)-3)/2; % number of points
  else
    nvrd = 0;
    Nrd = 0;
  end
elseif strcmp(o.solver,'method3')
  nv = size(o.bdiagVes.schur,3); % number of vesicles
  N = size(o.bdiagVes.schur,1)/2; % number of vesicles
elseif strcmp(o.solver,'method4')
  nv = size(o.bdiagVes.schur,3); % number of vesicles
  N = size(o.bdiagVes.schur,1); % number of vesicles
end
nvbd = size(o.wallDLP,3); % number of solid walls
Nbd = size(o.wallDLP,1)/2; % number of points per solid wall

zves = z(1:3*N*nv);
% extract the position and tension part.  Solid walls is
% handled in the next section of this routine
valVes = zeros(3*N*nv,1);
valRigid = zeros(2*Nrd+3,nvrd);

if (strcmp(o.solver,'method1') || strcmp(o.solver,'method2'))
  for k=1:nv
    valVes((k-1)*3*N+1:3*k*N) = o.bdiagVes.U(:,:,k)\...
      (o.bdiagVes.L(:,:,k)\zves((k-1)*3*N+1:3*k*N));
  end % k
  % precondition with the block diagonal preconditioner for the
  % vesicle position and tension
elseif strcmp(o.solver,'method3')
  for k = 1:nv
    istart = (k-1)*3*N+1;
    iend = istart + 2*N - 1;
    bx = zves(istart:iend); 
    istart = iend + 1;
    iend = istart + N - 1;
    bsig = zves(istart:iend);
    % seperate position and tension terms

    rhs = bx + o.dt*o.bdiagVes.GT(:,:,k)*o.bdiagVes.DGT(:,:,k)*bsig;
    pos = o.bdiagVes.schur(:,:,k)*rhs;
    rhs = bsig + o.bdiagVes.DGB(:,:,k)*pos;
    sig = o.bdiagVes.DGT(:,:,k)*rhs;
    valVes((k-1)*3*N+1:3*k*N) = [pos;sig];
    % use schur decomposition to operate preconditioner
  end % k
elseif strcmp(o.solver,'method4')
  for k = 1:nv
    istart = (k-1)*3*N+1;
    iend = istart + 2*N - 1;
    bx = zves(istart:iend); 
    istart = iend + 1;
    iend = istart + N - 1;
    bsig = zves(istart:iend);
    % seperate position and tension terms

    rhs = bsig + o.bdiagVes.DGB(:,:,k)*o.bdiagVes.IpBen(:,:,k)*bx;
    sig = o.bdiagVes.schur(:,:,k)*rhs;
    rhs = bx + o.dt*o.bdiagVes.GT(:,:,k)*sig;
    pos = o.bdiagVes.IpBen(:,:,k)*rhs;
    valVes((k-1)*3*N+1:3*k*N) = [pos;sig];
    % use schur decomposition to operate preconditioner
  end % k

end % o.solver

zwalls = z(3*N*nv+1:3*N*nv+2*Nbd*nvbd+3*max([(nvbd-1);0]));
% part of z from solid walls
valWalls = o.bdiagWall * zwalls;
zrigid = z(3*N*nv+2*Nbd*nvbd+3*max([(nvbd-1);0])+1:end);
for k=1:nvrd
  valRigid(:,k) = (o.bdiagRigid.U(:,:,k))\...
      (o.bdiagRigid.L(:,:,k)\zrigid((k-1)*(2*Nrd+3)+1:k*(2*Nrd+3)));
end
% this matrix is well-coditioned since it is in the form I + DLP
%if (any(vesicle.viscCont ~= 1) && ...
%      o.confined)
  %strcmp(o.vesves,'implicit') && o.confined)
%  valWalls = valWalls / o.dt;
%end


val = [valVes;valWalls;valRigid(:)];
% stack the two componenets of the preconditioner

end % preconditionerBD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Mat = wallsPrecond(o,walls)
% wallsPrecond(walls) computes the matrix which is the exact inverse of
% the double-layer potential for stokes flow in a bounded domain.  Used
% in the preconditioner for vesicle simulations

Nbd = walls.N;
nvbd = walls.nv;
oc = curve;
[x,y] = oc.getXY(walls.X);
[nory,norx] = oc.getXY(walls.xt);
nory = -nory;
sa = walls.sa;
[cx,cy] = oc.getXY(walls.center);

M11 = zeros(2*Nbd*nvbd,2*Nbd*nvbd);
M12 = zeros(2*Nbd*nvbd,3*(nvbd-1));
M21 = zeros(3*(nvbd-1),2*Nbd*nvbd);
% Allocate space for blocks of matrix that carries the double- layer
% potential, rotlets, and stokeslets to the velocity and the conditions
% in (A4) and (A5) in Rahimian et al.

M11(1:2*Nbd,1:2*Nbd) = M11(1:2*Nbd,1:2*Nbd) + o.wallN0(:,:,1);
jump = - 1/2; 
for k = 1:nvbd
  istart = (k-1)*2*Nbd+1;
  iend = 2*k*Nbd;
  M11(istart:iend,istart:iend) = M11(istart:iend,istart:iend) + ...
      jump*eye(2*Nbd) + o.wallDLP(:,:,k);
end
% Self interaction terms with the jump coming from the double layer
% potential

D = zeros(2*Nbd,2*Nbd);
% temporary space for while we build each off-diagonal component of the
% double-layer potential
for ktar = 1:nvbd % loop over targets
  itar = 2*(ktar-1)*Nbd + 1;
  jtar = 2*ktar*Nbd;
  K = [(1:ktar-1) (ktar+1:nvbd)];
  for ksou = K % loop over all other walls
    isou = 2*(ksou-1)*Nbd + 1;
    jsou = 2*ksou*Nbd;
    for j = 1:Nbd % loop over targets
      rho2 = (x(j,ktar) - x(:,ksou)).^2 + (y(j,ktar)-y(:,ksou)).^2;

      coeff = 1/pi*((x(j,ktar) - x(:,ksou)).*norx(:,ksou) + ...
        (y(j,ktar) - y(:,ksou)).*nory(:,ksou)) .* ...
        sa(:,ksou)./rho2.^2;

      D(j,:) = 2*pi/Nbd*[coeff.*(x(j,ktar) - x(:,ksou)).^2; ...
          coeff.*(x(j,ktar) - x(:,ksou)).*(y(j,ktar) - y(:,ksou))]';
      D(j+Nbd,:) = 2*pi/Nbd*[...
          coeff.*(y(j,ktar) - y(:,ksou)).*(x(j,ktar) - x(:,ksou)); ...
          coeff.*(y(j,ktar) - y(:,ksou)).^2];

      M11(itar:jtar,isou:jsou) = D;
      % place off-diagonal term of wall to wall interactions in
      % appropriate location in top left corner of matrix
    end % j
  end % ksou
end % ktar


for k = 1:nvbd-1
  icol = 3*(k-1)+1;
  istart = 2*k*Nbd+1;
  iend = istart + Nbd - 1;
  M21(icol,istart:iend) = 2*pi/Nbd*sa(:,k+1)';
  M21(icol+2,istart:iend) = 2*pi/Nbd*sa(:,k+1)'.*y(:,k+1)';
  istart = istart + Nbd;
  iend = iend + Nbd;
  M21(icol+1,istart:iend) = 2*pi/Nbd*sa(:,k+1)';
  M21(icol+2,istart:iend) = -2*pi/Nbd*sa(:,k+1)'.*x(:,k+1)';
end % k
% These compute the integral of the density function around each of the
% inner componenents of the geometry

for k = 1:nvbd - 1
  for ktar = 1:nvbd
    rho2 = (x(:,ktar) - cx(k+1)).^2 + (y(:,ktar) - cy(k+1)).^2;
    istart = (ktar-1)*2*Nbd + 1;
    iend = istart + Nbd - 1;

    icol = 3*(k-1)+1;
    M12(istart:iend,icol) = ...
      M12(istart:iend,icol) + ...
      1/4/pi*(-0.5*log(rho2) + (x(:,ktar)-cx(k+1))./rho2.*...
          (x(:,ktar)-cx(k+1)));
    M12(istart + Nbd:iend + Nbd,icol) = ...
      M12(istart + Nbd:iend + Nbd,icol) + ...
      1/4/pi*((x(:,ktar)-cx(k+1))./rho2.*(y(:,ktar)-cy(k+1)));

    icol = 3*(k-1)+2;
    M12(istart:iend,icol) = ...
      M12(istart:iend,icol) + ...
      1/4/pi*((y(:,ktar)-cy(k+1))./rho2.*(x(:,ktar)-cx(k+1)));
    M12(istart + Nbd:iend + Nbd,icol) = ...
      M12(istart + Nbd:iend + Nbd,icol) + ...
      1/4/pi*(-0.5*log(rho2) + (y(:,ktar)-cy(k+1))./rho2.*...
          (y(:,ktar)-cy(k+1)));

    icol = 3*(k-1)+3;
    M12(istart:iend,icol) = ...
      M12(istart:iend,icol) + ...
      (y(:,ktar)-cy(k+1))./rho2;
    M12(istart + Nbd:iend + Nbd,icol) = ...
      M12(istart + Nbd:iend + Nbd,icol) - ...
      (x(:,ktar)-cx(k+1))./rho2;
  end
end
% This is the evaluation of the velocity field due to the stokeslet
% and rotlet terms

M22 = -2*pi*eye(3*(nvbd-1));
% different combinations of the density functions have to multiply to
% 2*pi multiplied by rotlet or stokeslet terms

Mat = ([M11 M12; M21 M22])\eye(2*nvbd*Nbd + 3*(nvbd-1));
% invert the matrix

end % wallsPrecond

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = precoIminusD(o,z,vesicle)
% z = o.precoIminusD(z,vesicle) is the block-diagonal preconditioner
% for the system alpha*I + DLP.  This did not significantly reduce the
% number of GMRES steps, so it is currently not being called

N = vesicle.N;
nv = vesicle.nv;
alpha = 0.5*(1+vesicle.viscCont);
val = zeros(2*N*nv,1);

for k = 1:nv
  istart = (k-1)*2*N + 1;
  iend = k*2*N;
  z(istart:iend) = (alpha(k)*eye(2*N) - o.D(:,:,k))\z(istart:iend);
end

end % precoIminusD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = preconditionerTen(o,z)
% val = preconditionerTen(z) applies the preconditioner to the tension
% term for when we are solving for just the tension and density
% function given a configuration.  Configuration is eliminated using
% the Schur complelent

nv = size(o.bdiagTen,3); % number of vesicles
N = size(o.bdiagTen,1); 
% number of points per vesicle

nvbd = size(o.wallDLP,3); % number of solid walls
Nbd = size(o.wallDLP,1)/2; 
% number of points per solid wall

zves = z(1:N*nv);
% part of z correpsonding to the vesicles
zwall = z(N*nv+1:end-3*(nvbd-1));
% part of z corresponding to the solid walls
zrot = z(end-3*(nvbd-1)+1:end);
% part of z corresonding to rotlets and stokeslets

valVes = zeros(N*nv,1);
for k=1:nv
  valVes((k-1)*N+1:k*N) = o.bdiagTen(:,:,k)*...
    zves((k-1)*N+1:k*N);
end

valWall = o.bdiagWall*[zwall;zrot];

val = [valVes;valWall];
% stack the preconditioned values due to the tension term and
% the solid wall term

end % preconditionerTen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = preconditionerSpectral(o,vesicle,walls,wallsCoarse,z)
% val = preconditionerSpectral(vesicle,walls,wallsCoarse,z) applies the
% block diagonal preconditioner where each block is the inverse
% corresponding to a circle of the same length as the vesicle

N = vesicle.N;
nv = vesicle.nv;
if o.confined
  Nbd = walls.N;
  nvbd = walls.nv;
  valWall = zeros(2*Nbd,nvbd);
end
bx = zeros(2*N,nv);
bsig = zeros(N,nv);
for k = 1:nv
  istart = (k-1)*3*N+1;
  iend = istart + 2*N - 1;
  bx(:,k) = z(istart:iend);
  istart = iend + 1;
  iend = istart + N - 1;
  bsig(:,k) = z(istart:iend);
end
% seperate the right hand side terms due to tension and position

if strcmp(o.solver,'method1')
%  fprintf('Preconditioner not yet developed for this solver\n')
  bschur = bsig - vesicle.surfaceDiv(o.invIplusSB(vesicle,bx));

  sig = o.invSchur11(vesicle,bschur);
  % solve the Schur equation using gmres to find the tension

  x = vesicle.tensionTerm(sig);
  for k = 1:nv
    x(:,k) = o.Galpert(:,:,k)*x(:,k);
  end
  x = bx + o.dt * x;
  % build the right hand side term for solving fot the positon

  x = o.invIplusSB(vesicle,x);
  % solve for the position

elseif strcmp(o.solver,'method2')
  bschur = o.invDST(vesicle,bsig);
  bschur = vesicle.tensionTerm(bschur);
  Sbschur = zeros(2*N,nv);
  for k = 1:nv
    Gbschur(:,k) = o.Galpert(:,:,k) * bschur(:,k);
  end
  bschur = bx + o.dt*Gbschur;
  % right-hand side of the Schur equation involving the Schur complement

  x = o.invSchur22(vesicle,bschur);
  % solve the schur equation using gmres to find the position

  Benx = -vesicle.bendingTerm(x);
  for k = 1:nv
    Benx(:,k) = o.Galpert(:,:,k)*Benx(:,k);
  end
  sigRHS = bsig + vesicle.surfaceDiv(Benx);
  % build the right hand side term for solving for the tension
  sig = o.invDST(vesicle,sigRHS);
  % solve for the tension
end
% solved for position and tension


valVes = zeros(3*N*nv,1);
for k = 1:nv
  istart = (k-1)*3*N+1;
  iend = istart + 3*N - 1;
  valVes(istart:iend) = [x(:,k);sig(:,k)];
end
% stack the position and tension as required by the outermost matvec
% routine

if o.confined
  for k = 1:nvbd
    istart = 3*N*nv + (k-1)*2*Nbd + 1;
    iend = istart + 2*Nbd - 1;
    valWall(:,k) = z(istart:iend);
  end
  valLets = z(iend+1:end);
  nvcycle = 1;
  for k = 1:nvcycle
    [valWall,valLets] = o.vcycle(walls,wallsCoarse,valWall,valLets);
  end
  valWall = valWall(:);
else
  valWall = [];
  valLets = [];
end

val = [valVes;valWall;valLets];
% stack the vesicle data followed by the solid wall data followed by
% the rotlets and stokeslets

end % preconditionerSpectral

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [valWalls,valLets] = vcycle(o,walls,wallsCoarse,zWalls,zLets)
% [valWalls,valLets] = vcycle(walls,wallsCoarse,zWalls,zLets) applies
% two-grid V-cycle to the integral equation defined on the solid
% walls.  Post-smoothing seems to not work correctly so use V(1,0)
% cicyles seems best.  Coarse grid solve is done with the precomputed
% double-layer potential on the coarse grid

N = walls.N;
Ncoarse = wallsCoarse.N; 
nv = walls.nv;
resWalls = zeros(2*N,nv);
valWalls = zeros(2*N,nv);
errWalls = zeros(2*N,nv);
resCoarse = zeros(2*Ncoarse,nv);
errCoarse = zeros(2*Ncoarse,nv);

op = o.op;

valWalls = zeros(2*N,nv);
valLets = zeros(3*(nv-1),1);
% initial guess

for npre = 1:1
  if norm(valWalls) == 0
    valWalls = -2*zWalls;
  else
    valWalls = -2*(zWalls - op.exactStokesDLdiag(walls,o.wallDLP,valWalls));
  end
  % Picard smoother applied to the density function
%  residual = ([zWalls(:);zLets] - ...
%     o.denMatVec([valWalls(:);valLets],walls,[]));
end
% pre-smoothing

residual = ([zWalls(:);zLets] - ...
    o.denMatVec([valWalls(:);valLets],walls,[]));
%clf;
%fprintf('After pre-smoothing\n')

for k = 1:nv
  istart = (k-1)*2*N + 1;
  iend = istart + 2*N - 1;
  resWalls(:,k) = residual(istart:iend);
end
resLets = residual(iend+1:end);
% form the residual

resCoarse = curve.restrict(resWalls,N,Ncoarse);
% restrict the residual

rhsCoarse = [resCoarse(:);resLets];

err = o.bdiagWall*rhsCoarse;
% solve for the coarse grid error

for k = 1:nv
  istart = (k-1)*2*Ncoarse + 1;
  iend = istart + 2*Ncoarse - 1;
  errCoarse(:,k) = err(istart:iend);
end
errLets = err(iend+1:end);

errWalls = curve.prolong(errCoarse,Ncoarse,N);
% prolong the error

valWalls = valWalls + errWalls;
valLets = valLets + errLets;
% add in the error 

for npost = 1:0
  valWalls = -2*(zWalls - ...
      op.exactStokesDLdiag(walls,o.wallDLP,valWalls));
  % Picard smoother applied to the density function
end
% post-smoothing step

end % vcycle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = invSchur22(o,vesicle,z)
% x = invSchur22(vesicle,z) solves the schur complement equation S*x =
% z where S is the operator I + dt*G*Ben - dt*S*Ten*invDST*Div*G*Ben

N = vesicle.N; nv = vesicle.nv;

[z,flag,relres,iter,relresvec] = gmres(@(X) ...
    o.Schur22MatVec(vesicle,X),...
    z(:),[],o.gmresTol,min(o.gmresMaxIter,2*N));

x = zeros(2*N,nv);
for k = 1:nv
  x(:,k) = z((k-1)*2*N+1:k*2*N);
end

end % invSchur22

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = Schur22MatVec(o,vesicle,x)
% val = Schur22MatVec(vesicle,x) applies the Schur complement operator
% to x.  Schur complement is given by I + dt*G*Ben -
% dt*S*Ten*invDST*Div*G*Ben

N = vesicle.N; nv = vesicle.nv;
val = zeros(2*N,nv);
xCols = zeros(2*N,nv);
for k = 1:nv
  xCols(:,k) = x((k-1)*2*N+1:k*2*N); 
end
x = -vesicle.bendingTerm(xCols);
for k = 1:nv
  x(:,k) = o.Galpert(:,:,k) * x(:,k);
end

sig = o.invDST(vesicle,vesicle.surfaceDiv(x));
Tsig = vesicle.tensionTerm(sig);
for k = 1:nv
  val(:,k) = o.Galpert(:,:,k) * Tsig(:,k);
end
val = xCols + o.dt*(x - val);

val = val(:);

end % Schur22MatVec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = precoSchur22(o,vesicle,alpha,scaling,x)
% val = precoSchur22(vesicle,alpha,scaling,x) solves the equation (I +
% dt*kappa*SLP*Ben)*val = x.  This can be used as a preconditioner when
% doing the Schur complement of the (1,1) block in method1 or when
% doing the Schur complement of the (2,2) block in method2
% TODO: This is wrong.  Should do the analytic preconditioner of the
% actual Schur block.

N = vesicle.N; nv = vesicle.nv;
val = zeros(2*N,nv);

for k = 1:nv
  val(:,k) = x((k-1)*2*N+1:k*2*N);
end
val(1:N,:) = fft(val(1:N,:));
val(N+1:2*N,:) = fft(val(N+1:2*N,:));

for k = 1:nv
  val(3:N-1,k) = val(3:N-1,k).*scaling(3:N-1);
  val(N+3:2*N-1,k) = val(N+3:2*N-1,k).*scaling(3:N-1);
end
% all modes with modulus greater than 1 are diagonal

plusOneMode = [val(2,:) ; val(N+2,:)];
minusOneMode = [val(N,:) ; val(2*N,:)];
% need to save these values since the plus and minus one modes
% communicate ie. isn't diagonal

val(N,:) = 1/(4*alpha-1)/(2*alpha-1) *...
  ((2*alpha^2-4*alpha+1)*val(N,:) - ...
  2*1i*alpha^2*val(2*N,:)) + ...
  alpha/(4*alpha-1) * ...
  (plusOneMode(1,:) + 1i*plusOneMode(2,:));
% -1 mode of the x component

val(2*N,:) = 1/(4*alpha-1)/(2*alpha-1) *...
  (2*1i*alpha^2*val(N,:) + ...
  (2*alpha^2-4*alpha+1)*val(2*N,:)) + ...
  alpha/(4*alpha-1) * ...
  (1i*plusOneMode(1,:) - plusOneMode(2,:));
% -1 mode of the y component

val(2,:) = 1/(4*alpha-1)/(2*alpha-1) *...
  ((2*alpha^2-4*alpha+1)*val(2,:) + ...
  2*1i*alpha^2*val(N+2,:)) + ...
  alpha/(4*alpha-1) * ...
  (minusOneMode(1,:) - 1i*minusOneMode(2,:));
% 1 mode of the x component

val(N+2,:) = 1/(4*alpha-1)/(2*alpha-1) *...
  (-2*1i*alpha^2*val(2,:) + ...
  (2*alpha^2-4*alpha+1)*val(N+2,:)) + ...
  alpha/(4*alpha-1) * ...
  (-1i*minusOneMode(1,:) - minusOneMode(2,:));
% 1 mode of the y component


val(1:N,:) = real(ifft(val(1:N,:)));
val(N+1:2*N,:) = real(ifft(val(N+1:end,:)));
val = val(:);

end % precoSchur22

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig = invDST(o,vesicle,z)
% sig = invDST(o,vesicle,z) solves the equation DST*sig = z using pGMRES

N = vesicle.N; nv = vesicle.nv;

[z,flag,relres,iter] = gmres(@(X) o.DSTMatVec(vesicle,X),...
    z(:),[],1e-2*o.gmresTol,min(o.gmresMaxIter,N),...
    @(X) o.precoDST(vesicle,X));
% need a smaller tolerance here to hide that fact that this does not
% make a liner preconditioner.  Krylov stuff breaks down in theory, but
% this shouldn't come up until the error of the outer iteration is
% smaller than the requested tolerance

sig = zeros(N,nv);
for k = 1:nv
  sig(:,k) = z((k-1)*N+1:k*N);
end

end % invDST

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig = DSTMatVec(o,vesicle,sig)
% sig = DSTMatVec(vesicle,sig) applies the operator Div*SLP*Tension to
% sig

N = vesicle.N; nv = vesicle.nv;
sigCols = zeros(N,nv);
for k = 1:nv
  sigCols(:,k) = sig((k-1)*N+1:k*N); 
end

tension = vesicle.tensionTerm(sigCols);
% compute Tension * sig
for k = 1:nv
  tension(:,k) = o.Galpert(:,:,k)*tension(:,k);
end
% apply the single-layer potential
sig = vesicle.surfaceDiv(tension);
% compute the surface divergence
sig = sig(:);

end % DSTMatVec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig = precoDST(o,vesicle,sig)
% sig = precoDST(vesicle,sig) inverts the operator Div*SLP*Tension where
% the geometry is assumed to be a circle that has the same radius as the
% vesicle

N = vesicle.N; nv = vesicle.nv;
rad = vesicle.length/2/pi;
% radius of circle with same length as the vesicles

imodes = 1./[(1:N/2-1)';(N/2:-1:1)'];
for k = 1:nv
  sig((k-1)*N+1:k*N) = fft(sig((k-1)*N+1:k*N));
  sig((k-1)*N+2:k*N) = -4*rad*sig((k-1)*N+2:k*N).*imodes;
  sig((k-1)*N+1:k*N) = ifft(sig((k-1)*N+1:k*N));
end

sig = real(sig);

end % precoDST

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig = invSchur11(o,vesicle,z)
% sig = invSchur11(vesicle,z) solves the schur complement equation S*x
% = z where S is the operator defined by Schur11MatVec using pGMRES

N = vesicle.N; nv = vesicle.nv;

[z,flag,relres,iter,relresvec] = gmres(@(X) ...
    o.Schur11MatVec(vesicle,X),...
    z(:),[],o.gmresTol,min(o.gmresMaxIter,N),...
    @(X) o.precoSchur11(vesicle,X));

sig = zeros(N,nv);
for k = 1:nv
  sig(:,k) = z((k-1)*N+1:k*N);
end

end % invSchur22

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = Schur11MatVec(o,vesicle,sig)
% val = Schur11MatVec(vesicle,sig) applies the operator
% dt*Div*inv(eye+dt*SLP*Bending)*SLP*TEN which is the Schur complement
% of the (1,1) component
% TODO: THIS IS WRONG.  THE OPERATOR TO SOLVE SHOULD BE DIV*SLP*TEN -
% dt*DIV*SLP*BEN*inv(eye+dt*SLP*BEN)*SLP*TEN

N = vesicle.N; nv = vesicle.nv;
sigCols = zeros(N,nv);

for k = 1:nv
  sigCols(:,k) = sig((k-1)*N+1:k*N); 
end

sigCols = vesicle.tensionTerm(sigCols);
% This changes the size of sigCols to 2*N x nv.  This is no longer the
% tension, but the tension operator applied to the tension
for k = 1:nv
  sigCols(:,k) = o.Galpert(:,:,k) * sigCols(:,k);
end

sigCols = o.invIplusSB(vesicle,sigCols);
sigCols = o.dt*vesicle.surfaceDiv(sigCols);

val = sigCols(:);

end % Schur11MatVec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig = precoSchur11(o,vesicle,sig)
% sig = precoSchur11(vesicle,sig) applies the spectral preconditioner
% of the Schur complement matrix when using the Schur complement of
% method1

N = vesicle.N; nv = vesicle.nv;
modes = [(0:N/2-1)';(-N/2:-1)'];
rad = vesicle.length/2/pi;
% radius of circle with same length as the vesicles

sigColh = zeros(N,nv);
for k = 1:nv
  sigColh(:,k) = sig((k-1)*N+1:k*N);
end
sigColh = fft(sigColh);

alpha = -o.dt*vesicle.kappa/8/rad^3;

scaling = 1./(-1/8*(abs(modes+1)./(1-2*alpha*abs(modes+1).^3) + ...
         abs(modes-1)./(1-2*alpha*abs(modes-1).^3)));
scaling(1) = 1;
scaling(2) = -4*(1-16*alpha);
scaling(3) = 1./(1.25e-1/(2*alpha-1) - 3.75e-1/(1-54*alpha));
scaling(N) = -4*(1-16*alpha);
scaling(N-1) = 1./(1.25e-1/(2*alpha-1) - 3.75e-1/(1-54*alpha));

for k = 1:nv
  sigColh(:,k) = sigColh(:,k) .* scaling;
end

sigColh = real(ifft(sigColh));

sig = rad*sigColh(:)/o.dt;

end % precoSchur11

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = invIplusSB(o,vesicle,z)
% x = invIplusSB(vesicle,z) solves the equation (I + dt*SLP*Bending)x =
% z using pGMRES 

N = vesicle.N; nv = vesicle.nv;

rad = vesicle.length/2/pi;
% radius of circle with same length as the vesicles
alpha = -o.dt*vesicle.kappa/8/rad^3;
% re-occuring constant
scaling = 1./(1-2*alpha*[(0:N/2-1) (N/2:-1:1)]'.^3);
% scaling that happens to the modes from modulus greater than or equal
% to 2 in the preconditioner
[z,flag,relres,iter] = gmres(@(X) o.IplusSBMatVec(vesicle,X),...
    z(:),[],1e-2*o.gmresTol,min(o.gmresMaxIter,N),...
    @(X) o.precoIplusSB(vesicle,alpha,scaling,X));

x = zeros(2*N,nv);
for k = 1:nv
  x(:,k) = z(2*(k-1)*N+1:k*2*N);
end

end % invIplusSB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = IplusSBMatVec(o,vesicle,x)
% x = IplusSBMatVec(vesicle,x) applies the operator identity + dt *
% SLP*Bending where \kappa is absored into the term Bending

N = vesicle.N; nv = vesicle.nv;
xCols = zeros(2*N,nv);
for k = 1:nv
  xCols(:,k) = x(2*(k-1)*N+1:k*2*N); 
end

xCols = -vesicle.bendingTerm(xCols);
% compute Bending * x
for k = 1:nv
  xCols(:,k) = o.Galpert(:,:,k)*xCols(:,k);
end
x = x + o.dt*xCols(:);

end % IplusSBMatVec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = precoIplusSB(o,vesicle,alpha,scaling,x)
% val = precoIplusSB(vesicle,alpha,scaling,x) solves the equation (I +
% dt*kappa*SLP*Ben)*val = x analytically for a circle

N = vesicle.N; nv = vesicle.nv;
val = zeros(2*N,nv);

for k = 1:nv
  val(:,k) = x((k-1)*2*N+1:k*2*N);
end
val(1:N,:) = fft(val(1:N,:));
val(N+1:2*N,:) = fft(val(N+1:2*N,:));

for k = 1:nv
  val(3:N-1,k) = val(3:N-1,k).*scaling(3:N-1);
  val(N+3:2*N-1,k) = val(N+3:2*N-1,k).*scaling(3:N-1);
end
% all modes with modulus greater than 1 are diagonal

plusOneMode = [val(2,:) ; val(N+2,:)];
minusOneMode = [val(N,:) ; val(2*N,:)];
% need to save these values since the plus and minus one modes
% communicate ie. isn't diagonal

val(N,:) = 1/(4*alpha-1)/(2*alpha-1) *...
  ((2*alpha^2-4*alpha+1)*val(N,:) - ...
  2*1i*alpha^2*val(2*N,:)) + ...
  alpha/(4*alpha-1) * ...
  (plusOneMode(1,:) + 1i*plusOneMode(2,:));
% -1 mode of the x component

val(2*N,:) = 1/(4*alpha-1)/(2*alpha-1) *...
  (2*1i*alpha^2*val(N,:) + ...
  (2*alpha^2-4*alpha+1)*val(2*N,:)) + ...
  alpha/(4*alpha-1) * ...
  (1i*plusOneMode(1,:) - plusOneMode(2,:));
% -1 mode of the y component

val(2,:) = 1/(4*alpha-1)/(2*alpha-1) *...
  ((2*alpha^2-4*alpha+1)*val(2,:) + ...
  2*1i*alpha^2*val(N+2,:)) + ...
  alpha/(4*alpha-1) * ...
  (minusOneMode(1,:) - 1i*minusOneMode(2,:));
% 1 mode of the x component

val(N+2,:) = 1/(4*alpha-1)/(2*alpha-1) *...
  (-2*1i*alpha^2*val(2,:) + ...
  (2*alpha^2-4*alpha+1)*val(N+2,:)) + ...
  alpha/(4*alpha-1) * ...
  (-1i*minusOneMode(1,:) - minusOneMode(2,:));
% 1 mode of the y component

val(1:N,:) = real(ifft(val(1:N,:)));
val(N+1:2*N,:) = real(ifft(val(N+1:end,:)));
val = val(:);

end % precoIplusSB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [time,pos] = noVesTracers(o,walls,Xtra)
% [time,pos] = noVesTracers(walls,Xtra) computes the position of the
% tracers based only on the density function defined on the solid walls

op = o.op;
[eta,RS,~] = walls.computeEta(o);

tracers.N = numel(Xtra)/2;
tracers.nv = 1;
tracers.X = Xtra;
% build a structure for tracers which mimics the class capsules.
% Can't use capsules unless there is an even number of tracer
% points (because of FFTs), so we don't want to build the class
% as we do for vesicle and walls

if ~o.near
  if ~o.fmmDLP
    kernel = @o.op.exactStokesDL;
  else
    kernel = @o.op.exactStokesDLfmm;
  end
  [~,Fwall2Tra] = kernel(walls,eta,Xtra,1:walls.nv);
else
  DLP = o.wallDLP;
  jump = -1/2;
  for k = 1:walls.nv
    DLP(:,:,k) = DLP(:,:,k) + jump*eye(2*walls.N);
  end
end

tfinal = 50;
odefun = @(t,z) o.tracersVel4RK(t,tfinal,z,walls,tracers,DLP,eta,RS);
options.RelTol = 1e-6;
options.AbsTol = 1e-9;
[time,pos] = ode45(odefun,linspace(0,tfinal,501),Xtra,options);
fprintf('\n')

end % noVesTracers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vel = tracersVel4RK(o,t,T,Xtra,walls,tracers,DLP,eta,RS)
% vel = tracersVel4RK(Xtra,walls,tracers,DLP,eta,RS) evaluates the
% right-hand-side for doing time stepping of the tracers when using RK45
message = ['ode45 ' num2str(t/T*100,'%04.1f') ' %% completed '];
nmess = numel(message);
fprintf(repmat('\b',1,nmess));
fprintf(message);

op = o.op;

tracers.X = Xtra;
[~,NearW2T] = walls.getZone(tracers,2);
if ~o.fmmDLP
  kernel = @op.exactStokesDL;
else
  kernel = @op.exactStokesDLfmm;
end
if ~o.near
  [~,Fwall2Tra] = kernel(walls,eta,Xtra,1:walls.nv);
else
  kernelDirect = kernel;
  DLPfun = @(X) op.exactStokesDLdiag(walls,DLP,X);
  Fwall2Tra = op.nearSingInt(walls,eta,DLPfun,DLPfun,...
      NearW2T,kernel,kernelDirect,tracers,false);
end
 

FLets2Tra = zeros(size(Xtra));
for k = 2:walls.nv
  stokeslet = RS(1:2,k);
  rotlet = RS(3,k);
  FLets2Tra = FLets2Tra + o.RSlets(Xtra,walls.center(:,k),...
    stokeslet,rotlet);
end

vel = Fwall2Tra + FLets2Tra;
ind = find(Xtra(1:end/2) > 8);
vel(ind) = 0;
vel(ind + tracers.N) = 0;

end % tracersVel4RK

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vel = tracersVel(o,X,sigma,u,kappa,viscCont,...
    walls,eta,RS,Xtra,NearW2T)
% vel = tracersVel(X,sigma,kappa,walls,eta,RS,Xtra) computes the
% velocity at the set of points Xtra due to the vesicle, viscosity
% contrast, solid walls, and stokeslets and rotlets.  The particles Xtra
% are treated passively

vesicle = capsules(X,sigma,[],kappa,viscCont,o.antiAlias);
tracers.N = numel(Xtra)/2;
tracers.nv = 1;
tracers.X = Xtra;
% build a structure for tracers which mimics the class capsules.
% Can't use capsules unless there is an even number of tracer
% points (because of FFTs), so we don't want to build the class
% as we do for vesicle and walls

[~,NearV2T] = vesicle.getZone(tracers,2);
%[~,NearV2T] = vesicle.getZoneNew(tracers,2);
% build near-singular integration structures for vesicle to tracer and
% wall to tracer interactions

[~,NearW2T] = walls.getZone(tracers,2);
% save WallGetZoneCase2_256x1024_TAdap NearW2T

Ntra = size(Xtra,1)/2; % Number of tracers
N = size(X,1)/2; % Number of points
nv = size(X,2); % Number of vesicles
if o.confined
  Nbd = walls.N; % Number of points on walls
  nvbd = walls.nv; % Number of components to walls
else
  Nbd = 0;
  nvbd = 0;
end
vel = zeros(2*Ntra,1);

f = vesicle.tracJump(X,sigma);
if o.gravity
fgrav = o.gravityForce(X,vesicle,o.gCont);
f = f+fgrav;
end
% compute traction jump

op = poten(N);

%% START OF VELOCITY DUE TO VESICLES
if ~o.fmm
  kernel = @op.exactStokesSL;
else
  kernel = @op.exactStokesSLfmm;
end
if ~o.near
  [~,Fves2Tra] = kernel(vesicle,f,[],Xtra,(1:nv));
  % velocity due to vesicle with the N-point trapezoid rule
else
  o.Galpert = op.stokesSLmatrix(vesicle);
  SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);
  SLPtrap = SLP;
  % need single-layer potential matrix for doing near-singular
  % integration
  kernelDirect = kernel;
  Fves2Tra = op.nearSingInt(vesicle,f,SLP,SLPtrap,...
     NearV2T,kernel,kernelDirect,tracers,false);
  % Evaulate velocity due to vesicle on tracers with near-singular
  % integration
end
% END OF VELOCITY DUE TO VESICLES

% START OF VELOCITY DUE TO SOLID WALLS OR BACKGROUND VELOCITY
if o.confined
  if ~o.fmmDLP
    kernel = @op.exactStokesDL;
  else
    kernel = @op.exactStokesDLfmm;
  end
  if ~o.near
    [~,Fwall2Tra] = kernel(walls,eta,[],Xtra,1:nvbd);
  else
    DLP = o.wallDLP;
    jump = -1/2;
    for k = 1:nvbd
      DLP(:,:,k) = DLP(:,:,k) + jump*eye(2*Nbd);
    end
    DLPfun = @(X) op.exactStokesDLdiag(walls,DLP,X);
    kernelDirect = kernel;
    Fwall2Tra = op.nearSingInt(walls,eta,DLPfun,DLPfun,...
        NearW2T,kernel,kernelDirect,tracers,false);
  end
  % compute the velocity on the tracers due to the
  % solid walls

  FLets2Tra = zeros(2*Ntra,1);
  for k = 2:nvbd
    stokeslet = RS(1:2,k);
    rotlet = RS(3,k);
    FLets2Tra = FLets2Tra + o.RSlets(Xtra,walls.center(:,k),...
      stokeslet,rotlet);
  end
  % velocity due to the stokeslets and rotlets
else
  Fwall2Tra = o.farField(Xtra);
  FLets2Tra = zeros(2*Ntra,1);
  % if no solid walls, velocity field comes from background
  % velocity.  No stokeslet or rotlets so they induce no
  % velocity
end
% END OF VELOCITY DUE TO SOLID WALLS OR BACKGROUND VELOCITY


% START OF VELOCITY DUE TO VISCOSITY CONTRAST
if any(viscCont ~= 1)
  if ~o.fmmDLP
    kernel = @op.exactStokesDL;
  else
    kernel = @op.exactStokesDLfmm;
  end

  if o.near
    o.D = op.stokesDLmatrix(vesicle);
    % need double-layer potential matrix for doing near-singular
    % integration
    DLP = o.D;
    % Get the DLP for inside&outside
    DLP_out = DLP;
    DLP_in  = DLP;
    for k=1:nv
      % if tracer is outside a vesicle, jump:  
      jump_out = 1/2*(1-viscCont(k));
      % if tracer is inside a vesicle, jump:
      jump_in  = -jump_out;
      
      % Based on jump, DLP is computed for inside/outside
      DLP_out(:,:,k) = DLP_out(:,:,k) + jump_out*eye(2*N);
      DLP_in (:,:,k) = DLP_in (:,:,k) + jump_in *eye(2*N);
      
      DLPfun_out = @(X) op.exactStokesDLdiag(vesicle,DLP_out,X);
      DLPfun_in  = @(X) op.exactStokesDLdiag(vesicle,DLP_in ,X);
      
      % TODO: THIS JUMP NEEDS TO BE NEGATED FOR TRACERS INSIDE
      % THE VESICLE BECAUSE OF THE DISCONTINUITY OF THE DOUBLE-
      % LAYER POTENTIAL.
    end
  end

  if ~o.near
    [~,Fvisc2Tra] = kernel(vesicle,u,[],Xtra,(1:nv));
    % velocity due to vesicle's viscosity contrast on 
    % tracers with the N-point trapezoid rule
  else
    kernelDirect = kernel;
    % ---------------------------------------------------------------------
    % Not correct:
%     Fvisc2Tra = op.nearSingInt(vesicle,u,DLPfun,DLPfun,...
%       NearV2T,kernel,kernelDirect,tracers,false);
    % ---------------------------------------------------------------------
    
    % Initialize Fvisc2Tra
    Fvisc2Tra = zeros(2*Ntra,1);
    % Determine if the points inside or outside
    InOut = vesicle.sortPts(Xtra,o.fmm,NearV2T);
    InOut_ext = [InOut;InOut]; % since size(InOut) = size(Xtra)/2 (only for one component) 
    
    % OUTSIDE: build tracers and compute Fvisc2Tra
    %----------------------------------------------------------------------
    Xtra_out = Xtra(InOut_ext == 0);
    tracers.N = numel(Xtra_out)/2;
    tracers.X = Xtra_out;
    [~,NearV2T_out] = vesicle.getZone(tracers,2);
    
    Fvisc2Tra(InOut_ext==0) = op.nearSingInt(vesicle,u,DLPfun_out,DLPfun_out,...
      NearV2T_out,kernel,kernelDirect,tracers,false);
    % Evaulate velocity due to vesicle's viscosity contrast on 
    % tracers with near-singular integration
    
    % INSIDE: build tracers and compute Fvisc2Tra
    %----------------------------------------------------------------------
    Xtra_in = Xtra(InOut_ext == 1);
    tracers.N = numel(Xtra_in)/2;
    tracers.X = Xtra_in;
    [~,NearV2T_in] = vesicle.getZone(tracers,2);
    
    Fvisc2Tra(InOut_ext==1) = op.nearSingInt(vesicle,u,DLPfun_in,DLPfun_in,...
      NearV2T_in,kernel,kernelDirect,tracers,false);
    % Evaulate velocity due to vesicle's viscosity contrast on 
    % tracers with near-singular integration
    
    % NOTE: vesicle.getZone computed before might be used here.
    
    % Back to original
    tracers.N = numel(Xtra)/2;
    tracers.X = Xtra;
  end

else
  Fvisc2Tra = zeros(2*Ntra,1);
  % no velocity due to viscosity contrast
end
% END OF VELOCITY DUE TO VISCOSITY CONTRAST

vel = Fves2Tra + Fwall2Tra + FLets2Tra + Fvisc2Tra;
% velocity is the summy of velocities due to vesicles, solid walls,
% rotlet and stokeslets, and viscosity contrast

%vel = Fwall2Tra + FLets2Tra;
% DEBUG: FOR TESTING TRACERS WITH NO VESICLES
%r2 = Xtra(1).^2 + Xtra(2).^2;
%speed = -1/3 + 400/3/r2;
%velTrue = speed*[-Xtra(2);Xtra(1)];
% This is the true velocity for the couette apparatus when both
% boundaries are centered at the origin, the inner boundary rotates
% once every 2*pi time units and the outer boundary is fixed

end % tracersVel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vInf = bgFlow(o,X,varargin);
% vInf = bgFlow(X,varargin) computes the velocity field at the points X.
% Default flow is shear with magnitude 1.  Other flows are relaxation,
% extensional, parabolic, Taylor-Green, cubic, invParabolic.  Some flows
% have an option to input a strength of the flow (ie. shear rate,
% extensional rate, etc).  Flows are given by:
%     cubic:          (y^3,x)
%     relaxation:     (0,0)
%     extensional:    (-x,y)
%     parabolic:      (k(1-y/r)^2,0)
%     invParabolic:   (ky^2,0)
%     rotate:         (2y,-x)
%     taylorGreen:    (sin(x)cos(y),-cos(x)sin(y))
%     shear:          (ky,0)
%     choke:          poeusille-like flow at intake and outtake
%     doublechoke:    same as choke
%     couette:        rotating boundary
%     doubleCouette   two rotating boundaries
%     quadCouette     four rotation boundaries
%     doubleFlower    same as doubleCouette
%     figureEight     tangential boundary condition
%     tube            shear boundary condition
% k is a 'flow strength', and r controls the width of the parabola
% in the parabolic flow

N = size(X,1)/2; % number of points per vesicle
nv = size(X,2); % number of vesicles

oc = curve;
[x,y] = oc.getXY(X);
% Separate out x and y coordinates

if any(strcmp(varargin,'cubic'))
  vInf = [y.^3;zeros(N,nv)];

elseif any(strcmp(varargin,'relaxation'))
  vInf = zeros(2*N,nv); 

elseif any(strcmp(varargin,'extensional'))
  vInf = [-x;y];

elseif any(strcmp(varargin,'parabolic'))
  R = find(strcmp(varargin,'R'));
  if isempty(R)
%    R = 100;
      R = 20;
  else
    R = varargin{R+1};
  end
  UM = find(strcmp(varargin,'Umax'));
  if isempty(UM)
%    UM = 1e5; % default value for strength of flow
    UM = 1;
  else
    UM = varargin{UM+1};
  end
  vInf = [UM*(1-(y/R).^2);zeros(N,nv)];

elseif any(strcmp(varargin,'invParabolic'))
  k = find(strcmp(varargin,'k'));
  if isempty(k)
    k = 1; % default value for strength of flow
  else
    k = varargin{k+1};
  end
  vInf = [k*y.^2;zeros(N,nv)];

elseif any(strcmp(varargin,'rotate'))
  vInf = [2*y;-x];

elseif any(strcmp(varargin,'taylorGreen'))
  vInf = [sin(x).*cos(y);-cos(x).*sin(y)];

elseif any(strcmp(varargin,'shear'))
  k = find(strcmp(varargin,'k'));
  if isempty(k)
    k = 1; % default value for strength of flow
  else
    k = varargin{k+1};
  end
  vInf = [k*y;zeros(N,nv)];

elseif (any(strcmp(varargin,'choke')) || ...
      any(strcmp(varargin,'doublechoke')) || ...
      any(strcmp(varargin,'choke2')) || ...
      any(strcmp(varargin,'tube')))
  vInf = zeros(2*N,nv);
  ind = abs(x)>7;
  vx = exp(1./((y(ind)/max(y)).^2-1))/exp(-1);
  % typical mollifer so that velocity decays smoothly to 0
  vx(vx==Inf) = 0;
  vInf(ind,:) = vx;

elseif any(strcmp(varargin,'chokeerf'))
  vInf = zeros(2*N,nv);
  ind = abs(x)>12;
  vx = exp(1./((y(ind)/max(y)).^2-1))/exp(-1);
  % typical mollifer so that velocity decays smoothly to 0
  vx(vx==Inf) = 0;
  vInf(ind,:) = vx;

elseif any(strcmp(varargin,'stenosis'))
  vInf = [y+max(y);zeros(N,nv)];
  indx = (abs(x)<=pi);
  indy = (y<(max(y)-0.1));
  ind = and(indx,indy);
  vInf(ind,:) = 0;
  
elseif any(strcmp(varargin,'box'))
  vInf = zeros(2*N,nv);

elseif any(strcmp(varargin,'couette'));
  vInf = [zeros(2*N,1) 1*[-y(:,2)+mean(y(:,2));x(:,2)-mean(x(:,2))]];

elseif any(strcmp(varargin,'couette10'));
  vInf = [zeros(2*N,1) 10*[-y(:,2)+mean(y(:,2));x(:,2)-mean(x(:,2))]];
  
elseif any(strcmp(varargin,'couette100'));
  vInf = [zeros(2*N,1) 100*[-y(:,2)+mean(y(:,2));x(:,2)-mean(x(:,2))]];
  
elseif (any(strcmp(varargin,'doubleCouette')) || ...
      any(strcmp(varargin,'doubleFlower')));
  vInf = [zeros(2*N,1) 1*[-y(:,2)+mean(y(:,2));x(:,2)-mean(x(:,2))] ...
      -[y(:,3)-mean(y(:,3));-x(:,3)+mean(x(:,3))]];

elseif (any(strcmp(varargin,'quadCouette')));
  vInf = [zeros(2*N,1) ...
      -[y(:,2)-mean(y(:,2));-x(:,2)+mean(x(:,2))] ...
      -[y(:,3)-mean(y(:,3));-x(:,3)+mean(x(:,3))] ...
      +[y(:,4)-mean(y(:,4));-x(:,4)+mean(x(:,4))] ...
      +[y(:,5)-mean(y(:,5));-x(:,5)+mean(x(:,5))]];

elseif any(strcmp(varargin,'cylinder'))
%  theta = (0:N-1)'*2*pi/N;
%  vInf = [cos(10*theta);sin(2*theta)];
  vInf = 1*[-y+mean(y);x-mean(x)];

elseif any(strcmp(varargin,'figureEight'))
  oc = curve;
  [~,vInf,~] = oc.diffProp([x;y]);
%  vInf(1:end/2,:) = 3*vInf(1:end/2,:);

  sup = find(abs(x)<=1 & y>0);
  sdown = find(abs(x)<=1 & y<0);
  omega = linspace(-1,1,numel(sup)+2)';
  omega = omega(2:end-1);
  mollifier = 4*exp(1./(omega.^2 - 1))+1;
  vInf(sup,:) = vInf(sup,:) .* mollifier;
  vInf(sdown,:) = vInf(sdown,:) .* mollifier;
  % increase the velocity in a smooth fashion near the middle
  % of the solid walls

elseif any(strcmp(varargin,'shear'))
  vInf = [y;zeros(N,nv)];

elseif any(strcmp(varargin,'diffuser'));
  vInf = zeros(2*N,nv);
  ind = abs(x(:,1))>9;
  vx = exp(1./((y(ind,1)/max(y(ind,1))).^2-1))/exp(-1);
  % typical mollifer so that velocity decays smoothly to 0
  vx(vx==Inf) = 0;
  vInf(ind,1) = vx;

elseif any(strcmp(varargin,'microfluidic'));
  oc = curve;
  [~,tangent,~] = oc.diffProp(X); 
  vInf = tangent;
  vInf(:,1) = 0*vInf(:,1);
  vInf(:,2) = +1*vInf(:,2);
  vInf(:,3) = -1*vInf(:,3);
  vInf(:,4) = -1*vInf(:,4);
  vInf(:,5) = +1*vInf(:,5);

else
  vInf = [y;zeros(N,nv)];
  % default flow is shear
end

speed = varargin{2};
% speed of the background velocity
vInf = vInf * speed;

end % bgFlow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GLpts = gaussLobatto(o,orderGL)
% GLpts = gaussLobatto(orderGL) loads the Guass- Lobatto points for the
% desired number of points 

if (orderGL == 2 || orderGL == -2)
  GLpts = [-1 1];
elseif (orderGL == 3 || orderGL == -3)
  GLpts = [-1 0 1];
elseif orderGL == 4
  GLpts = [-1 -0.447213595499958 0.447213595499958 1];
elseif orderGL == -4
  GLpts = [-1 -0.333333333333333 0.333333333333333 1];
elseif orderGL == 5
  GLpts = [-1 -0.654653670707977 0 0.654653670707977 1];
elseif orderGL == -5
  GLpts = [-1 -0.5 0 0.5 1];
elseif orderGL == 6
  GLpts = [-1 -0.765055323929465 -0.285231516480645 ...
      0.285231516480645 0.765055323929465 1];
elseif orderGL == -6
  GLpts = [-1 -0.6 -0.2 0.2 0.6 1];
elseif orderGL == 7
  GLpts = [-1 -0.830223896278567 -0.468848793470714 ...
    0 0.468848793470714 0.830223896278567 1];
elseif orderGL == -7
  GLpts = [-1 -0.6666666666666666  -0.3333333333333333 ...
    0 0.3333333333333333 0.6666666666666667 1];
else
  fprintf('**************************************\n')
  fprintf('NO GAUSS-LOBATTO POINTS FOR THIS ORDER\n')
  fprintf('**************************************\n')
  pause
end

end % gaussLobatto

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fint = lobattoInt(o,f)
% fint = lobattoInt(f) returns the integral of f from [-1,t] where t
% ranges from -1 to 1.  If orderGL > 0, the interior points are
% Gauss-Lobato points.  Otherwise, they are equispaced points.  The
% quadrature rules were found by first interpolating f with a
% polynomial and then integrating the polynomial exactly.  All
% coefficients are precomputed f is of size (N,nv,order) where N is the
% number of points that need to be integrated, nv is a second index
% which in this case corresponds to number of vesicles, and order is
% the number of time steps

nv = size(f,2);
% number of vesicles
fint = zeros(size(f));
% initialize the integral to zero

if (o.orderGL == 2 || o.orderGL == -2)
  t = [-1 1];

  for n = 1:2
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*0.25*(-f(:,k,1)+f(:,k,2));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*0.5*(f(:,k,1)+f(:,k,2));

      fint(:,k,n) = fint(:,k,n) + ...
          (0.75*f(:,k,1) + 0.25*f(:,k,2));
    end
  end
  % order 1 Gauss-Lobatto or equispaced


elseif (o.orderGL == 3 || o.orderGL == -3)
  t = [-1 0 1];

  for n = 1:3
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(0.1666666666666667*(f(:,k,1)+f(:,k,3)) - ...
          0.3333333333333333*f(:,k,2));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(-0.25*(f(:,k,1)-f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*f(:,k,2);

      fint(:,k,n) = fint(:,k,n) + ...
          0.4166666666666667*f(:,k,1) + ...
          0.6666666666666667*f(:,k,2) - ...
          0.0833333333333333*f(:,k,3);
    end
  end
  % order 2 Gauss-Lobatto or equispaced


elseif (o.orderGL == 4)
  t = [-1 -0.447213595499958 0.447213595499958 1];

  for n = 1:4
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(-0.15625*(f(:,k,1)-f(:,k,4)) + ...
          0.3493856214843422*(f(:,k,2)-f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(0.2083333333333333*(f(:,k,1)+f(:,k,4)) - ...
          0.2083333333333333*(f(:,k,2)+f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(0.0625*(f(:,k,1)-f(:,k,4)) - ...
          0.6987712429686845*(f(:,k,2)-f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*(-0.125*(f(:,k,1)+f(:,k,4)) + ...
          0.625*(f(:,k,2)+f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          (0.1770833333333333*f(:,k,1) + ...
          0.7660522881510090*f(:,k,2) + ...
          0.0672810451823247*f(:,k,3) - ...
          0.0104166666666666*f(:,k,4));
    end
  end
  % order 3 Gauss-Lobatto

elseif (o.orderGL == -4)
  t = [-1 -0.333333333333333 0.333333333333333 1];

  for n = 1:4
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(-0.140625*(f(:,k,1)-f(:,k,4)) + ...
          0.421875*(f(:,k,2)-f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(0.1875*(f(:,k,1)+f(:,k,4)) - ...
          0.1875*(f(:,k,2)+f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(0.03125*(f(:,k,1)-f(:,k,4)) - ...
          0.84375*(f(:,k,2)-f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*(-0.0625*(f(:,k,1)+f(:,k,4)) + ...
          0.5625*(f(:,k,2)+f(:,k,3)));

      fint(:,k,n) = fint(:,k,n) + ...
          (0.234375*f(:,k,1) + ...
          0.796875*f(:,k,2) - ...
          0.046875*f(:,k,3) + ...
          0.015625*f(:,k,4));
    end
  end
  % order 3 equispaced

elseif (o.orderGL == 5)
  t = [-1 -0.654653670707977 0 0.654653670707977 1];

  for n = 1:5
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^5*(0.175*(f(:,k,1)+f(:,k,5)) - ...
          0.4083333333333333*(f(:,k,2)+f(:,k,4)) + ...
          0.4666666666666666*f(:,k,3));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(-0.21875*(f(:,k,1)-f(:,k,5)) + ...
          0.3341461444238633*(f(:,k,2)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(-0.125*(f(:,k,1)+f(:,k,5)) + ...
          0.6805555555555555*(f(:,k,2)+f(:,k,4)) - ...
          1.1111111111111111*f(:,k,3));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(0.1875*(f(:,k,1)-f(:,k,5)) - ...
          0.6682922888477265*(f(:,k,2)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*f(:,k,3);

      fint(:,k,n) = fint(:,k,n) + ...
          0.08125*f(:,k,1) + ...
          0.6063683666460855*f(:,k,2) + ...
          0.3555555555555555*f(:,k,3) - ...
          0.0619239222016411*f(:,k,4) + ...
          0.01875*f(:,k,5);
    end
  end
  % order 4 Gauss-Lobatto

elseif (o.orderGL == -5)
  t = [-1 -0.5 0 0.5 1];

  for n = 1:5
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^5*(0.1333333333333333*(f(:,k,1)+f(:,k,5)) - ...
          0.53333333333333333*(f(:,k,2)+f(:,k,4)) + ...
          0.8*f(:,k,3));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(-0.16666666666666667*(f(:,k,1)-f(:,k,5)) + ...
          0.3333333333333333*(f(:,k,2)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(-0.0555555555555556*(f(:,k,1)+f(:,k,5)) + ...
          0.8888888888888889*(f(:,k,2)+f(:,k,4)) - ...
          1.6666666666666667*f(:,k,3));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(0.0833333333333333*(f(:,k,1)-f(:,k,5)) - ...
          0.6666666666666667*(f(:,k,2)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*f(:,k,3);

      fint(:,k,n) = fint(:,k,n) + ...
          0.1611111111111111*f(:,k,1) + ...
          0.6888888888888889*f(:,k,2) + ...
          0.1333333333333333*f(:,k,3) + ...
          0.0222222222222222*f(:,k,4) - ...
          0.0055555555555556*f(:,k,5);
    end
  end
  % order 4 equi-spaced points

elseif (o.orderGL == 6)
  t = [-1 -0.765055323929465 -0.285231516480645 ...
      0.285231516480645 0.765055323929465 1];

  for n = 1:6
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^6*(-0.21875*(f(:,k,1)-f(:,k,6)) + ...
          0.5212094304495727*(f(:,k,2)-f(:,k,5)) - ...
          0.6310805056491861*(f(:,k,3)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^5*(0.2625*(f(:,k,1)+f(:,k,6)) - ...
          0.4785048595772274*(f(:,k,2)+f(:,k,5)) + ...
          0.2160048595772274*(f(:,k,3)+f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(0.21875*(f(:,k,1)-f(:,k,6)) - ...
          0.8454202131918329*(f(:,k,2)-f(:,k,5)) + ...
          1.500687022042463*(f(:,k,3)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(-0.2916666666666666*(f(:,k,1)+f(:,k,6)) + ...
          0.8623909800799931*(f(:,k,2)+f(:,k,5)) - ...
          0.5707243134133265*(f(:,k,3)+f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(-0.03125*(f(:,k,1)-f(:,k,6)) + ...
          0.1272121350349483*(f(:,k,2)-f(:,k,5)) - ...
          1.108132527137368*(f(:,k,3)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*(0.0625*(f(:,k,1)+f(:,k,6)) - ...
          0.1946486423538423*(f(:,k,2)+f(:,k,5)) + ...
          0.6321486423538425*(f(:,k,3)+f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          0.06458333333333333*f(:,k,1) + ...
          0.3862361258562357*f(:,k,2) + ...
          0.5159551992618339*f(:,k,3) + ...
          0.03890317777365201*f(:,k,4) - ...
          0.007761169558388661*f(:,k,5) + ...
          0.002083333333333333*f(:,k,6);
    end
  end
  % order 5 Gauss-Lobatto

elseif (o.orderGL == -6)
  t = [-1 -0.6 -0.2 0.2 0.6 1];

  for n = 1:6
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^6*(-0.1356336805555556*(f(:,k,1)-f(:,k,6)) + ...
          0.6781684027777778*(f(:,k,2)-f(:,k,5)) - ...
          1.3563368055555556*(f(:,k,3)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^5*(0.1627604166666667*(f(:,k,1)+f(:,k,6)) - ...
          0.48828125*(f(:,k,2)+f(:,k,5)) + ...
          0.3255208333333333*(f(:,k,3)+f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(0.0813802083333333*(f(:,k,1)-f(:,k,6)) - ...
          1.0579427083333333*(f(:,k,2)-f(:,k,5)) + ...
          2.7669270833333333*(f(:,k,3)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(-0.108506944444444*(f(:,k,1)+f(:,k,6)) + ...
          0.8463541666666667*(f(:,k,2)+f(:,k,5)) - ...
          0.7378472222222222*(f(:,k,3)+f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(-0.005859375*(f(:,k,1)-f(:,k,6)) + ...
          0.0813802083333333*(f(:,k,2)-f(:,k,5)) - ...
          1.46484375*(f(:,k,3)-f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*(0.01171875*(f(:,k,1)+f(:,k,6)) - ...
          0.09765625*(f(:,k,2)+f(:,k,5)) + ...
          0.5859375*(f(:,k,3)+f(:,k,4)));

      fint(:,k,n) = fint(:,k,n) + ...
          0.1260850694444444*f(:,k,1) + ...
          0.5588107638888889*f(:,k,2) + ...
          0.2278645833333333*f(:,k,3) + ...
          0.1193576388888889*f(:,k,4) - ...
          0.0379774305555556*f(:,k,5) + ...
          0.005859375*f(:,k,6);
    end
  end
  % order 5 equi-spaced points

elseif (o.orderGL == 7)
  t = [-1.000000000000000 -0.830223896278567 ...
      -0.468848793470714 0 0.468848793470714 ...
      0.830223896278567 1.000000000000000];

  for n = 1:7
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^7*(0.29464285714285718*(f(:,k,1)+f(:,k,7)) - ...
          0.71040995801242226*(f(:,k,2)+f(:,k,6)) + ...
          0.88719567229813686*(f(:,k,3)+f(:,k,5)) - ...
          0.94285714285714356*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^6*(-0.34375*(f(:,k,1)-f(:,k,7)) + ...
          0.68809921051219412*(f(:,k,2)-f(:,k,6)) - ...
          0.48528739061765716*(f(:,k,3)-f(:,k,5)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^5*(-0.375*(f(:,k,1)+f(:,k,7)) + ...
          1.21320038050366995*(f(:,k,2)+f(:,k,6)) - ...
          2.09820038050367082*(f(:,k,3)+f(:,k,5)) + ...
          2.52000000000000176*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(0.46875*(f(:,k,1)-f(:,k,7)) - ...
          1.25903493358549612*(f(:,k,2)-f(:,k,6)) + ...
          1.22967339607367386*(f(:,k,3)-f(:,k,5)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(0.10416666666666660*(f(:,k,1)+f(:,k,7)) - ...
          0.36437739881046465*(f(:,k,2)+f(:,k,6)) + ...
          1.42687739881046537*(f(:,k,3)+f(:,k,5)) - ...
          2.3333333333333346*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(-0.15625*(f(:,k,1)-f(:,k,7)) + ...
          0.45377223563440987*(f(:,k,2)-f(:,k,6)) - ...
          1.00348462029437623*(f(:,k,3)-f(:,k,5)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*(1.0*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          0.55059523809523700e-1*f(:,k,1) + ...
          0.25557651111967517*f(:,k,2) + ...
          0.47497130544329094*f(:,k,3) + ...
          0.24380952380952357*f(:,k,4) - ...
          0.4322592423342813e-1*f(:,k,5) + ...
          0.2124953624189091e-1*f(:,k,6) - ...
          0.7440476190476161e-2*f(:,k,7);
    end
  end
  % order 6 Gauss-Lobatto

elseif (o.orderGL == -7)
  t = [-1 -0.6666666666666666  ...
      -0.3333333333333333 0 0.3333333333333333 ...
      0.6666666666666667 1];

  for n = 1:7
    for k = 1:nv
      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^7*(0.1446428571428571*(f(:,k,1)+f(:,k,7)) - ...
          0.8678571428571429*(f(:,k,2)+f(:,k,6)) + ...
          2.1696428571428571*(f(:,k,3)+f(:,k,5)) - ...
          2.8928571428571429*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^6*(-0.16875*(f(:,k,1)-f(:,k,7)) + ...
          0.675*(f(:,k,2)-f(:,k,6)) - ...
          0.84375*(f(:,k,3)-f(:,k,5)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^5*(-0.1125*(f(:,k,1)+f(:,k,7)) + ...
          1.35*(f(:,k,2)+f(:,k,6)) - ...
          4.3875*(f(:,k,3)+f(:,k,5)) + ...
          6.3*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^4*(0.140625*(f(:,k,1)-f(:,k,7)) - ...
          1.125*(f(:,k,2)-f(:,k,6)) + ...
          1.828125*(f(:,k,3)-f(:,k,5)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^3*(0.0166666666666667*(f(:,k,1)+f(:,k,7)) - ...
          0.225*(f(:,k,2)+f(:,k,6)) + ...
          2.25*(f(:,k,3)+f(:,k,5)) - ...
          4.0833333333333333*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n).^2*(-0.025*(f(:,k,1)-f(:,k,7)) + ...
          0.225*(f(:,k,2)-f(:,k,6)) - ...
          1.125*(f(:,k,3)-f(:,k,5)));

      fint(:,k,n) = fint(:,k,n) + ...
          t(n)*(1.0*f(:,k,4));

      fint(:,k,n) = fint(:,k,n) + ...
          0.1019345238095238*f(:,k,1) + ...
          0.4821428571428571*f(:,k,2) + ...
          0.1727678571428571*f(:,k,3) + ...
          0.3238095238095238*f(:,k,4) - ...
          0.1084821428571429*f(:,k,5) + ...
          0.03214285714285714*f(:,k,6) - ...
          0.004315476190476190*f(:,k,7);
    end
  end
  % order 6 equi-spaced points
end % o.orderGL

end % lobattoInt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sdcallrhs
function [rhs,vesicle,Xo] = sdcallrhs(o,VES,WALL,PROP,SDCPROP)
% VES.X, VES.sig, VES.u, WALL.eta, WALL.RS, WALL.walls, PROP.kappa,
% PROP.viscCont, PROP.updatePreco,
% SDCPROP.deltaX, SDCPROP.deltaSig, SDCPROP.deltaEta, SDCPROP.deltaRS
% SDCPROP.diffResidual, SDCPROP.vesVel, SDCPROP.wallVel, SDCPROP.sa,
% SDCPROP.IK.

if o.SDCcorrect
  deltaX = SDCPROP.deltaX;
  deltaSig = SDCPROP.deltaSig;
  deltaEta = SDCPROP.deltaEta;
  deltaRS = SDCPROP.deltaRS;
  diffResidual = SDCPROP.diffResidual;
  vesVel = SDCPROP.vesVel;
  wallVel = SDCPROP.wallVel;
  sa = SDCPROP.sa;
  IK = SDCPROP.IK;
end

Xstore = VES.X;
sigStore = VES.sig;
uStore = VES.u;
walls = WALL.walls;
etaStore = WALL.eta;
RSstore = WALL.RS;
viscCont = PROP.viscCont;
kappa = PROP.kappa;
updatePreco = PROP.updatePreco;

N = size(Xstore,1)/2; % Number of points per vesicle
nv = size(Xstore,2); % Number of vesicles
if o.confined
  Xwalls = walls.X; % discretization points of solid walls
else
  Xwalls = [];
end
Nbd = size(Xwalls,1)/2; % Number of points on the solid walls
nvbd = size(Xwalls,2); % number of solid wall components
alpha = (1 + viscCont)/2; 
% constant that appears in front of time derivative in
% vesicle dynamical equations

Xm = zeros(2*N,nv);
sigmaM = zeros(N,nv);
uM = zeros(2*N,nv);
Xo = zeros(2*N,nv);
etaM = zeros(2*Nbd,nvbd);
RSm = zeros(3,nvbd);
for k = 1:o.order
  Xm = Xm + Xstore(:,:,k)*o.Xcoeff(k);
  sigmaM = sigmaM + sigStore(:,:,k)*o.Xcoeff(k);
  uM = uM + uStore(:,:,k)*o.Xcoeff(k);
  if o.confined
    etaM = etaM + etaStore(:,:,k)*o.Xcoeff(k);
    RSm = RSm + RSstore(:,:,k)*o.Xcoeff(k);
  end
  Xo = Xo + Xstore(:,:,k)*o.rhsCoeff(k);
end
% Form linear combinations of previous time steps needed
% for Ascher, Ruuth, and Wetton IMEX methods

vesicle = capsules(Xm,sigmaM,uM,kappa,viscCont,o.antiAlias);
% build an object vesicle that contains tangent vector, jacobian, etc.

op = o.op;
if ~o.SDCcorrect
  o.Galpert = op.stokesSLmatrix(vesicle);
%  o.Gbarnett = op.laplaceSLcomplexMatrix(vesicle);
end
% Build single layer potential matrix and put it in current object
% If we are doing an sdc update, this is already precomputed and 
% stored from when we formed the provisional solution

if ~o.SDCcorrect
  if any(viscCont ~= 1)
    o.D = op.stokesDLmatrix(vesicle);
  else
    o.D = [];
  end
end
% Compute double-layer potential matrix due to each vesicle
% independent of the others.  Matrix is zero if there is no
% viscosity contrast

if ~o.SDCcorrect
  if o.near
    if o.confined
      [o.NearV2V,o.NearV2W] = vesicle.getZone(walls,3);
      % Need vesicle to vesicle and vesicle to wall interactions
      if nvbd == 1 
        [o.NearW2W,o.NearW2V] = walls.getZone(vesicle,2);
      else
        [o.NearW2W,o.NearW2V] = walls.getZone(vesicle,3);
      end
      % Only need wall to vesicle interactions.  Wall to wall
      % interactions should also use near-singular integration since
      % they may be close to one another
    else
      o.NearV2V = vesicle.getZone([],1);
      % no solid walls, so only need vesicle-vesicle intearactions
      o.NearV2W = [];
      o.NearW2V = [];
      o.NearW2W = [];
    end
  else
    o.NearV2V = [];
    o.NearV2W = [];
    o.NearW2V = [];
    o.NearW2W = [];
  end
  % If using near-singular integration, need structures for deciding who
  % is close, how close it is, who is closest, etc.
end
% Only form near-singular integration structure if not doing an SDC
% update.  Otherwise this was formed and saved when forming the
% provisional solution

if o.collDete
  o.lapDLP = op.laplaceDLmatrix(vesicle);
  % build matrix that evaluates double-layer laplace potential for each
  % vesicle

  oc = curve;
  [icollisionVes,icollisionWall] = ...
    oc.collision(vesicle,walls,o.NearV2V,o.NearV2W,o.fmm,o.near);
  % Check for collisions 
  if icollisionVes
    fprintf('VESICLES HAVE CROSSED\n')
    pause
  end
  if icollisionWall
    fprintf('VESICLES HAVE CROSSED SOLID WALL\n')
    pause
  end
end
% check for collisions

if ~o.SDCcorrect
  rhs1 = Xo;
  rhs2 = zeros(N,nv);
  if o.confined
    rhs3 = walls.u;
  else
    rhs3 = [];
  end
else
  if strcmp(o.vesves,'explicit')
    rhs1 = diffResidual;
  else
    rhs1 = deltaX + diffResidual;
  end
  if any(vesicle.viscCont ~= 1)
    z = zeros(2*N*nv,1);
    for k = 1:nv
      z(2*(k-1)*N+1:2*k*N) = rhs1(:,k);
    end
    if strcmp(o.vesves,'explicit')
      z = o.IminusDexp(z,vesicle);
    else
      z = o.IminusD(z,vesicle);
    end
    for k = 1:nv
      rhs1(:,k) = z(2*(k-1)*N+1:2*k*N)/alpha(k);
    end
  end
  if strcmp(o.vesves,'explicit')
    rhs1 = rhs1 + deltaX;
  end
  if strcmp(o.solver,'method1')
    rhs2 = ones(N,nv);
  else
    rhs2 = -vesicle.surfaceDiv(vesVel);
  end
  if o.confined
    rhs3 = -wallVel + walls.u;
  else
    rhs3 = [];
  end
end
% Parts of rhs from previous solution.  The right-hand-side depends on
% whether we are doing an SDC correction or forming the provisional
% solution.

% START TO COMPUTE RIGHT-HAND SIDE DUE TO VESICLE TRACTION JUMP
if o.gravity
  % For now, treat gravity explicit same as backgroud flow
  fgrav = o.gravityForce(Xm,vesicle,o.gCont);
else
  fgrav = zeros(2*N,nv);
end
if strcmp(o.vesves,'explicit')
  if o.profile
    tic
  end
  if ~o.SDCcorrect
    f = vesicle.tracJump(Xm,sigmaM) + PROP.fc_tot + fgrav;
  else
    f = vesicle.tracJump(deltaX,deltaSig);
  end
  if ~o.fmm
    kernel = @op.exactStokesSL;
  else
    kernel = @op.exactStokesSLfmm;
  end

  if ~o.near
    Fslp = kernel(vesicle,f,[]);
    % Evaulate single-layer potential on all vesicles but itself without
    % near-singular integration
    if o.confined
      [~,FSLPwall] = kernel(vesicle,f,[],walls.X,(1:nv));
      % Evaluate single-layer potential on solid walls due to all
      % vesicles
    else
      FSLPwall = [];
    end

  else
    if o.nearStrat == 'cauchy'
      Fslp = op.nearSingStokesSLP(vesicle,f,vesicle,...
          o.Gbarnett,true,o.fmm);
    else
      SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);
      SLPtrap = SLP;
      kernelDirect = kernel;
      Fslp = op.nearSingInt(vesicle,f,SLP,SLPtrap,...
          o.NearV2V,kernel,kernelDirect,vesicle,true);
    end
    % Use near-singular integration to compute single-layer potential
    % due to all other vesicles.  Need to pass function op.exactStokesSL
    % so that the layer potential can be computed at far points and
    % Lagrange interpolation points

    if o.confined
      if o.nearStrat == 'cauchy'
        FSLPwall = op.nearSingStokesSLP(vesicle,f,walls,...
            o.Gbarnett,false,o.fmm);
      else
        FSLPwall = op.nearSingInt(vesicle,f,SLP,SLPtrap,...
            o.NearV2W,kernel,kernelDirect,walls,false);
        % Evaluate the velocity on the walls due to the vesicles
      end
    else
      FSLPwall = [];
    end

  end
  if o.profile
    fprintf('Build right-hand side V2V and V2W   %5.1e\n',toc);
  end
end
for k = 1:nv
  rhs1(:,k) = rhs1(:,k) + o.dt/alpha(k)*Fslp(:,k);
end

rhs3 = rhs3 - FSLPwall;
if (o.gravity && ~o.SDCcorrect)
Gfg = op.exactStokesSLdiag(vesicle,o.Galpert,fgrav);
rhs1 = rhs1 + o.dt*Gfg*diag(1./alpha);
end
% END TO COMPUTE RIGHT-HAND SIDE DUE TO VESICLE TRACTION JUMP

% START TO COMPUTE RIGHT-HAND SIDE DUE TO VISCOSITY CONTRAST
if any(vesicle.viscCont ~= 1)
  if o.near
    jump = 1/2*(1-vesicle.viscCont);
    DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle,o.D,X);
  end
  % Need to add jump to double-layer potential if using near-singular
  % integration so that we can compute boundary values for the
  % near-singular integration algorithm

  if ~o.fmmDLP
    kernel = @op.exactStokesDL;
  else
    kernel = @op.exactStokesDLfmm;
  end
  if strcmp(o.vesves,'explicit')
    if ~o.near
      Fdlp = -o.dt * kernel(vesicle,uM);
      % Evaulate the velocity due to the viscosity contrast
      % on all vesicles but itself WITHOUT near-singular
      % integration
      if o.confined
        [~,FDLPwall] = kernel(vesicle,uM,walls.X,(1:nv));
        FDLPwall = -o.dt*FDLPwall;
        % Evaulate the velocity due to the viscosity contrast on the
        % walls WITHOUT near-singulation integration
      else
        FDLPwall = [];
      end
    else
      kernelDirect = kernel;
      Fdlp = -o.dt * op.nearSingInt(vesicle,uM,DLP,DLP,...
          o.NearV2V,kernel,kernelDirect,vesicle,true);
      % Evaulate the velocity due to the viscosity contrast on all
      % vesicles but itself WITH near-singular integration

      if o.confined
        FDLPwall = -o.dt * op.nearSingInt(vesicle,uM,DLP,DLP,...
            o.NearV2W,kernel,kernelDirect,walls,false);
        % Evaulate the velocity due to the viscosity contrast on the
        % walls WITH near-singulation integration
      else
        FDLPwall = [];
      end
    end
  end
else
  Fdlp = zeros(2*N,nv);
  FDLPwall = zeros(2*Nbd,nvbd);
  % If no viscosity contrast, there is no velocity induced due to a
  % viscosity contrast
end

if (any(viscCont ~= 1) && o.SDCcorrect && strcmp(o.vesves,'explicit'))
  DXo = op.exactStokesDLdiag(vesicle,o.D,deltaX);
  rhs1 = rhs1 - (Fdlp + DXo) * diag(1./alpha);
end
if (any(viscCont ~= 1) && ~o.SDCcorrect)
  DXo = op.exactStokesDLdiag(vesicle,o.D,Xo);
  rhs1 = rhs1 - (Fdlp + DXo) * diag(1./alpha);
end
% add in viscosity contrast term due to each vesicle independent of the
% others (o.D * Xo) from the previous solution followed by the term due
% to all other vesicles (Fdlp)

rhs3 = rhs3 + FDLPwall/o.dt;
% compute the double-layer potential due to all other vesicles from the
% appropriate linear combination of previous time steps.  Depends on
% time stepping order and vesicle-vesicle discretization

if o.adhesion
  adhesion = vesicle.adhesionTerm(o.adStrength,o.adRange);
  if ~o.fmm
    kernel = @op.exactStokesSL;
  else
    kernel = @op.exactStokesSLfmm;
  end

  Fadhesion = op.exactStokesSLdiag(vesicle,o.Galpert,adhesion);
  % diagonal term of adhesion

  if ~o.near
    Fadhesion = Fadhesion + kernel(vesicle,adhesion);
    % Evaulate single-layer potential on all vesicles but itself without
    % near-singular integration
  else
    SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);
    SLPtrap = SLP;
    kernelDirect = kernel;
    if o.nearStrat == 'cauchy'
      Fadhesion = Fadhesion + ...
          op.nearSingStokesSLP(vesicle,adhesion,vesicle,...
          o.Gbarnett,true,o.fmm);
    else
      Fadhesion = Fadhesion + ...
          op.nearSingInt(vesicle,adhesion,SLP,SLPtrap,...
          o.NearV2V,kernel,kernelDirect,vesicle,true);
    end
    % Use near-singular integration to compute single-layer potential
    % due to all other vesicles.  Need to pass function op.exactStokesSL
    % so that the layer potential can be computed at far points and
    % Lagrange interpolation points
  end

  rhs1 = rhs1 + o.dt*Fadhesion*diag(1./alpha);
end

% START TO COMPUTE RIGHT-HAND SIDE DUE TO SOLID WALLS
if o.confined
  if strcmp(o.vesves,'explicit')
    if o.profile
      tic
    end
    if o.near
      jump = -1/2;
      DLP = @(X) jump*X + op.exactStokesDLdiag(walls,o.wallDLP,X);
    end
    % compute the matrix for doing evaluating the double-layer potential
    % on the solid walls required for near-singular integration

    if ~o.fmmDLP
      kernel = @op.exactStokesDL;
    else
      kernel = @op.exactStokesDLfmm;
    end
    if ~o.near
      if ~o.SDCcorrect
        charge = etaM;
      else
        charge = deltaEta;
      end
      [~,Fwall2Ves] = kernel(walls,charge,...
          vesicle.X,1:nvbd);
      % velocity field due to the walls evaluated on the vesicle
    else
      if ~o.SDCcorrect
        charge = etaM;
      else
        charge = deltaEta;
      end
      kernelDirect = kernel;
      Fwall2Ves = op.nearSingInt(walls,charge,DLP,DLP,...
          o.NearW2V,kernel,kernelDirect,vesicle,false);
    end

    for k = 2:nvbd
      if ~o.SDCcorrect
        stokeslet = RSm(1:2,k);
        rotlet = RSm(3,k);
      else
        stokeslet = deltaRS(1:2,k);
        rotlet = deltaRS(3,k);
      end
      Fwall2Ves = Fwall2Ves + ...
        o.RSlets(vesicle.X,walls.center(:,k),stokeslet,rotlet);
    end
    if o.profile
      fprintf('Build right-hand side W2V           %5.1e\n',toc);
    end

  end
else
  Fwall2Ves = zeros(2*N,nv);
  if ~o.SDCcorrect
    rhs1 = rhs1 + o.dt*o.farField(Xm)*diag(1./alpha);
    % Add in far-field condition (extensional, shear, etc.)
  end
end
rhs1 = rhs1 + o.dt*Fwall2Ves*diag(1./alpha);
% right-hand side of the velocity evaluated on the solid walls
% END TO COMPUTE RIGHT-HAND SIDE DUE TO SOLID WALLS


% START TO COMPUTE THE RIGHT-HAND SIDE FOR THE INEXTENSIBILITY CONDITION
if (strcmp(o.solver,'method1'))
  % If using method1, vesicle-vesicle interactions and the presence of a
  % viscosity contrast is irrelevent
  if ~o.SDCcorrect
    rhs2 = rhs2 + vesicle.surfaceDiv(Xo); 
  else
    divf = zeros(N,nv);
    for k = 1:nv
      divf(:,k) = curve.arcDeriv(Xm(1:N,k),1,1./sa(:,k),...
          IK(:,k)).^2 + ...
                  curve.arcDeriv(Xm(N+1:2*N,k),1,1./sa(:,k),...
          IK(:,k)).^2; 
    end
    rhs2 = 1/2*(rhs2 - divf);
  end
else 
  % If using method2, method3, or method4, vesicle-vesicle interaction
  % affects the right-hand side of the inextensibility condition

  if ~o.confined && ~o.SDCcorrect
    if any(viscCont ~= 1)
      if o.profile
        tic
      end
      rhs2 = rhs2 - vesicle.surfaceDiv(...
          o.solveIminusD(o.farField(Xm),vesicle));
      if o.profile
        fprintf('Solve system alpha*I - DLP          %5.1e\n',toc);
      end
    else
      rhs2 = rhs2 - vesicle.surfaceDiv(o.farField(Xm));
    end
  end
  % add in term from farfield

  if (strcmp(o.vesves,'explicit'))
    rhs2 = rhs2 - vesicle.surfaceDiv(Fslp); 
    rhs2 = rhs2 - vesicle.surfaceDiv(Fwall2Ves);
    if o.adhesion
      rhs2 = rhs2 + vesicle.surfaceDiv(Fadhesion);
    end
  end
end
% rhs2 is the right-hand side for the inextensibility condition
% END TO COMPUTE THE RIGHT-HAND SIDE FOR THE INEXTENSIBILITY CONDITION

if (any(vesicle.viscCont ~= 1) && ...
      strcmp(o.vesves,'implicit') && o.confined)
  rhs3 = rhs3 * o.dt;
end
% This makes sure that the rhs is all order one rather than have rhs3
% being order 1/o.dt and other two parts (rhs1 and rhs2) being order 1.
% This of course needs to be compensated in the TimeMatVec routine

% START TO COMPUTE RIGHT-HAND SIDE DUE TO NEW BENDING TERM
if o.bending
  f = vesicle.newBending;
  % New part of traction jump if we have variable bending
  rhs1 = rhs + o.dt * op.exactStokesSLdiag(vesicle,o.Galpert,f);
end
% This term is always treated fully explicitly.  The highest derivative
% it contains is second-order and it appears through a curvature
% END TO COMPUTE RIGHT-HAND SIDE DUE TO NEW BENDING TERM.  THIS WILL
% ONLY SUPPORT A SINGLE VESICLE IN AN UNBOUNDED FLOW


rhs = [rhs1; rhs2];
rhs = rhs(:);
rhs = [rhs; rhs3(:)];
% Stack the right-hand sides in an alternating with respect to the
% vesicle fashion
% Add on the no-slip boundary conditions on the solid walls
rhs = [rhs; zeros(3*(nvbd-1),1)];
% Rotlet and Stokeslet equations

usePreco = true;
%usePreco = false;

% START BUILDING BLOCK-DIAGONAL PRECONDITIONER
if usePreco
  if updatePreco
    % only build the preconditioner if updatePreco == true
    if o.profile
      tic
    end
    [Ben,Ten,Div] = vesicle.computeDerivs;
    if o.profile
      fprintf('Build differential operators        %5.1e\n',toc);
    end
    % Compute bending, tension, and surface divergence of current
    % vesicle configuration

    if o.profile
      tic
    end

    if (strcmp(o.solver,'method1') || strcmp(o.solver,'method2'))
      bdiagVes.L = zeros(3*N,3*N,nv);
      bdiagVes.U = zeros(3*N,3*N,nv);
    elseif strcmp(o.solver,'method3')
      bdiagVes.GT = zeros(2*N,N,nv); % SLP of tension
      bdiagVes.DGT = zeros(N,N,nv); % divergence of SLP of tension
      bdiagVes.DGB = zeros(N,2*N,nv); % divergence of SLP of bending
      bdiagVes.schur = zeros(2*N,2*N,nv); 
      % schur complement of the lower right N by N block
    elseif strcmp(o.solver,'method4')
      bdiagVes.GT = zeros(2*N,N,nv); % SLP of tension
      bdiagVes.IpBen = zeros(2*N,2*N,nv); 
      % inverse of identity plus SLP of Bending
      bdiagVes.DGB = zeros(N,2*N,nv); % divergence of SLP of bending
      bdiagVes.schur = zeros(N,N,nv); 
      % schur complement of the upper left 2*N by 2*N block
    end

    for k=1:nv
      if strcmp(o.solver,'method1')
        if any(vesicle.viscCont ~= 1)
          [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu(...
            [o.beta*(eye(2*N) - o.D(:,:,k)/alpha(k)) + ...
                o.dt/alpha(k)*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
            -o.dt/alpha(k)*o.Galpert(:,:,k)*Ten(:,:,k); ...
            o.beta*Div(:,:,k) zeros(N)]);
        else
          [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu(...
            [o.beta*eye(2*N) + ...
                o.dt*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
            -o.dt/alpha(k)*o.Galpert(:,:,k)*Ten(:,:,k); ...
            o.beta*Div(:,:,k) zeros(N)]);
        end
      elseif strcmp(o.solver,'method2')
        if any(vesicle.viscCont ~= 1)
          [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu( ...
            [o.beta*(eye(2*N) - o.D(:,:,k)/alpha(k)) + ...
                o.dt/alpha(k)*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
            -o.dt/alpha(k)*o.Galpert(:,:,k)*Ten(:,:,k); ...
            -vesicle.kappa*Div(:,:,k)*(alpha(k)*eye(2*N) - o.D(:,:,k))\...
                (o.Galpert(:,:,k)*Ben(:,:,k)) ...
            Div(:,:,k)*(alpha(k)*eye(2*N) - o.D(:,:,k))\...
                (o.Galpert(:,:,k)*Ten(:,:,k))]);
        else
          [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu(...
            [o.beta*eye(2*N) + ...
                o.dt*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
            -o.dt*o.Galpert(:,:,k)*Ten(:,:,k); ...
            -vesicle.kappa*Div(:,:,k)*o.Galpert(:,:,k)*Ben(:,:,k) ...
            Div(:,:,k)*o.Galpert(:,:,k)*Ten(:,:,k)]);
        end
      elseif strcmp(o.solver,'method3')
        % schur complement of lower right block of method2
        bdiagVes.GT(:,:,k) = o.Galpert(:,:,k)*Ten(:,:,k);
        bdiagVes.DGT(:,:,k) = (Div(:,:,k)*o.Galpert(:,:,k)*...
            Ten(:,:,k))\eye(N);
        bdiagVes.DGB(:,:,k) = ...
          vesicle.kappa*Div(:,:,k)*o.Galpert(:,:,k)*Ben(:,:,k);
        bdiagVes.schur(:,:,k) = ...
          inv((o.beta*eye(2*N) + vesicle.kappa*o.dt*...
          o.Galpert(:,:,k)*Ben(:,:,k)) - ...
          o.dt*o.Galpert(:,:,k)*Ten(:,:,k)*bdiagVes.DGT(:,:,k)*...
          vesicle.kappa*Div(:,:,k)*o.Galpert(:,:,k)*Ben(:,:,k));

      elseif strcmp(o.solver,'method4')
        % schur complement of upper left block of method2
        bdiagVes.GT(:,:,k) = o.Galpert(:,:,k)*Ten(:,:,k);
        bdiagVes.IpBen(:,:,k) = inv(o.beta*eye(2*N) + ...
            vesicle.kappa*o.dt*o.Galpert(:,:,k)*Ben(:,:,k));
        bdiagVes.DGB(:,:,k) = ...
          vesicle.kappa*Div(:,:,k)*o.Galpert(:,:,k)*Ben(:,:,k);
        bdiagVes.schur(:,:,k) = ...
          inv(Div(:,:,k)*(-vesicle.kappa*o.dt*o.Galpert(:,:,k)*...
          Ben(:,:,k)*bdiagVes.IpBen(:,:,k) + eye(2*N))*...
          bdiagVes.GT(:,:,k));
      end
    end
    o.bdiagVes = bdiagVes;
    % Build block-diagonal preconditioner of self-vesicle 
    % intearctions in matrix form
    if o.profile
      fprintf('Build block-diagonal preconditioner %5.1e\n\n',toc);
    end
  end % updatePreco
end % usePreco

end % sdcallrhs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GSSDCCOL
function [X,sigma,u,eta,RS,iter,accept,dtScale,normRes,iflag,FC,Xrig,etaRig] = ...
    GSSDCCOL(o,Xstore,sigStore,uStore,...
    etaStore,RSstore,kappa,viscCont,walls,wallsCoarse,om,time,accept,FC,XrigStore,etaRigStore)
% only works with explicit scheme and explicit with collision
timeCol = 0;
N = size(Xstore,1)/2; % number of points per vesicle
nv = size(Xstore,2); % number of vesicles
Nbd = size(etaStore,1)/2; % number of points per solid wall
nvbd = size(etaStore,2); % number of solid walls
Nrd = size(etaRigStore,1)/2;
nvrd = size(etaRigStore,2);
a0 = om.area; % input area
l0 = om.length; % input length
oc = curve;

Xprov = zeros(2*N,nv,abs(o.orderGL));
sigmaProv = zeros(N,nv,abs(o.orderGL));
uProv = zeros(2*N,nv,abs(o.orderGL));
etaProv = zeros(2*Nbd,nvbd,abs(o.orderGL));
RSprov = zeros(3,nvbd,abs(o.orderGL));
FCprov = zeros(2*N,nv,abs(o.orderGL));
XrigProv = zeros(2*Nrd,nvrd,abs(o.orderGL));
etaRigProv = zeros(2*Nrd,nvrd,abs(o.orderGL));
% need to save the position, tension, velocity, density
% functions, rotlets/stokeslets at the Gauss-Lobatto
% time steps
for k = 1:nv
  Xprov(:,k,1) = Xstore(:,k);
  sigmaProv(:,k,1) = sigStore(:,k);
  uProv(:,k,1) = uStore(:,k);
end
% Store the initial conditions in the provisional
% solution which will store the solution at all
% Gauss-Lobatto points

for k = 1:nvrd
  XrigProv(:,k,1) = XrigStore(:,k);
  etaRigProv(:,k,1) = etaRigStore(:,k);
end

for k = 1:nvbd
  etaProv(:,k,1) = etaStore(:,k);
  RSprov(:,k,1) = RSstore(:,k);
end

Galpert = zeros(2*N,2*N,nv,abs(o.orderGL));
%Gbarnett = zeros(N,N,nv,abs(o.orderGL));
Gbarnett = [];
% need the single-layer potential at all levels
% of the provisional solution
if any(viscCont ~= 1)
  D = zeros(2*N,2*N,nv,abs(o.orderGL));
else
  D = [];
end
if o.withRigid
  rigidDLP = zeros(2*Nrd,2*Nrd,nvrd,abs(o.orderGL));
else
  rigidDLP = [];
end
% need the double-layer potential at all levels
% of the provisional solution

deltaX = zeros(2*N,nv,abs(o.orderGL));
deltaSigma = zeros(N,nv,abs(o.orderGL));
deltaEta = zeros(2*Nbd,nvbd,abs(o.orderGL));
deltaRS = zeros(3,nvbd,abs(o.orderGL));
deltaFC = zeros(2*N,nv,abs(o.orderGL));
% Errors that are solved for in each iteration.  The first column of
% deltaX will always be zero since we assume the inital condition is
% exact

X = Xprov(:,:,1);
sigma = sigmaProv(:,:,1);
u = uProv(:,:,1);

vesicle(1) = capsules(X,sigma,u,kappa,viscCont,o.antiAlias);
% vesicle configuration at first Gauss-Lobatto point

Xrig = XrigProv(:,:,1);
etaRig = etaRigProv(:,:,1);
rigid(1) = capsules(Xrig,[],[],zeros(nvrd,1),zeros(nvrd,1),o.antiAlias);

if o.timeAdap
  [~,aInit,lInit] = oc.geomProp(X);
  eaInit = abs(aInit - a0)./a0;
  elInit = abs(lInit - l0)./l0;
end
% compute intial reduced area, area, and length if using adaptive
% time stepping

dt = o.dt;
% need to save the time step size
dtimes = diff(o.GLpts)*dt/2;
% time step sizes use for Gauss-Lobatto points

% START OF FORMING PROVISIONAL SOLUTION
iflag = 0; % initialize gmres flag as everything okay
iter = 0; % initialize number of gmres iterations
updatePreco = true; 
fc_tot = FC;
% Form preconditioner at initial configuartion
for k = 1:abs(o.orderGL)-1
  o.dt = dtimes(k);
  o.SDCcorrect = false;

  VES.X = Xprov(:,:,k-o.order+1:k); 
  VES.sig = sigmaProv(:,:,k-o.order+1:k);
  VES.u = uProv(:,:,k-o.order+1:k);

  WALL.walls = walls;
  WALL.eta = etaProv(:,:,k-o.order+1:k);
  WALL.RS = RSprov(:,:,k-o.order+1:k); 

  RIGID.rigid = rigid(k);
  RIGID.etaRig = etaRigProv(:,:,k-o.order+1:k);

  PROP.viscCont = viscCont;
  PROP.kappa = kappa;
  PROP.updatePreco = updatePreco;
  PROP.fc_tot = fc_tot;
  PROP.om = om;

  SDCPROP = [];

  if nv 
    ticsdcvesrhs = tic;
    [rhs, ves, Xo] = o.sdcvesrhs(VES,WALL,PROP,SDCPROP,RIGID); 
   % [rhs, ves, Xo,testen] = o.sdcvesrhs(VES,WALL,PROP,SDCPROP,RIGID); 
    tocsdcvesrhs = toc(ticsdcvesrhs);
    om.writeMessage(['time for sdcvesrhs: ' num2str(tocsdcvesrhs)],'%s\n');
    ticvesmat = tic;
    [Xn,iflagTemp,R,I,resvec] = gmres(@(X) o.vesmat(X,ves),...
          rhs,[],o.gmresTol,o.gmresMaxIter,...
          @o.preconditionerBDVes);
    tocvesmat = toc(ticvesmat);
    om.writeMessage(['time for vesmat: ' num2str(tocvesmat)],'%s\n');
    subIter = I(2);
    iter = iter + subIter;
    [X,sigma,u] = o.extractRHSVes(Xn,N,nv,Xo);
    if iflagTemp ~= 0
      iflag = iflagTemp;
    end
    % if any of the gmres iterations fail, assign the failure to
    % flag to iflag for monitor to output
%    Idx = ifftshift(-N/2:N/2+1);
%    Idx = abs(Idx) >= N/4;
%    coefften = fft(sigma);
%    sigmaHighFreqNorm = norm(coefften(Idx))
%    coeffx = fft(testen(1:N));
%    coeffy = fft(testen(N+1:2*N));
%    rhsHighFreqNorm = sqrt(norm(coeffx(Idx))^2+norm(coeffy(Idx))^2)
%
%    [Xntmp,iflagTemp,R,I,resvec] = gmres(@(X) o.vesmat(X,ves),...
%          [testen;zeros(N,1)],[],o.gmresTol*1e-2,o.gmresMaxIter,...
%          @o.preconditionerBDVes);
%    subIter = I(2);
%    iter = iter + subIter;
%    [Xten,sigmaten,uten] = o.extractRHSVes(Xntmp,N,nv,Xo);
%    coefften = fft(sigmaten);
%    sigmaFcHighFreqNorm = norm(coefften(Idx))
  else
     ves = [];
     %op = poten(Nrd);
     op = o.op;
     o.rigidDLP = op.stokesDLmatrix(rigid);
     o.NearR2R = rigid.getZone([],1);
     if o.confined
       [~,o.NearR2W] = rigid.getZone(walls,2);
       [~,o.NearW2R] = walls.getZone(rigid,2);
     end
     if o.withRigid
       bdiagRigid.L = zeros(2*Nrd+3,2*Nrd+3,nvrd);
       bdiagRigid.U = zeros(2*Nrd+3,2*Nrd+3,nvrd);
       for kk=1:nvrd
         cmR = rigid.center(:,kk);
         Xr = rigid.X(:,kk);
         saR = rigid.sa(:,kk);
         [bdiagRigid.L(:,:,kk),bdiagRigid.U(:,:,kk)] = ... 
	 lu([0.5*eye(2*Nrd)-o.rigidDLP(:,:,kk) ... 
	 [ones(Nrd,1);zeros(Nrd,1)] [zeros(Nrd,1);ones(Nrd,1)] ... 
	 [Xr(Nrd+1:end)-cmR(2);-(Xr(1:Nrd)-cmR(1))];...
         [saR'*2*pi/Nrd zeros(1,Nrd) 0 0 0];[zeros(1,Nrd) saR'*2*pi/Nrd 0 0 0];...
	 [2*pi/Nrd*(Xr(Nrd+1:end)'-cmR(2)).*saR' -2*pi/Nrd*(Xr(1:Nrd)'-cmR(1)).*saR' 0 0 0]]);
       end
       o.bdiagRigid = bdiagRigid;
     end
  end

  if o.withRigid
    VES.X = X; 
    VES.sig = sigma;
    VES.u = u;
    [rhsRig] = o.sdcrigidrhs(VES,WALL,PROP,SDCPROP,ves,RIGID);
    [Xn,iflagTemp,R,I,resvec] = gmres(@(X) o.rigidmat(X,rigid(k)),rhsRig,[],o.gmresTol,(2*Nrd+3)*nvrd,@o.preconditionerBDRigid);
    %[Xn,iflagTemp,R,I,resvec] = gmres(@(X) o.rigidmat(X,rigid(k)),rhsRig,[],o.gmresTol,(2*Nrd+3)*nvrd);
    [Xrig,etaRig,uRig,omgRig] = o.extractRHSrigid(Xn,rigid(k));
    if iflagTemp ~= 0
      iflag = iflagTemp;
    end
%    myx1 = o.bdiagRigid.U(:,:,1)\(o.bdiagRigid.L(:,:,1)\rhsRig(1:2*Nrd+3));
%    myx2 = o.bdiagRigid.U(:,:,2)\(o.bdiagRigid.L(:,:,2)\rhsRig(1+2*Nrd+3:end));
%    norm([myx1;myx2]-Xn)
    % if any of the gmres iterations fail, assign the failure to
    % flag to iflag for monitor to output
  end

  if o.confined
    VES.X = X; 
    VES.sig = sigma;
    VES.u = u;
    RIGID.etaRig = etaRig;
    ticsdcwallrhs = tic;
    [rhs] = o.sdcwallrhs(VES,WALL,PROP,SDCPROP,ves,RIGID);
    tocsdcwallrhs = toc(ticsdcwallrhs);
    om.writeMessage(['time for sdcwallrhs: ' num2str(tocsdcwallrhs)],'%s\n');
    ticwallmat = tic;
    [Xn,iflagTemp,R,I,resvec] = gmres(@(X) o.wallmat(X,walls),...
          rhs,[],o.gmresTol,o.gmresMaxIter,...
          @o.preconditionerBDWall);
    tocwallmat = toc(ticwallmat);
    om.writeMessage(['time for wallmat: ' num2str(tocwallmat)],'%s\n');
    [eta,RS] = o.extractRHSWall(Xn,Nbd,nvbd);
    if iflagTemp ~= 0
      iflag = iflagTemp;
    end
    % if any of the gmres iterations fail, assign the failure to
    % flag to iflag for monitor to output
  else
    eta = zeros(2*Nbd,nvbd);
    RS = zeros(3,nvbd);
  end

  if(o.resolveCol)
    %updatePreco = false;
    updatePreco = true;
  else
    updatePreco = false;
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% add collision resolve
  if(o.resolveCol)
    ticcol = tic;
    PROPCOL.viscCont = viscCont;
    PROPCOL.vesicle = ves;
    PROPCOL.N = N;
    PROPCOL.nv = nv;
    PROPCOL.Nbd = Nbd;
    PROPCOL.nvbd = nvbd;
    PROPCOL.Nrd = Nrd;
    PROPCOL.nvrd = nvrd;
    PROPCOL.X = X;
    PROPCOL.sigma = sigma;
    PROPCOL.u = u;
    PROPCOL.eta = eta;
    PROPCOL.RS = RS;
    PROPCOL.Xrig = Xrig;
    PROPCOL.etaRig = etaRig;
    PROPCOL.om = om;
    PROPCOL.time = time;
    if o.withRigid
      PROPCOL.uRig = uRig;
      PROPCOL.omgRig = omgRig;
      PROPCOL.rhsRig = rhsRig; 
    end
    [X,sigma,u,eta,RS,Xrig,etaRig,fc_tot] = o.resolveCollisionGSrigid([Xprov(:,:,k) XrigProv(:,:,k)],[X Xrig],walls,rigid(k),PROPCOL);
    toccol = toc(ticcol);
    timeCol = timeCol + toccol;
    om.writeMessage('After resolve collision.','%s\n');
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of add collision resolve
  if o.nsdc > 0
    Galpert(:,:,:,k) = o.Galpert;
%    Gbarnett(:,:,:,k) = o.Gbarnett;
    if any(viscCont ~= 1)
      D(:,:,:,k) = o.D;
    end
    if o.withRigid
      rigidDLP(:,:,:,k) = o.rigidDLP;
    end
    % want to save single-layer potential for computing the residual.
    % This will not save the last one, but this is computed in
    % computeResidual
    if nv
      NearV2V{k} = o.NearV2V;
      NearW2V{k} = o.NearW2V;
      NearV2W{k} = o.NearV2W;
      NearW2W{k} = o.NearW2W;
    end
    if o.withRigid
      NearV2R{k} = o.NearV2R;
      NearW2R{k} = o.NearW2R;
      NearR2W{k} = o.NearR2W;
      NearR2R{k} = o.NearR2R;
      NearR2V{k} = o.NearR2V;
    end
    vesicle(k+1) = capsules(X,sigma,u,kappa,viscCont,o.antiAlias);
    rigid(k+1) = capsules(Xrig,[],[],zeros(nvrd,1),zeros(nvrd,1),o.antiAlias);
    % need to save if doing SDC corrections 
  end

  Xprov(:,:,k+1) = X;
  sigmaProv(:,:,k+1) = sigma;
  uProv(:,:,k+1) = u;
  etaProv(:,:,k+1) = eta;
  RSprov(:,:,k+1) = RS;
  FCprov(:,:,k+1) = fc_tot;
  XrigProv(:,:,k+1) = Xrig;
  etaRigProv(:,:,k+1) = etaRig;

  % save the provisional solution at the intermediate 
  % Gauss-Lobatto point
end % k


% need to save the near-singular integration structure and
% layer-potential matricies at the final Gauss-Lobatto point for future
% SDC sweeps and computing the residual, if there are any SDC sweeps
if o.nsdc > 0
  if o.near
    if o.confined
      [NearV2V{abs(o.orderGL)},NearV2W{abs(o.orderGL)}] = ...
          vesicle(abs(o.orderGL)).getZone(walls,3);
      % Need vesicle to vesicle and vesicle to wall interactions
      if nvbd == 1 
        [NearW2W{abs(o.orderGL)},NearW2V{abs(o.orderGL)}] = ...
            walls.getZone(vesicle(abs(o.orderGL)),2);
      else
        [NearW2W{abs(o.orderGL)},NearW2V{abs(o.orderGL)}] = ...
            walls.getZone(vesicle(abs(o.orderGL)),3);
      end
      % Only need wall to vesicle interactions.  Wall to wall 
      % interactions should also use near-singular integration since
      % they may be close to one another
    else
      NearV2V{abs(o.orderGL)} = ...
          vesicle(abs(o.orderGL)).getZone([],1);
      % no solid walls, so only need vesicle-vesicle intearactions
      NearV2W{abs(o.orderGL)} = [];
      NearW2V{abs(o.orderGL)} = [];
      NearW2W{abs(o.orderGL)} = [];
    end
  else
    NearV2V{abs(o.orderGL)} = [];
    NearV2W{abs(o.orderGL)} = [];
    NearW2V{abs(o.orderGL)} = [];
    NearW2W{abs(o.orderGL)} = [];
  end
  % Need near-singular integration strucutures for final state
  % for sdc updates and computing the residual
  Galpert(:,:,:,abs(o.orderGL)) = ...
      o.op.stokesSLmatrix(vesicle(o.orderGL));
%  Gbarnett(:,:,:,abs(o.orderGL)) = ...
%      o.op.laplaceSLcomplexMatrix(vesicle(o.orderGL));
  if any(viscCont ~= 1)
    D(:,:,:,abs(o.orderGL)) = ...
        o.op.stokesDLmatrix(vesicle(o.orderGL));
  end
  % need single-layer potential at the final state of the 
  % provisional solution to form the residual
end
% END OF FORMING PROVISIONAL SOLUTION AND THE SINGLE-LAYER
% POTENTIAL AT THESE INTERMEDIATE STATES

o.dt = dt;
% change back to original time step
 
%color = ['r' 'g' 'b' 'k' 'c' 'm' 'y'];
%figure(2); clf; hold on
%for k = 1:abs(o.orderGL)
%  plot(sigmaProv(:,1,k),color(k))
%end
%figure(3); clf; hold on
%for k = 1:abs(o.orderGL)
%  plot(squeeze(etaProv(1:end/2,1,k)),color(k))
%end
%figure(5);clf;hold on;
%for k = 1:abs(o.orderGL)
%  plot(etaProv(end/2+1:end,1,k),color(k))
%end
%pause
%% DEBUG: TO MAKE SURE THAT THE TENSION IS CONTINUOUS FROM BETWEEN
%% THE DIFFERENT GAUSS-LOBATTO POINTS

Xinit = Xprov;
uInit = uProv;
sigmaInit = sigmaProv;
etaInit = etaProv;
RSinit = RSprov;
% save the very first provisional solution so that it can be taken
% as the solution if sdc corrections increase the error but this
% original provisional solution has the smallest error
if o.nsdc > 0
  [vesVel,divVesVel,wallVel,residual] = ...
      o.computeVelocitiessdccol(FCprov,vesicle,Galpert,Gbarnett,D,walls,...
        etaProv,RSprov,NearV2V,NearV2W,NearW2V,NearW2W,uProv);
  om.writeMessage('After compute residual.','%s\n');
  % form the residual as well as the velocities due to the different
  % components which are necessary to form the right-hand sides when
  % doing sdc corrections

  normRes = max(abs(residual(:,:,end)));
  normRes = max(normRes);
  % save the size of the residual at the final Gauss-Lobatto
  % point.  We may have more of these later if we do a full
  % deferred correction method.  Use the maximum L2 error where
  % the maximum is taken over all the vesicles
else
  normRes = 0;
end
% compute the integrand that is in the residual and the residual 
% of the picard integral formulation of the dynamic equation
% If the residual is not desired, then errors is set to 0 and
% this will never refine the time step size

[~,a,l] = oc.geomProp(Xprov(:,:,end));
ea = abs(a - a0)./abs(a0);
el = abs(l - l0)./abs(l0);
eaVec = ea;
elVec = el;
message = ['sdcCount   ' num2str(0,'%2d') ...
  ': Residual = ' num2str(normRes(end),'%5.2e') ...
  ', eA = ' num2str(max(ea),'%5.2e') ...
  ', eL = ' num2str(max(el),'%5.2e')];
om.writeMessage(message,'%s\n');
% print the provisional solution's residual and errors

% Start doing SDC corrections
for sdcCount = 1:o.nsdc
uDeltaX = 0*deltaX(:,:,1);
  for n = 1:abs(o.orderGL) - 1
    o.dt = dtimes(n);
    if(o.resolveCol)
      %updatePreco = false;
      updatePreco = true;
    else
      updatePreco = false;
    end
    o.SDCcorrect = true;
    o.Galpert = Galpert(:,:,:,n+1);
%    o.Gbarnett = Gbarnett(:,:,:,n+1);
    if any(viscCont ~= 1)
      o.D = D(:,:,:,n+1);
    end
    o.NearV2V = NearV2V{n+1};
    o.NearV2W = NearV2W{n+1};
    o.NearW2V = NearW2V{n+1};
    o.NearW2W = NearW2W{n+1};
    
    VES.X = Xprov(:,:,n+1); 
    VES.sig = sigmaProv(:,:,n+1);
    VES.u = uDeltaX;

    WALL.walls = walls;
    WALL.eta = etaProv(:,:,n+1);
    WALL.RS = RSprov(:,:,n+1); 
    
    RIGID.rigid = rigid(k);
    RIGID.etaRig = etaRigProv(:,:,k-o.order+1:k);

    PROP.viscCont = viscCont;
    PROP.kappa = kappa;
    PROP.updatePreco = updatePreco;
    PROP.fc_tot = deltaFC(:,:,n);

    SDCPROP.deltaX = deltaX(:,:,n);
    SDCPROP.deltaSig = deltaSigma(:,:,n);
    SDCPROP.deltaEta = deltaEta(:,:,n);
    SDCPROP.deltaRS = deltaRS(:,:,n);
    SDCPROP.diffResidual = residual(:,:,n+1) - residual(:,:,n);
    SDCPROP.vesVel = vesVel(:,:,n+1);
    SDCPROP.wallVel = wallVel(:,:,n+1);
    SDCPROP.sa = vesicle(1).sa;
    SDCPROP.IK = vesicle(1).IK;

    %{
    [rhs, ves, Xo] = o.sdcallrhs(VES,WALL,PROP,SDCPROP); 
    [Xn,iflagTemp,R,I,resvec] = gmres(@(X) o.TimeMatVec(X,ves,walls,0),...
        rhs,[],o.gmresTol,o.gmresMaxIter,...
        @o.preconditionerBD);
    subIter = I(2);
    [X,sigma,u,eta,RS] = o.extractRHS(Xn,N,nv,Nbd,nvbd,Xo);
    iter = iter + subIter;
    if iflagTemp ~= 0
      iflag = iflagTemp;
    end
    %}

    if nv
      [rhs, ves, Xo] = o.sdcvesrhs(VES,WALL,PROP,SDCPROP,RIGID); 
      warning off
      [Xn,iflagTemp,R,I,resvec] = gmres(@(X) o.vesmat(X,ves),...
	    rhs,[],o.gmresTol,o.gmresMaxIter,...
	    @o.preconditionerBDVes);
      subIter = I(2);
      iter = iter + subIter;
      [X,sigma,u] = o.extractRHSVes(Xn,N,nv,Xo);
      if iflagTemp ~= 0
        iflag = iflagTemp;
      end
      % if any of the gmres iterations fail, assign the failure to
      % flag to iflag for monitor to output
    else
       ves = [];
       %op = poten(Nrd);
       op = o.op;
       o.rigidDLP = op.stokesDLmatrix(rigid);
       o.NearR2R = rigid.getZone([],1);
       if o.confined
	      [~,o.NearR2W] = rigid.getZone(walls,2);
	      [~,o.NearW2R] = walls.getZone(rigid,2);
       end
    end

    if o.withRigid
      [rhsRig] = o.sdcrigidrhs(VES,WALL,PROP,SDCPROP,ves,RIGID);
      [Xn,iflagTemp,R,I,resvec] = gmres(@(X) o.rigidmat(X,rigid(k)),rhsRig,[],o.gmresTol,(2*Nrd+3)*nvrd);
      [Xrig,etaRig,uRig,omgRig] = o.extractRHSrigid(Xn,rigid(k));
      if iflagTemp ~= 0
	      iflag = iflagTemp;
      end
      % if any of the gmres iterations fail, assign the failure to
      % flag to iflag for monitor to output
    end

    if o.confined
      VES.X = X; 
      VES.sig = sigma;
      VES.u =(X - deltaX(:,:,n))/dtimes(n);
      SDCPROP.deltaX = X;
      SDCPROP.deltaSig = sigma;
      [rhs] = o.sdcwallrhs(VES,WALL,PROP,SDCPROP,ves,RIGID);
      [Xn,iflagTemp,R,I,resvec] = gmres(@(X) o.wallmat(X,walls),...
	    rhs,[],o.gmresTol,o.gmresMaxIter,...
	    @o.preconditionerBDWall);
      [eta,RS] = o.extractRHSWall(Xn,Nbd,nvbd);
      if iflagTemp ~= 0
	iflag = iflagTemp;
      end
      % if any of the gmres iterations fail, assign the failure to
      % flag to iflag for monitor to output
    else
      eta = zeros(2*Nbd,nvbd);
      RS = zeros(3,nvbd);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% add collision resolveoc = curve;
    if(o.resolveCol)
      ticcol = tic;
      PROPCOL.viscCont = viscCont;
      PROPCOL.vesicle = ves;
      PROPCOL.N = N;
      PROPCOL.nv = nv;
      PROPCOL.Nbd = Nbd;
      PROPCOL.nvbd = nvbd;
      PROPCOL.X = X;
      PROPCOL.sigma = sigma;
      PROPCOL.u = u;
      PROPCOL.eta = eta;
      PROPCOL.RS = RS;
      
      PROPCOL.Nrd = Nrd;
      PROPCOL.nvrd = nvrd;
      PROPCOL.Xrig = Xrig;
      PROPCOL.etaRig = etaRig;
      PROPCOL.om = om;
      if o.withRigid
        PROPCOL.uRig = uRig;
        PROPCOL.omgRig = omgRig;
        PROPCOL.rhsRig = rhsRig; 
      end

      [X,sigma,u,eta,RS,Xrig,etaRig,fc_tot] = o.resolveCollisionGSrigid(Xprov(:,:,n+1),Xprov(:,:,n+1)+X,walls,rigid(k),PROPCOL);
      toccol = toc(ticcol);
      om.writeMessage('After resolve collision.','%s\n');
      timeCol = timeCol + toccol;
      %[X,sigma,u,eta,RS,fc_tot] = o.resolveCollision(Xprov(:,:,n+1),Xprov(:,:,n+1)+X,walls,PROPCOL);
      
      %if om.saveData
      %  fcfile = [om.dataFile 'fc'];
      %  fid = fopen(fcfile,'a');
      %  fwrite(fid,[time(:);fc_tot(:)],'double');
      %  fclose(fid);
      %end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of add collision resolve
    
    o.SDCcorrect = false;
    % turn correct off since it is not used when forming the
    % provisional solution


    deltaX(:,:,n+1) = X;
    deltaSigma(:,:,n+1) = sigma;
    deltaEta(:,:,n+1) = eta;
    deltaRS(:,:,n+1) = RS;
    deltaFC(:,:,n+1) = fc_tot;
    % approximations of the error
    uDeltaX = (deltaX(:,:,n+1) - deltaX(:,:,n))/dtimes(n);
    uProv(:,:,n+1) = uProv(:,:,n+1) + uDeltaX;
  end
  o.dt = dt;
  % go back to original time step
  
  Xprov = Xprov + deltaX;
  sigmaProv = sigmaProv + deltaSigma;
  etaProv = etaProv + deltaEta;
  RSprov = RSprov + deltaRS;
  FCprov = FCprov + deltaFC;
  % update provision solution
  
  if sdcCount < o.nsdc
    for n = 1:abs(o.orderGL)
      vesicle(n) = capsules(Xprov(:,:,n),sigmaProv(:,:,n),...
          [],kappa,viscCont,o.antiAlias);
      Galpert(:,:,:,n) = o.op.stokesSLmatrix(vesicle(n));
%      Gbarnett(:,:,:,n) = o.op.laplaceSLcomplexMatrix(vesicle(n));
      if any(viscCont ~= 1)
        D(:,:,:,n) = o.op.stokesDLmatrix(vesicle(n));
      else
        D = [];
      end
    end
    % update the vesicle objects and the single-layer potentials

    [vesVel,divVesVel,wallVel,residual] = ...
      o.computeVelocitiessdccol(FCprov,vesicle,Galpert,Gbarnett,D,walls,...
          etaProv,RSprov,NearV2V,NearV2W,NearW2V,NearW2W,uProv);
    om.writeMessage('After compute residual.','%s\n');
  end
  % if we are only recording the error in area and length, don't need
  % to compute the final residual

  normRes = max(abs(residual(:,:,end)));
  normRes = max(normRes);

  [~,a,l] = oc.geomProp(Xprov(:,:,end));
  ea = abs(a - a0)./abs(a0);
  el = abs(l - l0)./abs(l0);
  eaVec = [eaVec ea];
  elVec = [elVec el];
  message = ['sdcCount   ' num2str(sdcCount,'%2d') ...
    ': Residual = ' num2str(normRes(end),'%5.2e') ...
    ', eA = ' num2str(max(ea),'%5.2e') ...
    ', eL = ' num2str(max(el),'%5.2e')];
  om.writeMessage(message,'%s\n');
end
% End of doing SDC corrections
% update solution with an SDC iteration

if o.timeAdap
  [accept,dtScale] = o.newTimeStepSize(a,l,...
      aInit,lInit,accept,om);
  % if doing adaptive time stepping, get new time step
  % size and time step scaling
else
  accept = true;
  dtScale = 1;
  % if not doing adaptive time stepping, keep the same
  % time step size and always accept solution
end

if accept
  % take the new solution
  X = Xprov(:,:,end);
  sigma = sigmaProv(:,:,end);
  if o.nsdc > 0
    u = vesVel(:,:,end);
  else
    u = uProv(:,:,end);
  end
  eta = etaProv(:,:,end);
  RS = RSprov(:,:,end);
  FC = FCprov(:,:,end);
  om.writeMessage(['total time for collision resolving: ' num2str(timeCol)],'%s\n');
  if om.saveData
    fcfile = [om.dataFile 'fc'];
    fid = fopen(fcfile,'a');
    fwrite(fid,[time(:);FC(:)],'double');
    fclose(fid);
  end
else
  % revert to the old solution
  X = Xstore;
  sigma = sigStore;
  u = uStore;
  eta = etaStore;
  RS = RSstore;
  Xrig = XrigStore;
  etaRig = etaRigStore;
end

if accept && o.periodic
  X = oc.addAndRemove(X,walls,o.near,o.fmm);
end
end % GSSDCCOL


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sdcvesrhs
%function [rhs,vesicle,Xo,testen] = sdcvesrhs(o,VES,WALL,PROP,SDCPROP,RIGID)
function [rhs,vesicle,Xo] = sdcvesrhs(o,VES,WALL,PROP,SDCPROP,RIGID)
% VES.X, VES.sig, VES.u, WALL.eta, WALL.RS, WALL.walls, PROP.kappa,
% PROP.viscCont, PROP.updatePreco,
% SDCPROP.deltaX, SDCPROP.deltaSig, SDCPROP.deltaEta, SDCPROP.deltaRS
% SDCPROP.diffResidual, SDCPROP.vesVel, SDCPROP.wallVel, SDCPROP.sa,
% SDCPROP.IK.

if o.SDCcorrect
  deltaX = SDCPROP.deltaX;
  deltaSig = SDCPROP.deltaSig;
  deltaEta = SDCPROP.deltaEta;
  deltaRS = SDCPROP.deltaRS;
  diffResidual = SDCPROP.diffResidual;
  vesVel = SDCPROP.vesVel;
  wallVel = SDCPROP.wallVel;
  sa = SDCPROP.sa;
  IK = SDCPROP.IK;
end

Xstore = VES.X;
sigStore = VES.sig;
uStore = VES.u;

walls = WALL.walls;
etaStore = WALL.eta;
RSstore = WALL.RS;


rigid = RIGID.rigid;
etaRigStore = RIGID.etaRig;

viscCont = PROP.viscCont;
kappa = PROP.kappa;
updatePreco = PROP.updatePreco;
om = PROP.om;

N = size(Xstore,1)/2; % Number of points per vesicle
nv = size(Xstore,2); % Number of vesicles
if o.confined
  Xwalls = walls.X; % discretization points of solid walls
else
  Xwalls = [];
end
Nbd = size(Xwalls,1)/2; % Number of points on the solid walls
nvbd = size(Xwalls,2); % number of solid wall components
if o.withRigid
  Xrigid = rigid.X;
else
  Xrigid = [];
end
Nrd = size(Xrigid,1)/2;
nvrd = size(Xrigid,2);

alpha = (1 + viscCont)/2;
% constant that appears in front of time derivative in
% vesicle dynamical equations

Xm = zeros(2*N,nv);
sigmaM = zeros(N,nv);
uM = zeros(2*N,nv);
Xo = zeros(2*N,nv);
etaM = zeros(2*Nbd,nvbd);
RSm = zeros(3,nvbd);
etaRigM = zeros(2*Nrd,nvrd);
for k = 1:o.order
  Xm = Xm + Xstore(:,:,k)*o.Xcoeff(k);
  sigmaM = sigmaM + sigStore(:,:,k)*o.Xcoeff(k);
  uM = uM + uStore(:,:,k)*o.Xcoeff(k);
  if o.confined
    etaM = etaM + etaStore(:,:,k)*o.Xcoeff(k);
    RSm = RSm + RSstore(:,:,k)*o.Xcoeff(k);
  end
  if o.withRigid
    etaRigM = etaRigM + etaRigStore(:,:,k)*o.Xcoeff(k);
  end
  Xo = Xo + Xstore(:,:,k)*o.rhsCoeff(k);
end
% Form linear combinations of previous time steps needed
% for Ascher, Ruuth, and Wetton IMEX methods

vesicle = capsules(Xm,sigmaM,uM,kappa,viscCont,o.antiAlias);
% build an object vesicle that contains tangent vector, jacobian, etc.

op = o.op;
if ~o.SDCcorrect
  o.Galpert = op.stokesSLmatrix(vesicle);
  %o.Gbarnett = op.laplaceSLcomplexMatrix(vesicle);
end
% Build single layer potential matrix and put it in current object
% If we are doing an sdc update, this is already precomputed and 
% stored from when we formed the provisional solution

if ~o.SDCcorrect
  if any(viscCont ~= 1)
    o.D = op.stokesDLmatrix(vesicle);
  else
    o.D = [];
  end
end
% Compute double-layer potential matrix due to each vesicle
% independent of the others.  Matrix is zero if there is no
% viscosity contrast

if ~o.SDCcorrect
  if o.withRigid
    o.rigidDLP = op.stokesDLmatrix(rigid);
  end
end

if ~o.SDCcorrect
  if o.near
    if o.confined
      [o.NearV2V,o.NearV2W] = vesicle.getZone(walls,3);
      % Need vesicle to vesicle and vesicle to wall interactions
      if nvbd == 1 
        [o.NearW2W,o.NearW2V] = walls.getZone(vesicle,2);
      else
        [o.NearW2W,o.NearW2V] = walls.getZone(vesicle,3);
      end
      % Only need wall to vesicle interactions.  Wall to wall
      % interactions should also use near-singular integration since
      % they may be close to one another
    else
      o.NearV2V = vesicle.getZone([],1);
      % no solid walls, so only need vesicle-vesicle intearactions
      o.NearV2W = [];
      o.NearW2V = [];
      o.NearW2W = [];
    end
    if o.withRigid
      [~,o.NearV2R] = vesicle.getZone(rigid,2); 
      [o.NearR2R,o.NearR2V] = rigid.getZone(vesicle,3);
      if o.confined
        [~,o.NearR2W] = rigid.getZone(walls,2);
        [~,o.NearW2R] = walls.getZone(rigid,2);
      end
    end
  else
    o.NearV2V = [];
    o.NearV2W = [];
    o.NearW2V = [];
    o.NearW2W = [];
  end
  % If using near-singular integration, need structures for deciding who
  % is close, how close it is, who is closest, etc.
end
% Only form near-singular integration structure if not doing an SDC
% update.  Otherwise this was formed and saved when forming the
% provisional solution

if o.collDete
  o.lapDLP = op.laplaceDLmatrix(vesicle);
  % build matrix that evaluates double-layer laplace potential for each
  % vesicle

  oc = curve;
  [icollisionVes,icollisionWall] = ...
    oc.collision(vesicle,walls,o.NearV2V,o.NearV2W,o.fmm,o.near);
  % Check for collisions 
  if icollisionVes
    fprintf('VESICLES HAVE CROSSED\n')
    pause
  end
  if icollisionWall
    fprintf('VESICLES HAVE CROSSED SOLID WALL\n')
    pause
  end
end
% check for collisions

if ~o.SDCcorrect
  rhs1 = Xo;
  rhs2 = zeros(N,nv);
else
  if strcmp(o.vesves,'explicit')
    rhs1 = diffResidual;
  else
    rhs1 = deltaX + diffResidual;
  end
  if any(vesicle.viscCont ~= 1)
    z = zeros(2*N*nv,1);
    for k = 1:nv
      z(2*(k-1)*N+1:2*k*N) = rhs1(:,k);
    end
    if strcmp(o.vesves,'explicit')
      z = o.IminusDexp(z,vesicle);
    else
      z = o.IminusD(z,vesicle);
    end
    for k = 1:nv
      rhs1(:,k) = z(2*(k-1)*N+1:2*k*N)/alpha(k);
    end
  end
  if strcmp(o.vesves,'explicit')
    rhs1 = rhs1 + deltaX;
  end
  if strcmp(o.solver,'method1')
    rhs2 = ones(N,nv);
  else
    rhs2 = -vesicle.surfaceDiv(vesVel);
  end
end
% Parts of rhs from previous solution.  The right-hand-side depends on
% whether we are doing an SDC correction or forming the provisional
% solution.

% START TO COMPUTE RIGHT-HAND SIDE DUE TO VESICLE TRACTION JUMP
if o.gravity
  % For now, treat gravity explicit same as backgroud flow
  fgrav = o.gravityForce(Xm,vesicle,o.gCont);
else
  fgrav = zeros(2*N,nv);
%  itest = N/4+1;
%  ctest = 1000;
%  myrand = rand*1000;
%  fgrav(itest) = -1*(myrand+ctest)*vesicle.normal(itest);
%  fgrav(N+itest) = -1*(myrand+ctest)*vesicle.normal(N+itest);
%  Idx = ifftshift(-N/2:N/2+1);
%  Idx = abs(Idx) >= N/4;
%  fgrav = o.f_smooth(fgrav,N,1);
%  coeffx = fft(fgrav(1:N));
%  coeffy = fft(fgrav(N+1:2*N));
%  fcHighFreqNorm = sqrt(norm(coeffx(Idx))^2+norm(coeffy(Idx))^2)
end
if strcmp(o.vesves,'explicit')
  if o.profile
    tic
  end
  if ~o.SDCcorrect
    tictrac = tic;
    f = vesicle.tracJump(Xm,sigmaM) + PROP.fc_tot + fgrav;
    toctrac = toc(tictrac);
    om.writeMessage(['time for calculate traction jump: ' num2str(toctrac)],'%s\n');
  else
    f = vesicle.tracJump(deltaX,deltaSig);
  end
  if ~o.fmm
    kernel = @op.exactStokesSL;
  else
    kernel = @op.exactStokesSLfmm;
  end

  ticSLP = tic;
  if ~o.near
    Fslp = kernel(vesicle,f,[]);
    % Evaulate single-layer potential on all vesicles but itself without
    % near-singular integration
  else
    if o.nearStrat == 'cauchy'
      Fslp = op.nearSingStokesSLP(vesicle,f,vesicle,...
          o.Gbarnett,true,o.fmm);
    else
      SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);
      SLPtrap = SLP;
      kernelDirect = kernel;
      Fslp = op.nearSingInt(vesicle,f,SLP,SLPtrap,...
          o.NearV2V,kernel,kernelDirect,vesicle,true);
     end
    % Use near-singular integration to compute single-layer potential
    % due to all other vesicles.  Need to pass function op.exactStokesSL
    % so that the layer potential can be computed at far points and
    % Lagrange interpolation points
  end
  tocSLP = toc(ticSLP);
  om.writeMessage(['time for ves-ves explicit SLP rhs part: ' num2str(tocSLP)],'%s\n');
  if o.profile
    fprintf('Build right-hand side V2V and V2W   %5.1e\n',toc);
  end
end
for k = 1:nv
  rhs1(:,k) = rhs1(:,k) + o.dt/alpha(k)*Fslp(:,k);
end

if (o.gravity && ~o.SDCcorrect)
Gfg = op.exactStokesSLdiag(vesicle,o.Galpert,fgrav);
rhs1 = rhs1 + o.dt*Gfg*diag(1./alpha);
end
%if 1
%Gfg = op.exactStokesSLdiag(vesicle,o.Galpert,fgrav);
%testen = Gfg;
%coeffx = fft(testen(1:N));
%coeffy = fft(testen(N+1:2*N));
%GfcHighFreqNorm = sqrt(norm(coeffx(Idx))^2+norm(coeffy(Idx))^2)
%rhs1 = rhs1 + o.dt*Gfg*diag(1./alpha);
%testen = o.dt*Gfg*diag(1./alpha);
%coeffx = fft(testen(1:N));
%coeffy = fft(testen(N+1:2*N));
%abs(coeffy) < 1e-4;
%abs(coeffx) < 1e-4;
%sqrt(norm(coeffx(Idx))^2+norm(coeffy(Idx))^2);
%end

% END TO COMPUTE RIGHT-HAND SIDE DUE TO VESICLE TRACTION JUMP

% START TO COMPUTE RIGHT-HAND SIDE DUE TO VISCOSITY CONTRAST
if any(vesicle.viscCont ~= 1)
  if o.near
    jump = 1/2*(1-vesicle.viscCont);
    DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle,o.D,X);
  end
  % Need to add jump to double-layer potential if using near-singular
  % integration so that we can compute boundary values for the
  % near-singular integration algorithm

  if ~o.fmmDLP
    kernel = @op.exactStokesDL;
  else
    kernel = @op.exactStokesDLfmm;
  end
  if strcmp(o.vesves,'explicit')
    ticDLP = tic;
    if ~o.near
      Fdlp = -o.dt * kernel(vesicle,uM,[]);
      % Evaulate the velocity due to the viscosity contrast
      % on all vesicles but itself WITHOUT near-singular
      % integration
    else
      kernelDirect = kernel;
      Fdlp = -o.dt * op.nearSingInt(vesicle,uM,DLP,DLP,...
          o.NearV2V,kernel,kernelDirect,vesicle,true);
      % Evaulate the velocity due to the viscosity contrast on all
      % vesicles but itself WITH near-singular integration
    end
    tocDLP = toc(ticDLP);
    om.writeMessage(['time for ves-ves explicit DLP rhs part: ' num2str(tocDLP)],'%s\n');
  end
else
  Fdlp = zeros(2*N,nv);
  % If no viscosity contrast, there is no velocity induced due to a
  % viscosity contrast
end

if (any(viscCont ~= 1) && o.SDCcorrect && strcmp(o.vesves,'explicit'))
  DXo = op.exactStokesDLdiag(vesicle,o.D,deltaX);
  rhs1 = rhs1 - (Fdlp + DXo) * diag(1./alpha);
end
if (any(viscCont ~= 1) && ~o.SDCcorrect)
  DXo = op.exactStokesDLdiag(vesicle,o.D,Xo);
  rhs1 = rhs1 - (Fdlp + DXo) * diag(1./alpha);
end
% add in viscosity contrast term due to each vesicle independent of the
% others (o.D * Xo) from the previous solution followed by the term due
% to all other vesicles (Fdlp)

if o.adhesion
  adhesion = vesicle.adhesionTerm(o.adStrength,o.adRange);
  if ~o.fmm
    kernel = @op.exactStokesSL;
  else
    kernel = @op.exactStokesSLfmm;
  end

  Fadhesion = op.exactStokesSLdiag(vesicle,o.Galpert,adhesion);
  % diagonal term of adhesion

  if ~o.near
    Fadhesion = Fadhesion + kernel(vesicle,adhesion);
    % Evaulate single-layer potential on all vesicles but itself without
    % near-singular integration
  else
    SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);
    SLPtrap = SLP;
    kernelDirect = kernel;
    if o.nearStrat == 'cauchy'
      Fadhesion = Fadhesion + ...
          op.nearSingStokesSLP(vesicle,adhesion,vesicle,...
          o.Gbarnett,true,o.fmm);
    else
      Fadhesion = Fadhesion + ...
          op.nearSingInt(vesicle,adhesion,SLP,SLPtrap,...
          o.NearV2V,kernel,kernelDirect,vesicle,true);
    end
    % Use near-singular integration to compute single-layer potential
    % due to all other vesicles.  Need to pass function op.exactStokesSL
    % so that the layer potential can be computed at far points and
    % Lagrange interpolation points
  end

  rhs1 = rhs1 + o.dt*Fadhesion*diag(1./alpha);
end

% START TO COMPUTE RIGHT-HAND SIDE DUE TO SOLID WALLS
if o.confined
  if strcmp(o.vesves,'explicit')
    if o.profile
      tic
    end
    if o.near
      jump = -1/2;
      DLP = @(X) jump*X + op.exactStokesDLdiag(walls,o.wallDLP,X);
    end
    % compute the matrix for doing evaluating the double-layer potential
    % on the solid walls required for near-singular integration

    if ~o.fmmDLP
      kernel = @op.exactStokesDL;
    else
      kernel = @op.exactStokesDLfmm;
    end
    ticwall2ves = tic;
    if ~o.near
      if ~o.SDCcorrect
        charge = etaM;
      else
        charge = deltaEta;
      end
      [~,Fwall2Ves] = kernel(walls,charge,[],...
          vesicle.X,1:nvbd);
      % velocity field due to the walls evaluated on the vesicle
    else
      if ~o.SDCcorrect
        charge = etaM;
      else
        charge = deltaEta;
      end
      kernelDirect = kernel;
      Fwall2Ves = op.nearSingInt(walls,charge,DLP,DLP,...
          o.NearW2V,kernel,kernelDirect,vesicle,false);
    end
    tocwall2ves = toc(ticwall2ves);
    om.writeMessage(['time for ves explicit wall2ves rhs part: ' num2str(tocwall2ves)],'%s\n');

    for k = 2:nvbd
      if ~o.SDCcorrect
        stokeslet = RSm(1:2,k);
        rotlet = RSm(3,k);
      else
        stokeslet = deltaRS(1:2,k);
        rotlet = deltaRS(3,k);
      end
      Fwall2Ves = Fwall2Ves + ...
        o.RSlets(vesicle.X,walls.center(:,k),stokeslet,rotlet);
    end
    if o.profile
      fprintf('Build right-hand side W2V           %5.1e\n',toc);
    end

  end
else
  Fwall2Ves = zeros(2*N,nv);
  if ~o.SDCcorrect
    rhs1 = rhs1 + o.dt*o.farField(Xm)*diag(1./alpha);
    % Add in far-field condition (extensional, shear, etc.)
  end
end
rhs1 = rhs1 + o.dt*Fwall2Ves*diag(1./alpha);
% right-hand side of the velocity evaluated on the solid walls
% END TO COMPUTE RIGHT-HAND SIDE DUE TO SOLID WALLS


% START TO COMPUTE RIGHT-HAND SIDE DUE TO RIGIDBODY
if o.withRigid
  if strcmp(o.vesves,'explicit')
    if o.profile
      tic
    end
    if o.near
      jump = -1/2;
      DLP = @(X) jump*X + op.exactStokesDLdiag(rigid,o.rigidDLP,X);
      %DLP = @(X) jump*X + op.exactStokesRigidDLdiag(rigid,o.rigidDLP,theta,X);
    end
    % compute the matrix for doing evaluating the double-layer potential
    % on the solid walls required for near-singular integration

    if ~o.fmmDLP
      kernel = @op.exactStokesDL;
    else
      kernel = @op.exactStokesDLfmm;
    end
    if ~o.near
      if ~o.SDCcorrect
        charge = etaRigM;
      else
        charge = deltaRigEta;
      end
      [~,Frigid2Ves] = kernel(rigid,charge,...
          vesicle.X,1:nvrd);
      % velocity field due to the walls evaluated on the vesicle
    else
      if ~o.SDCcorrect
        charge = etaRigM;
      else
        charge = deltaRigEta;
      end
      kernelDirect = kernel;
      Frigid2Ves = op.nearSingInt(rigid,charge,DLP,DLP,...
          o.NearR2V,kernel,kernelDirect,vesicle,false);
    end

    for k = 1:nvrd
      if ~o.SDCcorrect
        [stokeslet, rotlet] = o.getRS(rigid.X(:,k),rigid.sa(:,k),rigid.center(:,k),etaRigM(:,k));
      else
        [stokeslet, rotlet] = o.getRS(rigid.X(:,k),rigid.sa(:,k),rigid.center(:,k),deltaRigEta(:,k));
      end
      Frigid2Ves = Frigid2Ves + ...
        o.RSlets(vesicle.X,rigid.center(:,k),stokeslet,rotlet);
    end
    if o.profile
      fprintf('Build right-hand side W2V           %5.1e\n',toc);
    end

  end
else
  Frigid2Ves = zeros(2*N,nv);
end
rhs1 = rhs1 + o.dt*Frigid2Ves*diag(1./alpha);
% right-hand side of the velocity evaluated on the solid walls
% END TO COMPUTE RIGHT-HAND SIDE DUE TO RIGIDBODY

% START TO COMPUTE THE RIGHT-HAND SIDE FOR THE INEXTENSIBILITY CONDITION
if (strcmp(o.solver,'method1'))
  % If using method1, vesicle-vesicle interactions and the presence of a
  % viscosity contrast is irrelevent
  if ~o.SDCcorrect
    rhs2 = rhs2 + vesicle.surfaceDiv(Xo); 
  else
    divf = zeros(N,nv);
    for k = 1:nv
      divf(:,k) = curve.arcDeriv(Xm(1:N,k),1,1./sa(:,k),...
          IK(:,k)).^2 + ...
                  curve.arcDeriv(Xm(N+1:2*N,k),1,1./sa(:,k),...
          IK(:,k)).^2; 
    end
    rhs2 = 1/2*(rhs2 - divf);
  end
else 
  % If using method2, method3, or method4, vesicle-vesicle interaction
  % affects the right-hand side of the inextensibility condition

  if ~o.confined && ~o.SDCcorrect
    if any(viscCont ~= 1)
      if o.profile
        tic
      end
      rhs2 = rhs2 - vesicle.surfaceDiv(...
          o.solveIminusD(o.farField(Xm),vesicle));
      if o.profile
        fprintf('Solve system alpha*I - DLP          %5.1e\n',toc);
      end
    else
      rhs2 = rhs2 - vesicle.surfaceDiv(o.farField(Xm));
    end
  end
  % add in term from farfield

  if (strcmp(o.vesves,'explicit'))
    rhs2 = rhs2 - vesicle.surfaceDiv(Fslp); 
    rhs2 = rhs2 - vesicle.surfaceDiv(Fwall2Ves);
    if o.adhesion
      rhs2 = rhs2 + vesicle.surfaceDiv(Fadhesion);
    end
  end
end
% rhs2 is the right-hand side for the inextensibility condition
% END TO COMPUTE THE RIGHT-HAND SIDE FOR THE INEXTENSIBILITY CONDITION

% START TO COMPUTE RIGHT-HAND SIDE DUE TO NEW BENDING TERM
if o.bending
  f = vesicle.newBending;
  % New part of traction jump if we have variable bending
  rhs1 = rhs + o.dt * op.exactStokesSLdiag(vesicle,o.Galpert,f);
end
% This term is always treated fully explicitly.  The highest derivative
% it contains is second-order and it appears through a curvature
% END TO COMPUTE RIGHT-HAND SIDE DUE TO NEW BENDING TERM.  THIS WILL
% ONLY SUPPORT A SINGLE VESICLE IN AN UNBOUNDED FLOW


rhs = [rhs1; rhs2];
rhs = rhs(:);
% Stack the right-hand sides in an alternating with respect to the

usePreco = true;
%usePreco = false;

% START BUILDING BLOCK-DIAGONAL PRECONDITIONER
if usePreco
  if updatePreco
    % only build the preconditioner if updatePreco == true
    if o.profile
      tic
    end
    %[Ben,Ten,Div] = vesicle.computeDerivs;
    ticdiffop = tic;
    [Ben1,Ten1,Div1] = o.computeUpDerivs(vesicle);
    [Ben,Ten,Div] = o.upSampleDerivs(Ben1,Ten1,Div1);
    tocdiffop = toc(ticdiffop);
    om.writeMessage(['time for build diff op for precond: ' num2str(tocdiffop)],'%s\n');
    if o.profile
      fprintf('Build differential operators        %5.1e\n',toc);
    end
    % Compute bending, tension, and surface divergence of current
    % vesicle configuration

    if o.profile
      tic
    end

    if (strcmp(o.solver,'method1') || strcmp(o.solver,'method2'))
      bdiagVes.L = zeros(3*N,3*N,nv);
      bdiagVes.U = zeros(3*N,3*N,nv);
    elseif strcmp(o.solver,'method3')
      bdiagVes.GT = zeros(2*N,N,nv); % SLP of tension
      bdiagVes.DGT = zeros(N,N,nv); % divergence of SLP of tension
      bdiagVes.DGB = zeros(N,2*N,nv); % divergence of SLP of bending
      bdiagVes.schur = zeros(2*N,2*N,nv); 
      % schur complement of the lower right N by N block
    elseif strcmp(o.solver,'method4')
      bdiagVes.GT = zeros(2*N,N,nv); % SLP of tension
      bdiagVes.IpBen = zeros(2*N,2*N,nv); 
      % inverse of identity plus SLP of Bending
      bdiagVes.DGB = zeros(N,2*N,nv); % divergence of SLP of bending
      bdiagVes.schur = zeros(N,N,nv); 
      % schur complement of the upper left 2*N by 2*N block
    end

    for k=1:nv
      if strcmp(o.solver,'method1')
        if any(vesicle.viscCont ~= 1)
          [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu(...
            [o.beta*(eye(2*N) - o.D(:,:,k)/alpha(k)) + ...
                o.dt/alpha(k)*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
            -o.dt/alpha(k)*o.Galpert(:,:,k)*Ten(:,:,k); ...
            o.beta*Div(:,:,k) zeros(N)]);
        else
          [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu(...
            [o.beta*eye(2*N) + ...
                o.dt*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
            -o.dt/alpha(k)*o.Galpert(:,:,k)*Ten(:,:,k); ...
            o.beta*Div(:,:,k) zeros(N)]);
        end
      elseif strcmp(o.solver,'method2')
        if any(vesicle.viscCont ~= 1)
          [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu( ...
            [o.beta*(eye(2*N) - o.D(:,:,k)/alpha(k)) + ...
                o.dt/alpha(k)*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
            -o.dt/alpha(k)*o.Galpert(:,:,k)*Ten(:,:,k); ...
            -vesicle.kappa*Div(:,:,k)*(alpha(k)*eye(2*N) - o.D(:,:,k))\...
                (o.Galpert(:,:,k)*Ben(:,:,k)) ...
            Div(:,:,k)*(alpha(k)*eye(2*N) - o.D(:,:,k))\...
                (o.Galpert(:,:,k)*Ten(:,:,k))]);
        else
          [bdiagVes.L(:,:,k),bdiagVes.U(:,:,k)] = lu(...
            [o.beta*eye(2*N) + ...
                o.dt*vesicle.kappa*o.Galpert(:,:,k)*Ben(:,:,k) ...
            -o.dt*o.Galpert(:,:,k)*Ten(:,:,k); ...
            -vesicle.kappa*Div(:,:,k)*o.Galpert(:,:,k)*Ben(:,:,k) ...
            Div(:,:,k)*o.Galpert(:,:,k)*Ten(:,:,k)]);
        end
      elseif strcmp(o.solver,'method3')
        % schur complement of lower right block of method2
        bdiagVes.GT(:,:,k) = o.Galpert(:,:,k)*Ten(:,:,k);
        bdiagVes.DGT(:,:,k) = (Div(:,:,k)*o.Galpert(:,:,k)*...
            Ten(:,:,k))\eye(N);
        bdiagVes.DGB(:,:,k) = ...
          vesicle.kappa*Div(:,:,k)*o.Galpert(:,:,k)*Ben(:,:,k);
        bdiagVes.schur(:,:,k) = ...
          inv((o.beta*eye(2*N) + vesicle.kappa*o.dt*...
          o.Galpert(:,:,k)*Ben(:,:,k)) - ...
          o.dt*o.Galpert(:,:,k)*Ten(:,:,k)*bdiagVes.DGT(:,:,k)*...
          vesicle.kappa*Div(:,:,k)*o.Galpert(:,:,k)*Ben(:,:,k));

      elseif strcmp(o.solver,'method4')
        % schur complement of upper left block of method2
        bdiagVes.GT(:,:,k) = o.Galpert(:,:,k)*Ten(:,:,k);
        bdiagVes.IpBen(:,:,k) = inv(o.beta*eye(2*N) + ...
            vesicle.kappa*o.dt*o.Galpert(:,:,k)*Ben(:,:,k));
        bdiagVes.DGB(:,:,k) = ...
          vesicle.kappa*Div(:,:,k)*o.Galpert(:,:,k)*Ben(:,:,k);
        bdiagVes.schur(:,:,k) = ...
          inv(Div(:,:,k)*(-vesicle.kappa*o.dt*o.Galpert(:,:,k)*...
          Ben(:,:,k)*bdiagVes.IpBen(:,:,k) + eye(2*N))*...
          bdiagVes.GT(:,:,k));
      end
    end
    o.bdiagVes = bdiagVes;
    % Build block-diagonal preconditioner of self-vesicle 
    % intearctions in matrix form
    if o.withRigid
      bdiagRigid.L = zeros(2*Nrd+3,2*Nrd+3,nvrd);
      bdiagRigid.U = zeros(2*Nrd+3,2*Nrd+3,nvrd);
      for k=1:nvrd
        cmR = rigid.center(:,k);
        Xr = rigid.X(:,k);
        saR = rigid.sa(:,k);
        [bdiagRigid.L(:,:,k),bdiagRigid.U(:,:,k)] = ... 
	lu([0.5*eye(2*Nrd)-o.rigidDLP(:,:,k) ... 
	[ones(Nrd,1);zeros(Nrd,1)] [zeros(Nrd,1);ones(Nrd,1)] ... 
	[Xr(Nrd+1:end)-cmR(2);-(Xr(1:Nrd)-cmR(1))];...
        [saR'*2*pi/Nrd zeros(1,Nrd) 0 0 0];[zeros(1,Nrd) saR'*2*pi/Nrd 0 0 0];...
	[2*pi/Nrd*(Xr(Nrd+1:end)'-cmR(2)).*saR' -2*pi/Nrd*(Xr(1:Nrd)'-cmR(1)).*saR' 0 0 0]]);
      end
      o.bdiagRigid = bdiagRigid;
    end

    if o.profile
      fprintf('Build block-diagonal preconditioner %5.1e\n\n',toc);
    end
  end % updatePreco
end % usePreco
end % sdcvesrhs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sdcwallrhs
function [rhs] = sdcwallrhs(o,VES,WALL,PROP,SDCPROP,vesicle,RIGID)

Xstore = VES.X;
sigStore = VES.sig;
uStore = VES.u;

walls = WALL.walls;

rigid = RIGID.rigid;
etaRigStore = RIGID.etaRig;

viscCont = PROP.viscCont;
om = PROP.om;

if o.SDCcorrect
  if size(SDCPROP) 
    deltaX = SDCPROP.deltaX;
    deltaSig = SDCPROP.deltaSig;
    wallVel = SDCPROP.wallVel;
  else
    deltaX = VES.X;
    deltaSig = VES.sig;
    wallVel = walls.u;
  end
end

N = size(Xstore,1)/2; % Number of points per vesicle
nv = size(Xstore,2); % Number of vesicles
if o.confined
  Xwalls = walls.X; % discretization points of solid walls
else
  Xwalls = [];
end
Nbd = size(Xwalls,1)/2; % Number of points on the solid walls
nvbd = size(Xwalls,2); % number of solid wall components
if o.withRigid
  Xrigid = rigid.X;
else
  Xrigid = [];
end
Nrd = size(Xrigid,1)/2;
nvrd = size(Xrigid,2);

alpha = (1 + viscCont)/2; 
% constant that appears in front of time derivative in
% vesicle dynamical equations

Xm = zeros(2*N,nv);
sigmaM = zeros(N,nv);
uM = zeros(2*N,nv);
Xo = zeros(2*N,nv);
etaRigM = zeros(2*Nrd,nvrd);
for k = 1:o.order
  Xm = Xm + Xstore(:,:,k)*o.Xcoeff(k);
  sigmaM = sigmaM + sigStore(:,:,k)*o.Xcoeff(k);
  uM = uM + uStore(:,:,k)*o.Xcoeff(k);
  Xo = Xo + Xstore(:,:,k)*o.rhsCoeff(k);
  if o.withRigid
    etaRigM = etaRigM + etaRigStore(:,:,k)*o.Xcoeff(k);
  end
end
% Form linear combinations of previous time steps needed
% for Ascher, Ruuth, and Wetton IMEX methods

op = o.op;

if o.collDete
  o.lapDLP = op.laplaceDLmatrix(vesicle);
  % build matrix that evaluates double-layer laplace potential for each
  % vesicle

  oc = curve;
  [icollisionVes,icollisionWall] = ...
    oc.collision(vesicle,walls,o.NearV2V,o.NearV2W,o.fmm,o.near);
  % Check for collisions 
  if icollisionVes
    fprintf('VESICLES HAVE CROSSED\n')
    pause
  end
  if icollisionWall
    fprintf('VESICLES HAVE CROSSED SOLID WALL\n')
    pause
  end
end
% check for collisions

if ~o.SDCcorrect
  if o.confined
    rhs3 = walls.u;
  else
    rhs3 = [];
  end
else
  if o.confined
    rhs3 = -wallVel + walls.u;
  else
    rhs3 = [];
  end
end
% Parts of rhs from previous solution.  The right-hand-side depends on
% whether we are doing an SDC correction or forming the provisional
% solution.

if nv
% START TO COMPUTE RIGHT-HAND SIDE DUE TO VESICLE TRACTION JUMP
if o.gravity
  % For now, treat gravity explicit same as backgroud flow
  fgrav = o.gravityForce(Xm,vesicle,o.gCont);
else
  fgrav = zeros(2*N,nv);
end
if strcmp(o.vesves,'explicit')
  if o.profile
    tic
  end
  if ~o.SDCcorrect
    f = vesicle.tracJump(Xm,sigmaM) + PROP.fc_tot + fgrav;
  else
    f = vesicle.tracJump(deltaX,deltaSig);
  end
  if ~o.fmm
    kernel = @op.exactStokesSL;
  else
    kernel = @op.exactStokesSLfmm;
  end

  ticwallslp = tic;
  if ~o.near
    if o.confined
      [~,FSLPwall] = kernel(vesicle,f,[],walls.X,(1:nv));
      % Evaluate single-layer potential on solid walls due to all
      % vesicles
    else
      FSLPwall = [];
    end
  else
    if o.nearStrat == 'cauchy'
    else
      SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);
      SLPtrap = SLP;
      kernelDirect = kernel;
    end
    if o.confined
      if o.nearStrat == 'cauchy'
        FSLPwall = op.nearSingStokesSLP(vesicle,f,walls,...
            o.Gbarnett,false,o.fmm);
      else
        FSLPwall = op.nearSingInt(vesicle,f,SLP,SLPtrap,...
            o.NearV2W,kernel,kernelDirect,walls,false);
        % Evaluate the velocity on the walls due to the vesicles
      end
    else
      FSLPwall = [];
    end
  end
  tocwallslp = toc(ticwallslp);
  om.writeMessage(['time for wall slp: ' num2str(tocwallslp)],'%s\n');
  if o.profile
    fprintf('Build right-hand side V2V and V2W   %5.1e\n',toc);
  end
end
rhs3 = rhs3 - FSLPwall;

% START TO COMPUTE RIGHT-HAND SIDE DUE TO VISCOSITY CONTRAST
if any(vesicle.viscCont ~= 1)
  if o.near
    jump = 1/2*(1-vesicle.viscCont);
    DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle,o.D,X);
  end
  % Need to add jump to double-layer potential if using near-singular
  % integration so that we can compute boundary values for the
  % near-singular integration algorithm

  if ~o.fmmDLP
    kernel = @op.exactStokesDL;
  else
    kernel = @op.exactStokesDLfmm;
  end
  if strcmp(o.vesves,'explicit')
    ticwalldlp = tic;
    if ~o.near
      if o.confined
        [~,FDLPwall] = kernel(vesicle,uM,[],walls.X,(1:nv));
        % Evaulate the velocity due to the viscosity contrast on the
        % walls WITHOUT near-singulation integration
      else
        FDLPwall = [];
      end
    else
      kernelDirect = kernel;

      if o.confined
        FDLPwall = op.nearSingInt(vesicle,uM,DLP,DLP,...
            o.NearV2W,kernel,kernelDirect,walls,false);
        % Evaulate the velocity due to the viscosity contrast on the
        % walls WITH near-singulation integration
      else
        FDLPwall = [];
      end
    end
    tocwalldlp = toc(ticwalldlp);
    om.writeMessage(['time for wall dlp: ' num2str(tocwalldlp)],'%s\n');
  end
else
  FDLPwall = zeros(2*Nbd,nvbd);
  % If no viscosity contrast, there is no velocity induced due to a
  % viscosity contrast
end
rhs3 = rhs3 - FDLPwall;
% compute the double-layer potential due to all other vesicles from the
% appropriate linear combination of previous time steps.  Depends on
% time stepping order and vesicle-vesicle discretization
end

% START TO COMPUTE RIGHT-HAND SIDE DUE TO RIGIDBODY
if o.withRigid
  if strcmp(o.vesves,'explicit')
    if o.profile
      tic
    end
    if o.near
      jump = -1/2;
      DLP = @(X) jump*X + op.exactStokesDLdiag(rigid,o.rigidDLP,X);
      %DLP = @(X) jump*X + op.exactStokesRigidDLdiag(rigid,o.rigidDLP,theta,X);
    end
    % compute the matrix for doing evaluating the double-layer potential
    % on the solid walls required for near-singular integration

    if ~o.fmmDLP
      kernel = @op.exactStokesDL;
    else
      kernel = @op.exactStokesDLfmm;
    end
    if ~o.near
      if ~o.SDCcorrect
        charge = etaRigM;
      else
        charge = deltaRigEta;
      end
      [~,Frigid2Wall] = kernel(rigid,charge,...
          walls.X,1:nvrd);
      % velocity field due to the walls evaluated on the vesicle
    else
      if ~o.SDCcorrect
        charge = etaRigM;
      else
        charge = deltaRigEta;
      end
      kernelDirect = kernel;
      Frigid2Wall = op.nearSingInt(rigid,charge,DLP,DLP,...
          o.NearR2W,kernel,kernelDirect,walls,false);
    end

    for k = 1:nvrd
      if ~o.SDCcorrect
        [stokeslet, rotlet] = o.getRS(rigid.X(:,k),rigid.sa(:,k),rigid.center(:,k),etaRigM(:,k));
      else
        [stokeslet, rotlet] = o.getRS(rigid.X(:,k),rigid.sa(:,k),rigid.center(:,k),deltaRigEta(:,k));
      end
      Frigid2Wall = Frigid2Wall + ...
        o.RSlets(walls.X,rigid.center(:,k),stokeslet,rotlet);
    end
    if o.profile
      fprintf('Build right-hand side R2W           %5.1e\n',toc);
    end

  end
else
  Frigid2Wall = zeros(2*Nbd,nvbd);
end
rhs3 = rhs3 - Frigid2Wall;

rhs = [rhs3(:)];
rhs = [rhs; zeros(3*(nvbd-1),1)];
% Rotlet and Stokeslet equations

end % sdcwallrhs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vesmat
function [val] = vesmat(o,Xn,vesicle)
global matvecs

matvecs = matvecs + 1;
% counter for the number of matrix-vector multiplications
% that are required for the entire simulation

op = o.op; % poten class
N = vesicle.N; % Number of points
nv = vesicle.nv; % Number of vesicles

valPos = zeros(2*N,nv);
% right-hand side that corresponds to position equation
valTen = zeros(N,nv);
% right-hand side that corresponds to inextensibilty equation

Xm = zeros(2*N,nv);
sigmaM = zeros(N,nv);
for k=1:nv
  Xm(1:2*N,k) = Xn((3*k-3)*N+1:(3*k-1)*N);
  sigmaM(:,k) = Xn((3*k-1)*N+1:3*k*N);
end
% Unstack the position and tension from the input


f = vesicle.tracJump(Xm,sigmaM);
% f is the traction jump stored as a 2N x nv matrix

alpha = (1+vesicle.viscCont)/2; 
% constant that multiplies the time derivative in the 
% vesicle position equation

Gf = op.exactStokesSLdiag(vesicle,o.Galpert,f);
% Gf is the single-layer potential applied to the traction jump. 

if any(vesicle.viscCont ~= 1)
  DXm = op.exactStokesDLdiag(vesicle,o.D,Xm);
else
  DXm = zeros(2*N,nv);
end
% DXm is the double-layer potential applied to the position




% START OF EVALUATING VELOCITY ON VESICLES
valPos = valPos - o.dt*Gf*diag(1./alpha);
% self-bending and self-tension terms
valPos = valPos - o.beta*DXm*diag(1./alpha);
% self-viscosity contrast term
% END OF EVALUATING VELOCITY ON VESICLES


% START OF EVALUATING INEXTENSIBILITY CONDITION
if (strcmp(o.solver,'method1'))
  valTen = o.beta * vesicle.surfaceDiv(Xm);
  % compute surface divergence of the current GMRES iterate
  % method1 sets this equal to the surface divergence of
  % the previous time step
else
  if any(vesicle.viscCont ~= 1)
    if o.profile
      tic
    end
    valTen = vesicle.surfaceDiv(...
        o.solveIminusD(Gf+Fslp+Fwall2Ves+LetsVes,vesicle));
    if o.profile
      fprintf('Solve system alpha*I - DLP          %5.1e\n',toc);
    end
  else
    valTen = -1/o.dt*vesicle.surfaceDiv(valPos);
  end
  % method2, method3, and method4 sets the surface divergence of the sum
  % of single-layer potentials due to bending and tension plus the
  % farField to zero.  The only difference between the two methods is
  % how the preconditioner is used.  method2 uses the possibly
  % ill-conditioned full matrix where as method3 and method 4 use the
  % two schur complements.  Eventually will phase out method2 and it
  % will be fully replaced by method3 and method4
end
% Two possible discretizations of the inextensibility condition
% END OF EVALUATING INEXTENSIBILITY CONDITION

valPos = valPos + o.beta*Xm;
% beta times solution coming from time derivative

val = zeros(3*N*nv,1);
% Initialize output from vesicle and inextensibility equations to zero
for k=1:nv
  val((k-1)*3*N+1:3*k*N) = [valPos(:,k);valTen(:,k)];
end
% Stack val as [xcoordinate;ycoordinate;tension] repeated
% nv times for each vesicle
end % vesmat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vesmatSingle
function [val] = vesmatSingle(o,Xn,vesicle,i)
op = o.op; % poten class
N = vesicle.N; % Number of points
nv = vesicle.nv; % Number of vesicles

valPos = zeros(2*N,nv);
% right-hand side that corresponds to position equation
valTen = zeros(N,nv);
% right-hand side that corresponds to inextensibilty equation

Xm = zeros(2*N,nv);
sigmaM = zeros(N,nv);
for k=1:nv
  Xm(1:2*N,k) = Xn((3*k-3)*N+1:(3*k-1)*N);
  sigmaM(:,k) = Xn((3*k-1)*N+1:3*k*N);
end
% Unstack the position and tension from the input


f = vesicle.tracJump(Xm,sigmaM);
% f is the traction jump stored as a 2N x nv matrix

alpha = (1+vesicle.viscCont)/2; 
% constant that multiplies the time derivative in the 
% vesicle position equation

Gf = op.exactStokesSLdiag(vesicle,o.Galpert(:,:,i),f);
% Gf is the single-layer potential applied to the traction jump. 

if any(vesicle.viscCont ~= 1)
  DXm = op.exactStokesDLdiag(vesicle,o.D(:,:,i),Xm);
else
  DXm = zeros(2*N,nv);
end
% DXm is the double-layer potential applied to the position

% START OF EVALUATING VELOCITY ON VESICLES
valPos = valPos - o.dt*Gf*diag(1./alpha);
% self-bending and self-tension terms
valPos = valPos - o.beta*DXm*diag(1./alpha);
% self-viscosity contrast term
% END OF EVALUATING VELOCITY ON VESICLES


% START OF EVALUATING INEXTENSIBILITY CONDITION
if (strcmp(o.solver,'method1'))
  valTen = o.beta * vesicle.surfaceDiv(Xm);
  % compute surface divergence of the current GMRES iterate
  % method1 sets this equal to the surface divergence of
  % the previous time step
else
  if any(vesicle.viscCont ~= 1)
    if o.profile
      tic
    end
    valTen = vesicle.surfaceDiv(...
        o.solveIminusD(Gf+Fslp+Fwall2Ves+LetsVes,vesicle));
    if o.profile
      fprintf('Solve system alpha*I - DLP          %5.1e\n',toc);
    end
  else
    valTen = -1/o.dt*vesicle.surfaceDiv(valPos);
  end
  % method2, method3, and method4 sets the surface divergence of the sum
  % of single-layer potentials due to bending and tension plus the
  % farField to zero.  The only difference between the two methods is
  % how the preconditioner is used.  method2 uses the possibly
  % ill-conditioned full matrix where as method3 and method 4 use the
  % two schur complements.  Eventually will phase out method2 and it
  % will be fully replaced by method3 and method4
end
% Two possible discretizations of the inextensibility condition
% END OF EVALUATING INEXTENSIBILITY CONDITION

valPos = valPos + o.beta*Xm;
% beta times solution coming from time derivative

val = zeros(3*N*nv,1);
% Initialize output from vesicle and inextensibility equations to zero
for k=1:nv
  val((k-1)*3*N+1:3*k*N) = [valPos(:,k);valTen(:,k)];
end
% Stack val as [xcoordinate;ycoordinate;tension] repeated
% nv times for each vesicle
end % vesmatSingle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wallmat
function [val] = wallmat(o,Xn,walls)

op = o.op; % poten class
Nbd = walls.N; % Number of points on walls
nvbd = walls.nv; % Number of components to walls

valWalls = zeros(2*Nbd,nvbd);
% right-hand side that corresponds to solid wall equation
valLets = zeros(3*(nvbd-1),1);
% right-hand side corresponding to the rotlets and stokeslets


eta = Xn(1:end);
etaM = zeros(2*Nbd,nvbd);
for k = 1:nvbd
  etaM(:,k) = eta((k-1)*2*Nbd+1:2*k*Nbd);
end
% Unstack the density function from the input
otlets = Xn(2*nvbd*Nbd+1:end);
% stokeslets and rotlets of each vesicle.  Ordered as
% [stokeslet1(component 1);stokeslet1(component 2);rotlet1;...
%  stokeslet2(component 1);stokeslet2(component 2);rotlet2;...];

% START OF EVALUATING WALL TO WALL INTERACTIONS
if (o.confined && nvbd > 1)
  if o.profile
    tic
  end
  % only need to do wall to wall interactions if the domain is multiply
  % connected
  if ~o.fmmDLP
    kernel = @op.exactStokesDL;
    FDLPwall2wall = kernel(walls,etaM,[]);
  else
    kernel = @op.exactStokesDLfmm;
    FDLPwall2wall = zeros(2*walls.N,walls.nv);
    for k = 1:nvbd
      isou = [(1:k-1) (k+1:nvbd)];
      [~,FDLPwall2wall(:,k)] = kernel(walls,etaM,[],walls.X(:,k),isou);
    end
    % Since the double-layer potential is still expensive even with the
    % FMM, this eliminates the need to do one big FMM followed by a
    % bunch of small ones to subtract off the self-interaction term
    % which is calculated using the precomputed matrix
  end

  if o.profile
    fprintf('Apply W2W interaction               %5.1e\n\n',toc)
  end
end
% END OF EVALUATING WALL TO WALL INTERACTIONS

% START OF EVALUATING POTENTIAL DUE TO STOKESLETS AND ROTLETS
if nvbd > 1
  LetsWalls = zeros(2*Nbd,nvbd);
  for k = 2:nvbd
    stokeslet = otlets(3*(k-2)+1:3*(k-2)+2);
    rotlet = otlets(3*(k-1));
    LetsWalls = LetsWalls + o.RSlets(walls.X,walls.center(:,k),...
        stokeslet,rotlet);
    % compute velocity due to rotlets and stokeslets on the solid walls
  end
  valLets = o.letsIntegrals(otlets,etaM,walls);
  % Integral constraints on the density function eta related
  % to the weights of the stokeslets and rotlets
else
  LetsWalls = zeros(2*Nbd,nvbd);
  FDLPwall2wall = zeros(2*Nbd,nvbd);
end
% END OF EVALUATING POTENTIAL DUE TO STOKESLETS AND ROTLETS

% START OF EVALUATING VELOCITY ON WALLS
if o.confined
  valWalls = valWalls - 1/2*etaM + ...
    op.exactStokesDLdiag(walls,o.wallDLP,etaM);
  valWalls(:,1) = valWalls(:,1) + ...
    op.exactStokesN0diag(walls,o.wallN0,etaM(:,1));
end
% evaluate velocity on solid walls due to the density function.
% self solid wall interaction

valWalls = valWalls + FDLPwall2wall;
% velocity on walls due to all other walls
valWalls = valWalls + LetsWalls;
% velocity on walls due to the rotlets and stokeslets
% END OF EVALUATING VELOCITY ON WALLS

val = [valWalls(:);valLets];
% Stack velocity along the solid walls in same manner as above
% Stack the stokeslets and rotlet componenets at the end

end % wallmat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extractRHS
function [X,sigma,u,eta,RS] = extractRHS(o,Xn,N,nv,Nbd,nvbd,Xo)
X = zeros(2*N,nv);
sigma = zeros(N,nv);
eta = zeros(2*Nbd,nvbd);
RS = zeros(3,nvbd);
% allocate space for positions, tension, and density function

for k=1:nv
  X(:,k) = Xn((3*k-3)*N+1:(3*k-1)*N);
  sigma(:,k) = Xn((3*k-1)*N+1:3*k*N);
end
% unstack the positions and tensions

% Unstack the positions and tensions
Xn = Xn(3*nv*N+1:end);
for k = 1:nvbd
  eta(:,k) = Xn((k-1)*2*Nbd+1:2*k*Nbd);
end
% unstack the density function
otlets = Xn(2*nvbd*Nbd+1:end);
for k = 2:nvbd
  istart = (k-2)*3+1;
  iend = 3*(k-1);
  RS(:,k) = otlets(istart:iend);
end
% unstack the rotlets and stokeslets

u = (o.beta*X - Xo)/o.dt;
% Compute the velocity using the differencing stencil

end % extractRHS


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ns,totalnv,Xstart,Xend,totalPts] = preColCheck(o,X0,X1,walls,upSampleFactor)
N = size(X0,1)/2;
nv = size(X0,2);
%Ns = ones(nv,1)*N*upSampleFactor;
Ns = ones(nv,1)*N;
totalnv = nv;

% upsample positions
Xv0 = reshape(X0,N,2*nv);
%Xv0 = interpft(Xv0,N,1);
%Xv0 = interpft(Xv0,N*upSampleFactor,1);
Xv1 = reshape(X1,N,2*nv);
%Xv1 = interpft(Xv1,N*upSampleFactor,1);
%Xv1 = interpft(Xv1,N,1);

Xstart = Xv0(:);
Xend = Xv1(:);

if(size(walls,1))
  Ns = [Ns;ones(walls.nv,1)*walls.N*upSampleFactor];
  totalnv = totalnv + walls.nv;
  nb = walls.nv;
  Npb = walls.N;
  Xb = reshape(walls.X,Npb,2*nb);
  Xb = interpft(Xb,Npb*upSampleFactor,1);
  Xstart = [Xstart;Xb(:)];
  Xend = [Xend;Xb(:)];
end

totalPts = sum(Ns);

end % preColCheck

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = downSampleForce()



end % downSampleForce

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,jaco,ivs,listnv,jacoSmooth,vtoiv] = preprocessRigid(o,vgrad,ids,vols,G,vblock,N,nv,Nrd,nvrd,alpha,dt,gmres_tol,rigid,vesicle)
%function [A,jaco,ivs,listnv] = preprocessRigid(o,vgrad,ids,vols,G,vblock,N,nv,Nrd,nvrd,alpha,dt,gmres_tol,rigid,rhsRig,rhsRigChange,Xrig)
  nivs = max(ids);
  A = zeros(nivs,nivs);
  ivs = zeros(nivs,1);
  vtoiv = zeros(nv+nvrd,nivs);
  if o.withRigid
    %rhsRig = reshape(rhsRig,[],nvrd);
    %rhsRigChange = reshape(rhsRigChange,[],nvrd);
  end
  jI = [];
  jJ = [];
  jV = [];
  listnv = [];
  Nmax = max([N;Nrd]);
  for i = 1:2*N*nv+2*Nrd*nvrd
    if(ids(i)~=0)
      k = ceil(i/(2*Nmax));
      listnv = [listnv;k];
      jI = [ids(i);jI];
      jJ = [i;jJ];
      jV = [vgrad(i);jV];
      vtoiv(k,ids(i)) = 1;
      ivs(ids(i)) = vols(i);
    end
  end
  listnv = unique(listnv);
  ivs(ivs>-1e-12) = -1e-12;
  jaco = sparse(jI,jJ,jV,nivs,2*N*nv+2*Nrd*nvrd);
  jacoSmooth = jaco*0;

  
  for i = 1:nivs
    S = find(vtoiv(:,i)~=0);
    for j = 1:numel(S)
      k = S(j);
      f = jaco(i,1+(k-1)*2*Nmax:2*Nmax+(k-1)*2*Nmax)';
      %add smooth force
      %f = o.f_smooth(full(f),N,1);
      jacoSmooth(i,1+(k-1)*2*Nmax:2*Nmax+(k-1)*2*Nmax) = f';
      if k>nv
        iRig = k-nv;
        cm = rigid.center(:,iRig);
        X = rigid.X(:,iRig);
        [stokeslet,rotlet] = o.getRS(rigid.X(:,iRig),rigid.sa(:,iRig),rigid.center(:,iRig),f);
        fColRig = stokeslet*2*pi;
        torqueColRig = rotlet*2*pi; 
        rhsUp = [zeros(2*Nrd,1);fColRig(:);torqueColRig(:)];
        %[Xn,iflagTemp,R,I,resvec] = gmres(@(X) o.rigidmatk(X,rigid,iRig),rhsRig(:,iRig)+rhsRigChange(:,iRig)+rhsUp,[],1e-2*o.gmresTol,2*Nrd+3);
        %Xrigk = o.extractRHSrigidk(Xn,rigid,iRig);
        %b = (Xrigk - Xrig(:,iRig))/dt;

        [Xn,iflagTemp,R,I,resvec] = gmres(@(X) o.rigidmatk(X,rigid,iRig),rhsUp,[],o.gmresTol,2*Nrd+3,@(x)o.preconditionerBDRigidSingle(x,iRig));
	%myXn = o.bdiagRigid.U(:,:,iRig)\(o.bdiagRigid.L(:,:,iRig)\rhsUp);
        b = repmat(Xn(end-2:end-1)',Nrd,1) + Xn(end)*[X(Nrd+1:end)-cm(2),-(X(1:Nrd)-cm(1))];
        b = b(:);
      else
	rhsUpdatek = [1/alpha(k)*G(:,:,k)*f;zeros(N,1)];
        vesiclek = capsules(vesicle.X(:,k),vesicle.sig(:,k),vesicle.u(:,k),vesicle.kappa,vesicle.viscCont(k),o.antiAlias); 
        [Xn,iflagTemp,R,I,resvec] = gmres(@(X) o.vesmatSingle(X,vesiclek,k),...
            rhsUpdatek,[],o.gmresTol*1e-2,o.gmresMaxIter,...
            @(x)o.preconditionerBDVesSingle(x,k));
	b = Xn(1:2*N);

        %b = o.bdiagVes.U(:,:,k)\(o.bdiagVes.L(:,:,k)\[1/alpha(k)*G(:,:,k)*f;zeros(N,1)]);
        %b = b(1:2*N);
        %b = o.filterHighFreq(b);
      end

      SS = find(vtoiv(k,:)~=0);
      for l = 1:numel(SS)
        A(SS(l),i) = A(SS(l),i) + dot(jaco(SS(l),1+(k-1)*2*Nmax:2*Nmax+(k-1)*2*Nmax),b);
      end
    end
  end
end % preprocessRigid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,jaco,ivs,listnv] = preprocess(o,vgrad,ids,vols,G,vblock,N,nv,alpha,dt,gmres_tol)
  nivs = max(ids);
  A = zeros(nivs,nivs);
  ivs = zeros(nivs,1);
  vtoiv = zeros(nv,nivs);
  jI = [];
  jJ = [];
  jV = [];
  listnv = [];
  for i = 1:2*N*nv
    if(ids(i)~=0)
      k = ceil(i/(2*N));
      listnv = [listnv;k];
      jI = [ids(i);jI];
      jJ = [i;jJ];
      jV = [vgrad(i);jV];
      vtoiv(k,ids(i)) = 1;
      ivs(ids(i)) = vols(i);
    end
  end
  listnv = unique(listnv);
  ivs(ivs>-1e-12) = -1e-12;
  jaco = sparse(jI,jJ,jV,nivs,2*N*nv);

  
  for i = 1:nivs
    S = find(vtoiv(:,i)~=0);
    for j = 1:numel(S)
      k = S(j);
      f = jaco(i,1+(k-1)*2*N:2*N+(k-1)*2*N)';
      %rhs = [dt/alpha*G(:,:,k)*f;zeros(N,1)];
      %[b,iflagw,Rw,Iw,resvecw] = gmres(Av(:,:,k),rhs,[],gmres_tol,3*N,vblock(:,:,k));

      %add smooth force
      %f = o.f_smooth(full(f),N,1);
   
      b = o.bdiagVes.U(:,:,k)\(o.bdiagVes.L(:,:,k)\[1/alpha(k)*G(:,:,k)*f;zeros(N,1)]);
      %b=vblock(:,:,k)*[dt/alpha*G(:,:,k)*f;zeros(N,1)];
      b = b(1:2*N);


      SS = find(vtoiv(k,:)~=0);
      for l = 1:numel(SS)
        A(SS(l),i) = A(SS(l),i) + dot(jaco(SS(l),1+(k-1)*2*N:2*N+(k-1)*2*N),b);
      end
    end
  end
end % preprocess

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vgrad] = adjustNormal(o,vgrad,N,nv,vesicle,edgelength,rig,colCount)
% use the normal of Xn or X?
vgrad = reshape(vgrad,2*N,nv);

for nvi = 1:nv
  for vgradi = 1:N
    if(vgrad(vgradi,nvi)==0 && vgrad(N+vgradi,nvi)==0)
      continue;
    end
    
    n_tmp = zeros(2,1);
    n_tmp(1) = vgrad(vgradi,nvi);
    n_tmp(2) = vgrad(N+vgradi,nvi);

    % in rare case if space-time gradient is too small, 
    % collision resolving may get stuck
    % use surface normal instead
    %if(norm(n_tmp) < edgelength(nvi)*0.1 || colCount > 20)
    if(colCount > 20)
      bnnorm = edgelength(nvi);
      if rig
        vgrad(vgradi,nvi) = vesicle.normal(vgradi,nvi)*bnnorm;
        vgrad(N+vgradi,nvi) = vesicle.normal(vgradi+N,nvi)*bnnorm;
      else
        vgrad(vgradi,nvi) = -vesicle.normal(vgradi,nvi)*bnnorm;
        vgrad(N+vgradi,nvi) = -vesicle.normal(vgradi+N,nvi)*bnnorm;
      end
    end
  end
end

vgrad = vgrad(:);
end % adjustNormal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vgrad] = f_smooth(o,vgrad,N,nv)
  vgrad = reshape(vgrad,2*N,nv);
  M = ceil(N/4);
  gw = gausswin(N,10);
  gw = gw/sum(gw);
  for nvi = 1:nv
    fc_x = vgrad(1:N,nvi);
    fc_tmp = ifft(fft(fc_x).*fft(gw));
%    Idx = ifftshift(-N/2:N/2+1);
%    Idx = abs(Idx) >= N/4;
%    coeffx = fft(fc_x).*fft(gw);
%    norm(coeffx(Idx));
    fc_x = [fc_tmp(N/2:N);fc_tmp(1:N/2-1)];
%    coeffx = fft(fc_x);
%    norm(coeffx(Idx));

    %coef = fft(fc_x);
    %coef(M+1:N/2+1) = 0;
    %coef(N/2+2:N-M+1) = 0;
    %fc_x = ifft(coef);

    fc_y = vgrad(N+1:2*N,nvi);
    fc_tmp = ifft(fft(fc_y).*fft(gw));
%    coeffx = fft(fc_y).*fft(gw);
%    norm(coeffx(Idx));
    fc_y = [fc_tmp(N/2:N);fc_tmp(1:N/2-1)];
    %coef = fft(fc_y);
    %coef(M+1:N/2+1) = 0;
    %coef(N/2+2:N-M+1) = 0;
    %fc_y = ifft(coef);

    vgrad(1:N,nvi) = fc_x;
    vgrad(N+1:2*N,nvi) = fc_y;
  end
end % f_smooth

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fc,lambda] = getColForce(o,A,b,x0,jaco)
tol_rel = 0.0;
tol_abs = 0.0;
max_iter = 50;
[lambda, err, iter, flag, conv, msg] = fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'perturbation', false);
fc = jaco'*lambda;
end % getColForce

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rhsUpdate] = getColRHS(o,G,alpha,dt,N,nv,fc)
rhsUpdate = zeros(3*N*nv,1);
for i = 1:nv
  rhsUpdate((i-1)*3*N+1:i*3*N) = [dt/alpha(i)*G(:,:,i)*fc((i-1)*2*N+1:i*2*N);zeros(N,1)];
end
end % getColRHS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rhs] = getColRHSrigid(o,rigid,fc)
nv = rigid.nv;
N = rigid.N;
fc = reshape(fc,2*N,nv);
u0 = zeros(2*N,nv);
fColRig = zeros(2,nv);
torqueColRig = zeros(1,nv);
for i = 1:nv
   [stokeslet, rotlet] = o.getRS(rigid.X(:,i),rigid.sa(:,i),rigid.center(:,i),fc(:,i));
   fColRig(:,i) = stokeslet*2*pi;
   torqueColRig(i) = rotlet*2*pi;
end
rhs = [u0;fColRig;torqueColRig];
rhs = rhs(:);
end % getColRHSrigid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vesVel,divVesVel,wallVel,residual] = ...
      computeVelocitiessdccol(o,FC,vesicle,Galpert,Gbarnett,D,walls,...
          etaProv,RSprov,NearV2V,NearV2W,NearW2V,NearW2W,uProv)
% [vesVel,divVesVel,wallVel,residual] = ...
%      computeVelocities(vesicle,Galpert,Gbarnett,D,walls,...
%       etaProv,RSprov,NearV2V,NearV2W,NearW2V,NearW2W);
% computes the velocity induced by the provisional vesicle position,
% tension, and density function on the vesicles and the walls.  Also
% returns the vesicle divergence of the velocity field and the residual
% of the picard integral formulation of the vesicle velocity.  These
% quantities are needed to form the modificiations of the right-hand
% side when doing SDC updates

N = vesicle.N;
nv = vesicle.nv;
if o.confined
  Nbd = walls.N;
  nvbd = walls.nv;
else
  Nbd = 0;
  nvbd = 0;
end

vesVel = zeros(2*N,nv,abs(o.orderGL));
% velocity on the vesicles due to the provisional solution 
divVesVel = zeros(N,nv,abs(o.orderGL));
% vesicle divergence of the vesicles due to the provisional
% solution
wallVel = zeros(2*Nbd,nvbd,abs(o.orderGL));
% velocity on the solid walls due to the provisional solution
residual = zeros(2*N,nv,abs(o.orderGL));
% residual of the Picard integral coming from the time
% derivative term of the vesicle position

vesves = o.vesves;
o.vesves = 'implicit';
% need to use implicit so that other vesicles are used to compute
% the integrand z
order = o.order;
o.order = 1;
[o.Xcoeff,o.rhsCoeff,o.beta] = o.getCoeff(o.order);
% need to save the time stepping order
dt = o.dt;
% need to save current time step
o.dt = 1;
% avoid introducing numerical error by multiplying and dividing
% by a potentially small number

if ~o.confined
  z = zeros(3*N*nv,1);
else
  z = zeros(3*N*nv + 2*Nbd*nvbd + 3*(nvbd-1),1);
end
op = o.op;
for n = 1:abs(o.orderGL)
  % loop over the configurations
  o.Galpert = Galpert(:,:,:,n);
%  o.Gbarnett = Gbarnett(:,:,:,n);
  if any(vesicle(n).viscCont ~= 1)
    o.D = D(:,:,:,n);
  end
  % Already computed the single-layer potential when forming the
  % provisional solution Store the single- and double-layer potentials
  % in the object o
  o.NearV2V = NearV2V{n};
  o.NearV2W = NearV2W{n};
  o.NearW2V = NearW2V{n};
  o.NearW2W = NearW2W{n};
  % retrieve near-singular integration strucutre at Gauss-Lobatto time
  % step n

  for k = 1:nv
    z(3*(k-1)*N+1:3*k*N) = [vesicle(n).X(:,k);vesicle(n).sig(:,k)];
  end
  % put positions and tension into z the way that TimeMatVec requires

  for k = 1:nvbd
    istart = 3*nv*N+2*(k-1)*Nbd+1;
    iend = istart + 2*Nbd - 1;
    z(istart:iend) = etaProv(:,k,n);
  end
  % put in density function on solid walls into z the way that
  % TimeMatVec requires them

  for k = 2:nvbd
    istart = 3*nv*N+2*nvbd*Nbd+3*(k-2)+1;
    iend = istart + 2;
    z(istart:iend) = RSprov(:,k,n);
  end
  % put the rotlet and stokeslet coefficients into z

  viscCont = vesicle(n).viscCont;
  vesicle(n).viscCont = ones(1,nv);
  % don't want the double-layer contribution when computing the velocity
  % of the vesicle
  z = o.TimeMatVec(z,vesicle(n),walls,1);
  % use TimeMatVec to find velocity on vesicles and solid
  % walls due to the provisional solution
  vesicle(n).viscCont = viscCont;
  % set the viscosity contrast back to its original value

  rhs = zeros(2*N,nv);
  for k = 1:nv
    istart = 3*(k-1)*N+1;
    iend = istart + 2*N - 1;
    rhs(:,k) = -z(istart:iend);
  end
  if ~o.confined
    rhs = rhs + o.farField(vesicle(n).X);
  end
  
  %velocity by collsion force
  for k = 1:nv
    rhs(:,k) = rhs(:,k) + o.Galpert(:,:,k)*FC(:,k,n);
  end

  vesVel(:,:,n) = vesicle(n).X + rhs;
  % form the velocity on the vesicle due to the current provisional
  % solution

  if any(vesicle(n).viscCont ~= 1)
    if o.profile
      tic
    end
    if strcmp(vesves,'explicit')
      
      jump = 1/2*(1-vesicle(n).viscCont);
      DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle(n),o.D,X);
      if ~o.fmmDLP
        kernel = @op.exactStokesDL;
      else
        kernel = @op.exactStokesDLfmm;
      end
      kernelDirect = kernel;
      den = uProv(:,:,n);
      Fdlp = op.nearSingInt(vesicle(n),den,DLP,DLP,...
          o.NearV2V,kernel,kernelDirect,vesicle(n),true);
      vesVel(:,:,n) = vesVel(:,:,n) + Fdlp;
      vesVel(:,:,n) = o.solveIminusDexp(vesVel(:,:,n),vesicle(n));
      
      %vesVel(:,:,n) = o.solveIminusD(vesVel(:,:,n),vesicle(n));
    else
      vesVel(:,:,n) = o.solveIminusD(vesVel(:,:,n),vesicle(n));
    end
    if o.profile
      fprintf('Solve system alpha*I - DLP          %5.1e\n',toc);
    end
  end 
  % need to apply inv(alpha*I - DLP) if there is a viscosity contrast to
  % obtain the velocity of the vesicles.

  divVesVel(:,:,n) = vesicle(n).surfaceDiv(vesVel(:,:,n));
  % surface divergence of the velocity of the vesicle due
  % to the provisional solution

  for k = 1:nvbd
    istart = 3*nv*N + (k-1)*2*Nbd + 1;
    iend = istart + 2*Nbd - 1;
    wallVel(:,k,n) = z(istart:iend);
  end

  if any(vesicle(n).viscCont ~= 1) && o.confined
    jump = 1/2*(1-vesicle(n).viscCont);
    DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle(n),o.D,X);
    if ~o.fmmDLP
      kernel = @op.exactStokesDL;
    else
      kernel = @op.exactStokesDLfmm;
    end

    kernelDirect = kernel;
    den = vesVel(:,:,n);
    FDLPwall = op.nearSingInt(vesicle(n),den,DLP,DLP,...
        o.NearV2W,kernel,kernelDirect,walls,false);
    % Evaulate the velocity due to the viscosity contrast
    % on the walls WITH near-singulation integration
    wallVel(:,:,n) = wallVel(:,:,n) + FDLPwall;
  end
  % velocity on the solid walls due to the provisional solution
end

o.dt = dt;
o.order = order;
[o.Xcoeff,o.rhsCoeff,o.beta] = o.getCoeff(o.order);
o.vesves = vesves;
% change back to original time step and vesicle-vesicle interaction

IvesVel = o.lobattoInt(vesVel);
% integrate the vesicles velocity using quadrature rules 
% that are exact for polynomials defined at the 
% Gauss-Lobatto points
for n = 1:abs(o.orderGL)
  residual(:,:,n) = vesicle(1).X - vesicle(n).X + ...
      o.dt/2 * IvesVel(:,:,n);
end
% compute residual by adding the initial vesicle configuartion and
% subtracting the current vesicle configuartion

end % computeVelocitiessdccol

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x] = plotDisplacement(o,X0,X1,XchangeTot,Xup,vgrad)
oc = curve;
[x0,y0] = oc.getXY(X0);
[x1,y1] = oc.getXY(X1+XchangeTot);
[x2,y2] = oc.getXY(X1+(XchangeTot+Xup));
hold on;
plot([x0;x0(1,:)],[y0;y0(1,:)],'r-x');
plot([x1;x1(1,:)],[y1;y1(1,:)],'g-x');
plot([x2;x2(1,:)],[y2;y2(1,:)],'b-x');

end % plotDisplacement

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = IminusDexp(o,X,vesicle)

op = o.op;
N = vesicle.N;
nv = vesicle.nv;
alpha = 0.5*(1+vesicle.viscCont);

Xm = zeros(2*N,nv);
for k = 1:nv
  Xm(:,k) = X((k-1)*2*N+1:k*2*N);
end

val = Xm*diag(alpha);
% "jump" term since we are computing alpha * I - DLP

if vesicle.N > 128
  val = val - op.exactStokesDLdiag(vesicle,o.D,Xm);
else
  val = val - op.exactStokesDLdiag(vesicle,[],Xm);
end
% self-interaction term

val = val(:);
end % IminusDexp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zOut = solveIminusDexp(o,zIn,vesicle)

warning off
%tGMRES = tic;
%[zIn,flag,relres,iter] = gmres(@(X) o.IminusD(X,vesicle),zIn(:),...
%    [],1e-2*o.gmresTol,min(2*vesicle.N*vesicle.nv,o.gmresMaxIter));
[zIn,flag,relres,iter] = gmres(@(X) o.IminusDexp(X,vesicle),zIn(:),...
    [],1e-2*o.gmresTol,min(2*vesicle.N*vesicle.nv,o.gmresMaxIter),...
    @(z) o.precoIminusD(z,vesicle));
warning on
%fprintf('GMRES time is %4.2e\n',toc(tGMRES))
% solve with block-diagonal preconditioned GMRES.  Integral equation is
% of the form identity + compact.  Need a bit more accuracy in this
% gmres solve as it is an inner iteration within an outer GMRES
% iteration.  This hides the fact that the GMRES solver is not linear

zOut = zeros(2*vesicle.N,vesicle.nv);
for k = 1:vesicle.nv
  zOut(:,k) = zIn((k-1)*2*vesicle.N+1:2*k*vesicle.N);
end
% Sort output in appropriate format

end % solveIminusDexp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fgrav = gravityForce(o,Xm,vesicle,gCont)
N = size(Xm,1)/2;
nv = size(Xm,2);
fgrav = zeros(2*N,nv);
g = [0;-1];
for i = 1:nv
  for j = 1:N
    Xtmp = [Xm(j,i);Xm(N+j,i)];
    Xtmp = Xtmp - vesicle.center(:,i);
    gsize = gCont*sum(g.*Xtmp);
    fgrav(j,i) = gsize*vesicle.normal(j,i); 
    fgrav(N+j,i) = gsize*vesicle.normal(N+j,i);
  end
end
end % gravityForce

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = preconditionerBDVes(o,z)

if (strcmp(o.solver,'method1') || strcmp(o.solver,'method2'))
  nv = size(o.bdiagVes.L,3); % number of vesicles
  N = size(o.bdiagVes.L,1)/3; % number of points
elseif strcmp(o.solver,'method3')
  nv = size(o.bdiagVes.schur,3); % number of vesicles
  N = size(o.bdiagVes.schur,1)/2; % number of vesicles
elseif strcmp(o.solver,'method4')
  nv = size(o.bdiagVes.schur,3); % number of vesicles
  N = size(o.bdiagVes.schur,1); % number of vesicles
end

zves = z(1:3*N*nv);
% extract the position and tension part.  Solid walls is
% handled in the next section of this routine
valVes = zeros(3*N*nv,1);

if (strcmp(o.solver,'method1') || strcmp(o.solver,'method2'))
  for k=1:nv
    valVes((k-1)*3*N+1:3*k*N) = o.bdiagVes.U(:,:,k)\...
      (o.bdiagVes.L(:,:,k)\zves((k-1)*3*N+1:3*k*N));
  end % k
  % precondition with the block diagonal preconditioner for the
  % vesicle position and tension
elseif strcmp(o.solver,'method3')
  for k = 1:nv
    istart = (k-1)*3*N+1;
    iend = istart + 2*N - 1;
    bx = zves(istart:iend); 
    istart = iend + 1;
    iend = istart + N - 1;
    bsig = zves(istart:iend);
    % seperate position and tension terms

    rhs = bx + o.dt*o.bdiagVes.GT(:,:,k)*o.bdiagVes.DGT(:,:,k)*bsig;
    pos = o.bdiagVes.schur(:,:,k)*rhs;
    rhs = bsig + o.bdiagVes.DGB(:,:,k)*pos;
    sig = o.bdiagVes.DGT(:,:,k)*rhs;
    valVes((k-1)*3*N+1:3*k*N) = [pos;sig];
    % use schur decomposition to operate preconditioner
  end % k
elseif strcmp(o.solver,'method4')
  for k = 1:nv
    istart = (k-1)*3*N+1;
    iend = istart + 2*N - 1;
    bx = zves(istart:iend); 
    istart = iend + 1;
    iend = istart + N - 1;
    bsig = zves(istart:iend);
    % seperate position and tension terms

    rhs = bsig + o.bdiagVes.DGB(:,:,k)*o.bdiagVes.IpBen(:,:,k)*bx;
    sig = o.bdiagVes.schur(:,:,k)*rhs;
    rhs = bx + o.dt*o.bdiagVes.GT(:,:,k)*sig;
    pos = o.bdiagVes.IpBen(:,:,k)*rhs;
    valVes((k-1)*3*N+1:3*k*N) = [pos;sig];
    % use schur decomposition to operate preconditioner
  end % k
end % o.solver

val = valVes(:);
% stack the two componenets of the preconditioner
end % preconditionerBDVes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = preconditionerBDVesSingle(o,z,i)

if (strcmp(o.solver,'method1') || strcmp(o.solver,'method2'))
  %nv = size(o.bdiagVes.L,3); % number of vesicles
  N = size(o.bdiagVes.L,1)/3; % number of points
elseif strcmp(o.solver,'method3')
  %nv = size(o.bdiagVes.schur,3); % number of vesicles
  N = size(o.bdiagVes.schur,1)/2; % number of vesicles
elseif strcmp(o.solver,'method4')
  %nv = size(o.bdiagVes.schur,3); % number of vesicles
  N = size(o.bdiagVes.schur,1); % number of vesicles
end
nv = 1;

zves = z(1:3*N*nv);
% extract the position and tension part.  Solid walls is
% handled in the next section of this routine
valVes = zeros(3*N*nv,1);

if (strcmp(o.solver,'method1') || strcmp(o.solver,'method2'))
  for k=i:i
    valVes = o.bdiagVes.U(:,:,k)\...
      (o.bdiagVes.L(:,:,k)\zves(:));
  end % k
  % precondition with the block diagonal preconditioner for the
  % vesicle position and tension
elseif strcmp(o.solver,'method3')
elseif strcmp(o.solver,'method4')
end % o.solver

val = valVes(:);
% stack the two componenets of the preconditioner
end % preconditionerBDVesSingle


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = preconditionerBDWall(o,z)
zwalls = z(1:end);
% part of z from solid walls
valWalls = o.bdiagWall * zwalls;
% this matrix is well-coditioned since it is in the form I + DLP
val = valWalls(:);
% stack the two componenets of the preconditioner
end % preconditionerBDWall

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = preconditionerBDRigid(o,z)
nv = size(o.bdiagRigid.L,3); % number of vesicles
N = (size(o.bdiagRigid.L,1)-3)/2; % number of points
val = zeros(2*N+3,nv);
for i = 1:nv
  val(:,i) = (o.bdiagRigid.U(:,:,i))\(o.bdiagRigid.L(:,:,i)\z(1+(i-1)*(2*N+3):i*(2*N+3)));
end
val = val(:);
end % preconditionerBDRigid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = preconditionerBDRigidSingle(o,z,i)
nv = 1; % number of vesicles
N = (size(o.bdiagRigid.L,1)-3)/2; % number of points
val = zeros(2*N+3,nv);
val = (o.bdiagRigid.U(:,:,i))\(o.bdiagRigid.L(:,:,i)\z);
val = val(:);
end % preconditionerBDRigidSingle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extractRHSVes
function [X,sigma,u] = extractRHSVes(o,Xn,N,nv,Xo)
X = zeros(2*N,nv);
sigma = zeros(N,nv);
% allocate space for positions, tension, and density function

for k=1:nv
  X(:,k) = Xn((3*k-3)*N+1:(3*k-1)*N);
  sigma(:,k) = Xn((3*k-1)*N+1:3*k*N);
end
% unstack the positions and tensions
%X = o.filterHighFreq(X);
%sigma = o.filterHighFreqTen(sigma);

u = (o.beta*X - Xo)/o.dt;
% Compute the velocity using the differencing stencil
end % extractRHSVes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extractRHSWall
function [eta,RS] = extractRHSWall(o,Xn,Nbd,nvbd)
eta = zeros(2*Nbd,nvbd);
RS = zeros(3,nvbd);
% allocate space for positions, tension, and density function


for k = 1:nvbd
  eta(:,k) = Xn((k-1)*2*Nbd+1:2*k*Nbd);
end
% unstack the density function
otlets = Xn(2*nvbd*Nbd+1:end);
for k = 2:nvbd
  istart = (k-2)*3+1;
  iend = 3*(k-1);
  RS(:,k) = otlets(istart:iend);
end
% unstack the rotlets and stokeslets
end % extractRHSWall


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,sigma,u,eta,RS,Xrig,etaRig,fc_tot_store] = resolveCollisionGSrigid(o,X0,X1,walls,rigid,PROPCOL)

viscCont = PROPCOL.viscCont;
vesicle = PROPCOL.vesicle;
N = PROPCOL.N;
nv = PROPCOL.nv;
Nbd = PROPCOL.Nbd;
nvbd = PROPCOL.nvbd;
Nrd = PROPCOL.Nrd;
nvrd = PROPCOL.nvrd;
Nmax = max([N;Nrd]);
X = PROPCOL.X;
sigma = PROPCOL.sigma;
u = PROPCOL.u;
eta = PROPCOL.eta;
RS = PROPCOL.RS;
Xrig = PROPCOL.Xrig;
etaRig = PROPCOL.etaRig;
om = PROPCOL.om;
if o.withRigid
  uRig = PROPCOL.uRig;
  omgRig = PROPCOL.omgRig;
  rhsRig = PROPCOL.rhsRig;
  rhsRigChange = rhsRig*0;
end
oc = curve;
cellSize = 0;
if nv
  [ra,area,length] = oc.geomProp(X0(:,1:nv));
  edgelength = length/N;
  cellSize = max(cellSize,max(edgelength));
end
if o.withRigid
  [ra,area,length] = oc.geomProp(X0(:,nv+1:end));
  edgelengthRig = length/Nrd;
  cellSize = max(cellSize,max(edgelengthRig));
end
if o.confined
  [ra,area,length] = oc.geomProp(walls.X);
  walllength = max(length/Nbd);
  wallupSamp = ceil(walllength/cellSize);
else
  wallupSamp = 1;
end
alpha = (1+viscCont)/2;
gmres_tol = 1e-12;
upSampleFactor = 1;
nexten = 0;
c_tol = 1e-12;
minSep = o.minSep;
maxIter = 1000;

[Ns,totalnv,Xstart,Xend,totalPts] = o.preColCheck(X0,X1,walls,wallupSamp);
[vgrad, iv, ids, vols] = getCollision(Ns, totalnv, Xstart, Xend, minSep, maxIter, totalPts, c_tol, nv+nvrd, nvbd, Nmax*upSampleFactor, Nbd*wallupSamp, nexten, max(cellSize,minSep));
om.writeMessage(['ivolume: ' num2str(iv)],'%s\n');

fc_tot = zeros(2*N*nv+2*Nrd*nvrd,1);
fc_tot_store = zeros(2*N*nv+2*Nrd*nvrd,1);
XchangeTot = X*0;
sigmaChangeTot = sigma*0;
uChangeTot = u*0;
etaChangeTot = eta*0;
RSchangeTot = RS*0;
X1tmp = X1;

colCount = 0;
if(iv<0 && minSep==0)
om = PROPCOL.om;
om.writeMessage('collision, exit matlab.','%s\n');
exit;
end
while(iv<0)
  colCount = colCount + 1;
  vgrad = vgrad(1:2*N*nv+2*Nrd*nvrd);
  ids_tmp = ids;
  ids = ids(1:2*N*nv+2*Nrd*nvrd);
  vols = vols(1:2*N*nv+2*Nrd*nvrd);
  if nv
    vgrad(1:2*N*nv) = o.adjustNormal(vgrad(1:2*N*nv),N,nv,vesicle,edgelength,0,colCount);
  end
  if o.withRigid
    vgrad(1+2*N*nv:end) = o.adjustNormal(vgrad(1+2*N*nv:end),Nrd,nvrd,rigid,edgelengthRig,1,colCount);
  end

  [A,jaco,ivs,listnv,jacoSmooth,vtoiv] = o.preprocessRigid(vgrad,ids,vols,o.Galpert,o.bdiagVes,N,nv,Nrd,nvrd,alpha,o.dt,gmres_tol,rigid,vesicle);
  [fc, lambda] = o.getColForce(A,ivs/o.dt,ivs*0,jacoSmooth);
  fc_tot = fc_tot + fc;
  
  for k = 1:numel(listnv)
    i = listnv(k);
    if i > nv
      iRig = i - nv;
      fc_toti = fc_tot(2*N*nv+(iRig-1)*(2*Nrd)+1:2*N*nv+iRig*(2*Nrd));
      [stokeslet,rotlet] = o.getRS(rigid.X(:,iRig),rigid.sa(:,iRig),rigid.center(:,iRig),fc_toti);
      fColRig = stokeslet*2*pi;
      torqueColRig = rotlet*2*pi; 
      rhsRigChangei = [zeros(2*Nrd,1);fColRig(:);torqueColRig(:)];
      [Xn,iflagTemp,R,I,resvec] = gmres(@(X) o.rigidmatk(X,rigid,iRig),rhsRig((iRig-1)*(2*Nrd+3)+1:iRig*(2*Nrd+3))+rhsRigChangei,[],o.gmresTol,2*Nrd+3,@(x)o.preconditionerBDRigidSingle(x,iRig));
      [Xrigk, etaRigk, uRigk, omgRigk] = o.extractRHSrigidk(Xn,rigid,iRig);
      Xrig(:,iRig) = Xrigk;
      X1tmp(:,nv+iRig) = Xrig(:,iRig);
      etaRig(:,iRig) = etaRigk;
      etaRigknorm = norm(etaRigk)
      etaRigknormwithoutcontact = norm(PROPCOL.etaRig(:,iRig))
    else
      rhsUpdatei = [o.dt/alpha(i)*o.Galpert(:,:,i)*fc((i-1)*2*N+1:i*2*N);zeros(N,1)];
      vesiclei = capsules(vesicle.X(:,i),vesicle.sig(:,i),vesicle.u(:,i),vesicle.kappa,vesicle.viscCont(i),o.antiAlias); 
      [Xn,iflagTemp,R,I,resvec] = gmres(@(X) o.vesmatSingle(X,vesiclei,i),...
          rhsUpdatei,[],o.gmresTol*1e-2,o.gmresMaxIter,...
          @(x)o.preconditionerBDVesSingle(x,i));
      [Xup,sigmaup,uup] = o.extractRHSVes(Xn,N,1,X0(:,i));
      XchangeTot(:,i) = XchangeTot(:,i) + Xup;
      sigmaChangeTot(:,i) = sigmaChangeTot(:,i) + sigmaup;
      X1tmp(:,i) = X1(:,i) + XchangeTot(:,i);
    end
  end
  
  % ignore part of the contact force on vesicle if the contact force is vesicle-boundary contact
  bd_ids = [];
  for k = (2*N*nv+2*Nrd*nvrd+1):size(ids_tmp,1)
    if ids_tmp(k) ~= 0
      lambda(ids_tmp(k)) = 0;
      if o.withRigid
          bd_ids = [bd_ids;ids_tmp(k)];
      end
    end
  end
  
  % ignore part of the contact force on rigid if the contact force is rigid-boundary contact
  if (o.withRigid) && (size(bd_ids,1) > 0)
    bd_ids = unique(bd_ids);
    for k = 1:size(bd_ids);
      pt_ind = find(ids_tmp(2*N*nv+1:2*N*nv+2*Nrd*nvrd) == bd_ids(k));
      if size(pt_ind,1) > 0
        rig_ind = unique(ceil(pt_ind/(2*Nrd)));
        etaRig(:,rig_ind) = PROPCOL.etaRig(:,rig_ind);
        % since we ignore all the contact force on the rigid if the rigid contact with boundary
        % if this rigid also contact with some vesicle, this part of contact force is ignored 
        % too to make force balance
        for kk = 1:size(rig_ind)
            ves_rgd_ind = find(vtoiv(nv+rig_ind(kk),:) == 1)
            lambda(ves_rgd_ind) = 0;
        end
      end
    end
  end
  
  fc_tot_store = fc_tot_store + jacoSmooth'*lambda;
  
  [Ns,totalnv,Xstart,Xend,totalPts] = o.preColCheck(X0,X1tmp,walls,wallupSamp);
  [vgrad, iv, ids, vols] = getCollision(Ns, totalnv, Xstart, Xend, minSep, maxIter, totalPts, c_tol, nv+nvrd, nvbd, Nmax*upSampleFactor, Nbd*wallupSamp, nexten, max(cellSize,minSep));

  %{
  if (colCount > 50 && colCount < 2500)
    filename = 'debugtest';
    fid = fopen(filename,'a');
    dval = [PROPCOL.time(:);Ns(:);totalnv(:);Xstart(:);Xend(:);minSep(:);maxIter(:);totalPts(:);c_tol(:);nv+nvrd;nvbd;Nmax*upSampleFactor;Nbd*wallupSamp;nexten;max(cellSize,minSep)];
    fwrite(fid,dval,'double');
    fclose(fid);
  end
  if colCount > 8000
    om.writeMessage('too many collision iterations','%s\n');
    %exit;
  end
  %}

  om.writeMessage(['ivolume: ' num2str(iv)],'%s\n');
end

X = X + XchangeTot;
sigma = sigma + sigmaChangeTot;
u = u + XchangeTot/o.dt;

if(o.confined && colCount)
  RIGID.rigid = rigid;
  RIGID.etaRig = etaRig;

  WALL.walls = walls;
  
  PROP.viscCont = viscCont;
  PROP.fc_tot = zeros(2*N,nv);
  PROP.om = om;
  SDCPROP = [];
  
  if o.SDCcorrect
    VES.X = XchangeTot; 
    VES.sig = sigmaChangeTot;
    VES.u = XchangeTot/o.dt;
  else
    VES.X = X; 
    VES.sig = sigma;
    VES.u = u;
  end

  [rhs] = o.sdcwallrhs(VES,WALL,PROP,SDCPROP,vesicle,RIGID);

  [Xn,iflagTemp,R,I,resvec] = gmres(@(X) o.wallmat(X,walls),...
        rhs,[],o.gmresTol,o.gmresMaxIter,...
        @o.preconditionerBDWall);
  if ~o.SDCcorrect
    [eta,RS] = o.extractRHSWall(Xn,Nbd,nvbd);
  else
    [etaChange,RSChange] = o.extractRHSWall(Xn,Nbd,nvbd);
    eta = eta + etaChange;
    RS = RS + RSChange;
  end
end
if nv
  fc_tot_store = reshape(fc_tot_store(1:2*N*nv),2*N,nv);
else
  fc_tot_store = [];
end
end % resolveCollisionGSrigid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [val] = rigidmat(o,Xn,rigid)
N = rigid.N;
nv = rigid.nv;
op = o.op;
DLP = o.rigidDLP;

val1 = zeros(2*N,nv); % velocity
val2 = zeros(2,nv); % total force
val3 = zeros(1,nv); % total torque

tmp = reshape(Xn,2*N+2+1,nv);
etaAll = tmp(1:2*N,:);
u = tmp(2*N+1:2*N+2,:);
w = tmp(2*N+3,:);

for i = 1:nv
  X = rigid.X(:,i);
  sa = rigid.sa(:,i);
  cm = rigid.center(:,i);
  eta = etaAll(:,i);

  stokeslet = zeros(2,1);
  rotlet = 0;
  stokeslet(1) = sum(eta(1:N).*sa)*2*pi/N; % force_X
  stokeslet(2) = sum(eta(N+1:end).*sa)*2*pi/N; % force_y
  rotlet = sum(...
    ((X(N+1:end)-cm(2)).*eta(1:N) - ...
    (X(1:N)-cm(1)).*eta(N+1:end)).*...
    sa)*2*pi/N; % torque
  
  val2(:,i) = stokeslet;
  val3(:,i) = rotlet;
  
  stokeslet = stokeslet/(2*pi);
  rotlet = rotlet/(2*pi);

  letsvel = o.RSlets(X,cm,stokeslet,rotlet);
  
  U = repmat(u(:,i)',N,1) + w(i)*[X(N+1:end)-cm(2),-(X(1:N)-cm(1))];
  val1(:,i) = 1/2*eta - DLP(:,:,i)*eta + U(:);% - letsvel;
end

val = [val1;val2;val3];
val = val(:);

end % rigidmat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [val] = rigidmatk(o,Xn,rigid,k)
N = rigid.N;
nv = numel(k);
op = o.op;
DLP = o.rigidDLP;

val1 = zeros(2*N,nv); % velocity
val2 = zeros(2,nv); % total force
val3 = zeros(1,nv); % total torque

tmp = reshape(Xn,2*N+2+1,nv);
etaAll = tmp(1:2*N,:);
u = tmp(2*N+1:2*N+2,:);
w = tmp(2*N+3,:);

for j=1:nv
  i = k(j);
  X = rigid.X(:,i);
  sa = rigid.sa(:,i);
  cm = rigid.center(:,i);
  eta = etaAll(:,j);

  stokeslet = zeros(2,1);
  rotlet = 0;
  stokeslet(1) = sum(eta(1:N).*sa)*2*pi/N; % force_X
  stokeslet(2) = sum(eta(N+1:end).*sa)*2*pi/N; % force_y
  rotlet = sum(...
    ((X(N+1:end)-cm(2)).*eta(1:N) - ...
    (X(1:N)-cm(1)).*eta(N+1:end)).*...
    sa)*2*pi/N; % torque
  
  val2(:,j) = stokeslet;
  val3(:,j) = rotlet;
  
  stokeslet = stokeslet/(2*pi);
  rotlet = rotlet/(2*pi);

  letsvel = o.RSlets(X,cm,stokeslet,rotlet);
  
  U = repmat(u(:,j)',N,1) + w(j)*[X(N+1:end)-cm(2),-(X(1:N)-cm(1))];
  val1(:,j) = 1/2*eta - DLP(:,:,i)*eta + U(:) - letsvel;
end

val = [val1;val2;val3];
val = val(:);

end % rigidmatk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sdcrigidrhs
function [rhs] = sdcrigidrhs(o,VES,WALL,PROP,SDCPROP,vesicle,RIGID)

if o.SDCcorrect
  deltaX = SDCPROP.deltaX;
  deltaSig = SDCPROP.deltaSig;
  wallVel = SDCPROP.wallVel;
  deltaEta = SDCPROP.deltaEta;
  deltaRS = SDCPROP.deltaRS;
  diffResidual = SDCPROP.diffResidual;
  vesVel = SDCPROP.vesVel;
  wallVel = SDCPROP.wallVel;
  sa = SDCPROP.sa;
  IK = SDCPROP.IK;
end

Xstore = VES.X;
sigStore = VES.sig;
uStore = VES.u;

walls = WALL.walls;
etaStore = WALL.eta;
RSstore = WALL.RS;

rigid = RIGID.rigid;
etaRigStore = RIGID.etaRig;

viscCont = PROP.viscCont;
kappa = PROP.kappa;

N = size(Xstore,1)/2; % Number of points per vesicle
nv = size(Xstore,2); % Number of vesicles
if o.confined
  Xwalls = walls.X; % discretization points of solid walls
else
  Xwalls = [];
end
Nbd = size(Xwalls,1)/2; % Number of points on the solid walls
nvbd = size(Xwalls,2); % number of solid wall components
if o.withRigid
  Xrigid = rigid.X;
else
  Xrigid = [];
end
Nrd = size(Xrigid,1)/2;
nvrd = size(Xrigid,2);

alpha = (1 + viscCont)/2; 
% constant that appears in front of time derivative in
% vesicle dynamical equations

Xm = zeros(2*N,nv);
sigmaM = zeros(N,nv);
uM = zeros(2*N,nv);
Xo = zeros(2*N,nv);
etaM = zeros(2*Nbd,nvbd);
RSm = zeros(3,nvbd);
etaRigM = zeros(2*Nrd,nvrd);
for k = 1:o.order
  Xm = Xm + Xstore(:,:,k)*o.Xcoeff(k);
  sigmaM = sigmaM + sigStore(:,:,k)*o.Xcoeff(k);
  uM = uM + uStore(:,:,k)*o.Xcoeff(k);
  if o.confined
    etaM = etaM + etaStore(:,:,k)*o.Xcoeff(k);
    RSm = RSm + RSstore(:,:,k)*o.Xcoeff(k);
  end
  if o.withRigid
    etaRigM = etaRigM + etaRigStore(:,:,k)*o.Xcoeff(k);
  end
  Xo = Xo + Xstore(:,:,k)*o.rhsCoeff(k);
end

op = o.op;

if o.collDete
  o.lapDLP = op.laplaceDLmatrix(vesicle);
  % build matrix that evaluates double-layer laplace potential for each
  % vesicle

  oc = curve;
  [icollisionVes,icollisionWall] = ...
    oc.collision(vesicle,walls,o.NearV2V,o.NearV2W,o.fmm,o.near);
  % Check for collisions 
  if icollisionVes
    fprintf('VESICLES HAVE CROSSED\n')
    pause
  end
  if icollisionWall
    fprintf('VESICLES HAVE CROSSED SOLID WALL\n')
    pause
  end
end
% check for collisions

if ~o.SDCcorrect
  rhs1 = zeros(2*Nrd,nvrd); %velocity
  rhs2 = zeros(2,nvrd); %force
  rhs3 = zeros(1,nvrd); %torque
else
  rhs1 = zeros(2*Nrd,nvrd); %velocity
  rhs2 = zeros(2,nvrd); %force
  rhs3 = zeros(1,nvrd); %torque
end
% Parts of rhs from previous solution.  The right-hand-side depends on
% whether we are doing an SDC correction or forming the provisional
% solution.
if nv
  % START TO COMPUTE RIGHT-HAND SIDE DUE TO VESICLE TRACTION JUMP
  if o.gravity
    % For now, treat gravity explicit same as backgroud flow
    fgrav = o.gravityForce(Xm,vesicle,o.gCont);
  else
    fgrav = zeros(2*N,nv);
  end
  if strcmp(o.vesves,'explicit')
    if o.profile
      tic
    end
    if ~o.SDCcorrect
      f = vesicle.tracJump(Xm,sigmaM) + PROP.fc_tot + fgrav;
    else
      f = vesicle.tracJump(deltaX,deltaSig);
    end
    if ~o.fmm
      kernel = @op.exactStokesSL;
    else
      kernel = @op.exactStokesSLfmm;
    end

    if ~o.near
      if o.withRigid
	[~,FSLPrigid] = kernel(vesicle,f,[],rigid.X,(1:nv));
	% Evaluate single-layer potential on solid walls due to all
	% vesicles
      else
	FSLPrigid = [];
      end
    else
      if o.nearStrat == 'cauchy'
      else
	SLP = @(X) op.exactStokesSLdiag(vesicle,o.Galpert,X);
	SLPtrap = SLP;
	kernelDirect = kernel;
      end
      if o.withRigid
	if o.nearStrat == 'cauchy'
	  FSLPrigid = op.nearSingStokesSLP(vesicle,f,rigid,...
	      o.Gbarnett,false,o.fmm);
	else
	  FSLPrigid = op.nearSingInt(vesicle,f,SLP,SLPtrap,...
	      o.NearV2R,kernel,kernelDirect,rigid,false);
	  % Evaluate the velocity on the walls due to the vesicles
	end
      else
	FSLPrigid = [];
      end
    end
    if o.profile
      fprintf('Build right-hand side V2V and V2W   %5.1e\n',toc);
    end
  end
  rhs1 = rhs1 + FSLPrigid;

  % START TO COMPUTE RIGHT-HAND SIDE DUE TO VISCOSITY CONTRAST
  if any(vesicle.viscCont ~= 1)
    if o.near
      jump = 1/2*(1-vesicle.viscCont);
      DLP = @(X) X*diag(jump) + op.exactStokesDLdiag(vesicle,o.D,X);
    end
    % Need to add jump to double-layer potential if using near-singular
    % integration so that we can compute boundary values for the
    % near-singular integration algorithm

    if ~o.fmmDLP
      kernel = @op.exactStokesDL;
    else
      kernel = @op.exactStokesDLfmm;
    end
    if strcmp(o.vesves,'explicit')
      if ~o.near
	if o.withRigid
	  [~,FDLPrigid] = kernel(vesicle,uM,rigid.X,(1:nv));
	  % Evaulate the velocity due to the viscosity contrast on the
	  % walls WITHOUT near-singulation integration
	else
	  FDLPrigid = [];
	end
      else
	kernelDirect = @op.exactStokesDL;

	if o.withRigid
	  FDLPrigid = op.nearSingInt(vesicle,uM,DLP,DLP,...
	      o.NearV2R,kernel,kernelDirect,rigid,false);
	  % Evaulate the velocity due to the viscosity contrast on the
	  % walls WITH near-singulation integration
	else
	  FDLPrigid = [];
	end
      end
    end
  else
    FDLPrigid = zeros(2*Nrd,nvrd);
    % If no viscosity contrast, there is no velocity induced due to a
    % viscosity contrast
  end

  rhs1 = rhs1 + FDLPrigid;
  % compute the double-layer potential due to all other vesicles from the
  % appropriate linear combination of previous time steps.  Depends on
  % time stepping order and vesicle-vesicle discretization
end

% START TO COMPUTE RIGHT-HAND SIDE DUE TO SOLID WALLS
if o.confined
  if strcmp(o.vesves,'explicit')
    if o.profile
      tic
    end
    if o.near
      jump = -1/2;
      DLP = @(X) jump*X + op.exactStokesDLdiag(walls,o.wallDLP,X);
    end
    % compute the matrix for doing evaluating the double-layer potential
    % on the solid walls required for near-singular integration

    if ~o.fmmDLP
      kernel = @op.exactStokesDL;
    else
      kernel = @op.exactStokesDLfmm;
    end
    if ~o.near
      if ~o.SDCcorrect
        charge = etaM;
      else
        charge = deltaEta;
      end
      [~,Fwall2Rigid] = kernel(walls,charge,...
          rigid.X,1:nvbd);
      % velocity field due to the walls evaluated on the vesicle
    else
      if ~o.SDCcorrect
        charge = etaM;
      else
        charge = deltaEta;
      end
      kernelDirect = kernel;
      Fwall2Rigid = op.nearSingInt(walls,charge,DLP,DLP,...
          o.NearW2R,kernel,kernelDirect,rigid,false);
    end

    for k = 2:nvbd
      if ~o.SDCcorrect
        stokeslet = RSm(1:2,k);
        rotlet = RSm(3,k);
      else
        stokeslet = deltaRS(1:2,k);
        rotlet = deltaRS(3,k);
      end
      Fwall2Rigid = Fwall2Rigid + ...
        o.RSlets(rigid.X,walls.center(:,k),stokeslet,rotlet);
    end
    if o.profile
      fprintf('Build right-hand side W2V           %5.1e\n',toc);
    end

  end
else
  Fwall2Rigid = zeros(2*Nrd,nvrd);
  if ~o.SDCcorrect
    rhs1 = rhs1 + o.farField(rigid.X);
    % Add in far-field condition (extensional, shear, etc.)
  end
end
rhs1 = rhs1 + Fwall2Rigid;
% right-hand side of the velocity evaluated on the solid walls
% END TO COMPUTE RIGHT-HAND SIDE DUE TO SOLID WALLS

if o.withRigid
  if strcmp(o.vesves,'explicit')
    if o.profile
      tic
    end
    if o.near
      jump = -1/2;
      DLP = @(X) jump*X + op.exactStokesDLdiag(rigid,o.rigidDLP,X);
    end
    % compute the matrix for doing evaluating the double-layer potential
    % on the solid walls required for near-singular integration

    if ~o.fmmDLP
      kernel = @op.exactStokesDL;
    else
      kernel = @op.exactStokesDLfmm;
    end
    if ~o.near
      if ~o.SDCcorrect
        charge = etaRigM;
      else
        charge = deltaRigEta;
      end
      FDLPrigid2rigid = kernel(rigid,charge,[]);
    else
      if ~o.SDCcorrect
        charge = etaRigM;
      else
        charge = deltaRigEta;
      end
      kernelDirect = kernel;
      FDLPrigid2rigid = op.nearSingInt(rigid,charge,DLP,DLP,...
          o.NearR2R,kernel,kernelDirect,rigid,true);
    end

    LetsRigid = zeros(2*Nrd,nvrd);
    for k = 1:nvrd
      if ~o.SDCcorrect
        [stokeslet, rotlet] = o.getRS(rigid.X(:,k),rigid.sa(:,k),rigid.center(:,k),etaRigM(:,k));
      else
        [stokeslet, rotlet] = o.getRS(rigid.X(:,k),rigid.sa(:,k),rigid.center(:,k),deltaRigEta(:,k));
      end
      LetsRigid_tmp =  o.RSlets(rigid.X,rigid.center(:,k),stokeslet,rotlet);
      % only inter-rigid particle interactions are added to rhs
      LetsRigid_tmp(:,k) = 0;
      %LetsRigid = LetsRigid + o.RSlets(rigid.X,rigid.center(:,k),stokeslet,rotlet);
      LetsRigid = LetsRigid + LetsRigid_tmp;
    end
    
    rhs1 = rhs1 + FDLPrigid2rigid + LetsRigid;

    if o.profile
      fprintf('Build right-hand side W2V           %5.1e\n',toc);
    end
  end
end

rhs = [rhs1;rhs2;rhs3];
rhs = [rhs(:)];
% Rotlet and Stokeslet equations

end % sdcrigidrhs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stokeslet, rotlet] = getRS(o,X,sa,cm,eta)
  N = size(X,1)/2;
  stokeslet = zeros(2,1);
  rotlet = 0;
  stokeslet(1) = sum(eta(1:N).*sa)*2*pi/N;
  stokeslet(2) = sum(eta(N+1:end).*sa)*2*pi/N;
  rotlet = sum(...
    ((X(N+1:end)-cm(2)).*eta(1:N) - ...
    (X(1:N)-cm(1)).*eta(N+1:end)).*...
    sa)*2*pi/N;
  
  stokeslet = stokeslet/(2*pi);
  rotlet = rotlet/(2*pi);
end %getRS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xrig, etaRig, uRig, omgRig] = extractRHSrigid(o,Xn,rigid)
N = rigid.N;
nv = rigid.nv;

tmp = reshape(Xn,2*N+2+1,nv);
etaRig = tmp(1:2*N,:);
uRig = tmp(2*N+1:2*N+2,:);
omgRig = tmp(2*N+3,:);
Xrig = zeros(2*N,nv);
for i = 1:nv
  X  = rigid.X(:,i);
  cm = rigid.center(:,i);
  ui = uRig(:,i);
  wi = omgRig(:,i);

  X0 = reshape([X(1:N)-cm(1);X(N+1:end)-cm(2)],N,2);
  cm = cm + ui*o.dt;
  X1 = [cm(1)+cos(wi*o.dt)*X0(:,1)+sin(wi*o.dt)*X0(:,2);cm(2)-sin(wi*o.dt)*X0(:,1)+cos(wi*o.dt)*X0(:,2)];
  Xrig(:,i) = X1(:);
end

end% extractRHSrigid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xrig, etaRig, uRig, omgRig] = extractRHSrigidk(o,Xn,rigid,k)
N = rigid.N;
nv = numel(k);

tmp = reshape(Xn,2*N+2+1,nv);
etaRig = tmp(1:2*N,:);
uRig = tmp(2*N+1:2*N+2,:);
omgRig = tmp(2*N+3,:);
Xrig = zeros(2*N,nv);
for j = 1:nv
  i = k(j);
  X  = rigid.X(:,i);
  cm = rigid.center(:,i);
  ui = uRig(:,j);
  wi = omgRig(:,j);

  X0 = reshape([X(1:N)-cm(1);X(N+1:end)-cm(2)],N,2);
  cm = cm + ui*o.dt;
  X1 = [cm(1)+cos(wi*o.dt)*X0(:,1)+sin(wi*o.dt)*X0(:,2);cm(2)-sin(wi*o.dt)*X0(:,1)+cos(wi*o.dt)*X0(:,2)];
  Xrig(:,j) = X1(:);
end

end% extractRHSrigidk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X] = filterHighFreq(o,X)
N = size(X,1)/2;
nv = size(X,2);
idx = ifftshift(-N/2:N/2-1);
idx = (abs(idx)>=N/4);
for i = 1:nv
  Xi = X(:,i);
  x = Xi(1:N);
  y = Xi(N+1:end);
  
  coeff = fft(x);
  coeff(idx) = 0;
  x = ifft(coeff);

  coeff = fft(y);
  coeff(idx) = 0;
  y = ifft(coeff);

  Xi = [x;y];
  X(:,i) = Xi;
end
end % filterHighFreq

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sigma] = filterHighFreqTen(o,sigma)
N = size(sigma,1);
nv = size(sigma,2);
idx = ifftshift(-N/2:N/2-1);
idx = (abs(idx)>=N/4);
for i = 1:nv
  sigmai = sigma(:,i);
  
  coeff = fft(sigmai);
  coeff(idx) = 0;
  sigmai = ifft(coeff);


  sigma(:,i) = sigmai;
end
end % filterHighFreqTen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ben,Ten,Div] = computeUpDerivs(o,vesicle);

N = vesicle.N;
nv = vesicle.nv;
uprate = o.uprate;
Nup = N*uprate;

Ben = zeros(2*Nup,2*Nup,nv);
Ten = zeros(2*Nup,Nup,nv);
Div = zeros(Nup,2*Nup,nv);

oc = curve;
[x,y] = oc.getXY(vesicle.X);
xup = interpft(x,Nup); yup = interpft(y,Nup);
Xup = [xup;yup];
[sa,xt] = oc.diffProp(Xup);
for k = 1:nv
  % compute single arclength derivative matrix
%  arcDeriv = repmat(o.isa(:,k),1,o.N).*D1;
  isa = 1./sa(:,k);
  arcDeriv = isa(:,ones(Nup,1)).*o.D1up;
  
  % This is much faster than repmat as it doesn't have to do all the
  % checks in repmat

  D4 = arcDeriv*arcDeriv; D4 = D4*D4;
  Ben(:,:,k) = [D4 zeros(Nup); zeros(Nup) D4];

  Ten(:,:,k) = [arcDeriv*diag(xt(1:Nup,k));...
               arcDeriv*diag(xt(Nup+1:end,k))];

  Div(:,:,k) = [diag(xt(1:Nup,k))*arcDeriv ...
                diag(xt(Nup+1:end,k))*arcDeriv];
end

Ben = real(Ben);
Ten = real(Ten);
Div = real(Div);
end % computeUpDerivs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Benup,Tenup,Divup] = upSampleDerivs(o,Ben,Ten,Div)
Nup = size(Ben,1)/2;
nv = size(Ben,3);
uprate = o.uprate;
N = Nup/uprate;

Benup = zeros(2*N,2*N,nv);
Tenup = zeros(2*N,N,nv);
Divup = zeros(N,2*N,nv);
for k = 1:nv
  Benup(:,:,k) = [o.downMx zeros(N,Nup);zeros(N,Nup) o.downMx]*Ben(:,:,k)*[o.upMx zeros(Nup,N);zeros(Nup,N) o.upMx];
  Tenup(:,:,k) = [o.downMx zeros(N,Nup);zeros(N,Nup) o.downMx]*Ten(:,:,k)*o.upMx;
  Divup(:,:,k) =  o.downMx*Div(:,:,k)*[o.upMx zeros(Nup,N);zeros(Nup,N) o.upMx];
end
%Benup = real(Benup);
%Tenup = real(Tenup);
%Divup = real(Divup);

end % upSampleDerivs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [upSampleMx] = upSampleMatrix(o,N,uprate);
fftmx = dftmtx(N);
ifftmx = conj(dftmtx(N*uprate))/(N*uprate);
upCoeffMx = sparse(N*uprate,N);
for i = 1:N/2
  upCoeffMx(i,i) = uprate;
end
for i = (N/2+1):N
  upCoeffMx(N*uprate-N+i,i) = uprate;
end

upSampleMx = real(ifftmx*upCoeffMx*fftmx);
%upSampleMx = ifftmx*upCoeffMx*fftmx;
end %upSampleMatrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [downSampleMx] = downSampleMatrix(o,N,uprate);
fftmx = dftmtx(N*uprate);
ifftmx = conj(dftmtx(N))/N;
downCoeffMx = sparse(N,N*uprate);
for i = 1:N/2
downCoeffMx(i,i) = 1/uprate;
end
for i = (N/2+1):N
downCoeffMx(i,N*uprate-N+i) = 1/uprate;
end
downSampleMx = real(ifftmx*downCoeffMx*fftmx);
%downSampleMx = ifftmx*downCoeffMx*fftmx;
end %downSampleMatrix

end % methods

end % classdef
