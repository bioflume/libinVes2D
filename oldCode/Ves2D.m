function Xfinal = Ves2D(X,Xwalls,prams,options,Xtra,pressTar)
% Ves2D does time stepping on the intial configuration X
% with parameters and options defined in prams and options.
% Also can pass a set of initial tracer locations (Xtra) and 
% target points where one wants the pressure and stress (pressTar)

global matvecs; % number of matvecs
global derivs;  % number of times we compute differential
                % operators for preconditioning
global fmms     % number of fmm calls
global fileCount

if nargin == 4
  Xtra = [];       % Xtra is positions of tracers for visualization
  pressTar = [];   % points to compute pressure for postprocessing
elseif nargin == 5
  pressTar = [];   
end

matvecs = 0; % counter for the total number of time steps
derivs = 0;  % counter for the total number of times derivative
             % operators are computed
fmms = 0;    % counter for the totoal number of FMM calls
nreject = 0; % count for the number of rejected time steps
naccept = 0; % count for the number of accepted time steps
fileCount = 0;

N = prams.N; % Number of points per vesicle
nv = prams.nv; % Number of vesicles
Nbd = prams.Nbd; % number of points per solid wall
nvbd = prams.nvbd; % number of solid wall components

om = monitor(X,Xwalls,options,prams);
% Create class for doing input/output

if options.profile
  profile off; profile on;
end
% Turn on profiler

tstart = tic;
% Start a timer
c = clock;
% get current time in case tic and toc give funny anwers
message = ['Initial Month is ' num2str(c(2),'%d')];
om.writeMessage(message,'%s\n')
message = ['Initial Day is ' num2str(c(3),'%d')];;
om.writeMessage(message,'%s\n')
message = ['Initial Hour is ' num2str(c(4),'%d')];;
om.writeMessage(message,'%s\n')
message = ['Initial Minute is ' num2str(c(5),'%d')];;
om.writeMessage(message,'%s\n')
message = ['Initial Second is ' num2str(round(c(6)),'%d')];;
om.writeMessage(message,'%s\n')
message = ' ';
om.writeMessage(message,'%s\n')

tt = tstep(options,prams);
% build an object of class tstep with required options and parameters
Xstore = zeros(2*N,nv,options.order);
sigStore = zeros(N,nv,options.order);
uStore = zeros(2*N,nv,options.order);
etaStore = zeros(2*prams.Nbd,prams.nvbd,options.order);
RSstore = zeros(3,prams.nvbd,options.order);
% need to store options.order previous time steps to do
% higher-order time stepping
% RSstore is the rotlets and stokeslets stored as
% [stokeslet1;stokeslet2;rotlet]
Xstore(:,:,1) = X;
% initial configuration from user

if ~options.confined
  uStore(:,:,1) = tt.farField(X);
  etaStore = [];
  RSstore = [];
  walls = [];
else
  walls = tt.initialConfined(prams,Xwalls); 
end
% initial velocity on the vesicles is set to the background
% velocity.  If flow is unbounded, there is no density 
% function eta.  If bounded, compute a structure for the solid
% walls 

[Xstore,sigStore,uStore,etaStore,RSstore,Xtra] = ...
    tt.firstSteps(options,prams,...
    Xstore(:,:,1),sigStore(:,:,1),uStore(:,:,1),...
    walls,om,Xtra,pressTar);
% For higher-order methods (only 2nd order for now), need to 
% initially take smaller time steps so that we have two initial
% conditions

time = (options.order-1)*tt.dt;
% initial time.  firstSteps took the first time steps so that
% there is enough data to use the time stepping order that is
% desired

if strcmp(tt.vesves,'explicit') 
  sigStore = 0*sigStore;
  % this is just so that we can reproduce the results from the
  % 2dnear paper.  Without this, the next time step will use
  % the tension coming from doing the initial schur complement
  % but this was only introduced once sdc came into the
  % picture
end

% Main time stepping loop
Res = [];
accept = true;
while time < prams.T - 1e-10
  if time+tt.dt > prams.T
    tt.dt = prams.T - time;
  end
  dt = tt.dt;
  tt.currentTime = time;
  time = time + dt; % current time

  if options.order == 1
    [X,sigma,u,eta,RS,iter,accept,dtScale,res,iflag] = ...
        tt.timeStepGL(Xstore,sigStore,uStore,...
            etaStore,RSstore,prams.kappa,...
            prams.viscCont,walls,om,time,accept);
    % Do a single time step with internal substeps
    % corresponding to tt.orderGL
  else
    updatePreco = true;
    [X,sigma,u,eta,RS,iter,iflag] = tt.timeStep(...
        Xstore,sigStore,uStore,etaStore,RSstore,...
        [],[],[],[],[],[],[],...
        prams.kappa,prams.viscCont,walls,updatePreco);
    accept = true;
    dtScale = 1;
    res = 0;
    % if doing second or higher order time stepping,
    % can't use the timeStepGL, yet
  end

  if options.correctShape
    oc = curve;
    X = oc.correctAreaAndLength(X,om.area,om.length);
  end
  % Replace the shape with the closest one that has the same area and
  % length



  if options.tracers
    vel = tt.tracersVel(X,sigma,u,...
        prams.kappa,prams.viscCont,walls,eta,RS,Xtra);
    Xtra = Xtra + tt.dt*vel;
  end
  if ~accept
    time = time - dt;
  end
  % go back to old time

  for k = 1:options.order-1
    Xstore(:,:,k) = Xstore(:,:,k+1);
    sigStore(:,:,k) = sigStore(:,:,k+1);
    uStore(:,:,k) = uStore(:,:,k+1);
    etaStore(:,:,k) = etaStore(:,:,k+1);
    RSstore(:,:,k) = RSstore(:,:,k+1);
  end
  % Save the time steps still required if using second
  % order time stepping
  Xstore(:,:,options.order) = X;
  sigStore(:,:,options.order) = sigma;
  uStore(:,:,options.order) = u;
  etaStore(:,:,options.order) = eta;
  RSstore(:,:,options.order) = RS;
  % update the positions, tension, and velocity field of
  % the vesicles, and the density function and rotlet and
  % stokeslets

  if accept
    Res = [Res res(end)];
    % save the residual of the current time step
    terminate = om.outputInfo(X,sigma,u,Xwalls,Xtra,time,iter,...
        dtScale,res,iflag);
    % check if we have violated the error in area or length
    % also plot and save the current solution to dat file.
    % Print information to log file and console

    if terminate
      Xfinal = X;
      return;
    end
    % if error in area or length is too large, stop simulation
    naccept = naccept + 1;
   
    if options.pressure
      op = poten(N);
      [press,stress1,stress2] = op.pressAndStress(...
          X,sigma,u,prams.kappa,prams.viscCont,...
          walls,pressTar,eta,RS,options.confined,...
          options.near,options.fmm,om);
    end % options.pressure
    % compute the pressure and stress due to the vesicles
    % and the solid walls
  else
    nreject = nreject + 1;
  end % if accept
  % save data if solution was accepted, compute pressure and stress
end
% end of main time stepping loop


totTime = toc(tstart);
% Total time doing time stepping
c = clock;
% get current time in case tic and toc give funny anwers
message = ['Final Month is ' num2str(c(2),'%d')];
om.writeMessage(message,'%s\n')
message = ['Final Day is ' num2str(c(3),'%d')];;
om.writeMessage(message,'%s\n')
message = ['Final Hour is ' num2str(c(4),'%d')];;
om.writeMessage(message,'%s\n')
message = ['Final Minute is ' num2str(c(5),'%d')];;
om.writeMessage(message,'%s\n')
message = ['Final Second is ' num2str(round(c(6)),'%d')];;
om.writeMessage(message,'%s\n')
message = ' ';
om.writeMessage(message,'%s\n')

if options.profile
  profile off;
  filename = [options.logFile(1:end-4) 'Profile'];
  profsave(profile('info'),filename);
end
% Save profile

om.summary(matvecs,derivs,fmms,naccept,nreject,totTime);
% write a final summary with number of matvecs, accepts/rejects,
% and the CPU time

Xfinal = X;
% final configuation


