function [Xfinal] = Ves2D(X,Xwalls,Xrig,prams,options,Xtra,pressTar)
% Ves2D does time stepping on the intial configuration X
% with parameters and options defined in prams and options.
% Also can pass a set of initial tracer locations (Xtra) and 
% target points where one wants the pressure and stress (pressTar)

global matvecs; % number of matvecs
global derivs;  % number of times we compute differential
                % operators for preconditioning
global fmms     % number of fmm calls
warning off
if nargin == 5
  Xtra = [];       % Xtra is positions of tracers for visualization
  pressTar = [];   % points to compute pressure for postprocessing
elseif nargin == 6
  pressTar = [];   
end

matvecs = 0; % counter for the total number of time steps
derivs = 0;  % counter for the total number of times derivative
             % operators are computed
fmms = 0;    % counter for the totoal number of FMM calls
nreject = 0; % count for the number of rejected time steps
naccept = 0; % count for the number of accepted time steps

N = prams.N; % Number of points per vesicle
nv = prams.nv; % Number of vesicles
Nbd = prams.Nbd; % number of points per solid wall
nvbd = prams.nvbd; % number of solid wall components
Nrd = prams.Nrd;
nvrd = prams.nvrd;

om = monitor(X,Xwalls,options,prams);
% Create class for doing input/output

if options.profile
  profile off; profile on -timer real;
%  profile off; profile on
end
% Turn on profiler

tstart = tic;
% Start a timer
c = clock;
% get current time in case tic and toc give funny anwers
message = ['Initial Month is ' num2str(c(2),'%d')];
om.writeMessage(message,'%s\n')
message = ['Initial Day is ' num2str(c(3),'%d')];
om.writeMessage(message,'%s\n')
message = ['Initial Hour is ' num2str(c(4),'%d')];
om.writeMessage(message,'%s\n')
message = ['Initial Minute is ' num2str(c(5),'%d')];
om.writeMessage(message,'%s\n')
message = ['Initial Second is ' num2str(round(c(6)),'%d')];
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
XrigStore = zeros(2*Nrd,nvrd,options.order);
etaRigStore = zeros(2*Nrd,nvrd,options.order);
% need to store options.order previous time steps to do
% higher-order time stepping
% RSstore is the rotlets and stokeslets stored as
% [stokeslet1;stokeslet2;rotlet]
Xstore(:,:,1) = X;
XrigStore(:,:,1) = Xrig;
% initial configuration from user

if ~options.confined
  uStore(:,:,1) = tt.farField(X);
  etaStore = [];
  RSstore = [];
  walls = [];
  wallsCoarse = [];
else
  [walls,wallsCoarse] = tt.initialConfined(prams,Xwalls); 
end
% initial velocity on the vesicles is set to the background velocity.
% If flow is unbounded, there is no density function eta.  If bounded,
% compute a structure for the solid walls 

if nv
  tt.withRigid = false;
  [Xstore,sigStore,uStore,etaStore,RSstore,Xtra] = ...
      tt.firstSteps(options,prams,...
      Xstore(:,:,1),sigStore(:,:,1),uStore(:,:,1),Xrig,Xrig*0,...
      walls,wallsCoarse,om,Xtra,pressTar);
  tt.withRigid = options.withRigid;
else
  if options.saveData
    fid = fopen(options.dataFile,'w');
    fwrite(fid,[N;nv;Nbd;nvbd;Nrd;nvrd],'double');
    fclose(fid);
  end
end
% For higher-order methods (only 2nd order for now), need to initially
% take smaller time steps so that we have two initial conditions

time = (options.order-1)*tt.dt;
% initial time.  firstSteps took the first time steps so that there is
% enough data to use the time stepping order that is desired

% Main time stepping loop
accept = true;
while time < prams.T - 1e-10
  if time+tt.dt > prams.T
    tt.dt = prams.T - time;
  end
  % make sure we land exactly on the time horizon
  dt = tt.dt;
  tt.currentTime = time;
  time = time + dt; % current time
 
  if options.order == 1
    pertic = tic;
    [X,sigma,u,eta,RS,Xrig,etaRig,iter,accept,dtScale,res,iflag] = ...
        tt.timeStepGL(Xstore,sigStore,uStore,...
            etaStore,RSstore,XrigStore,etaRigStore,prams.kappa,...
            prams.viscCont,walls,wallsCoarse,om,time,accept);
    pertoc = toc(pertic);
    om.writeMessage(['Whole time step takes: ' num2str(pertoc)],'%s\n');
    % Do a single time step with internal substeps
    % corresponding to tt.orderGL
  else
    updatePreco = true;
    [X,sigma,u,eta,RS,Xrig,etaRig,iter,iflag] = tt.timeStep(...
        Xstore,sigStore,uStore,etaStore,RSstore,XrigStore,etaRigStore,...
        [],[],[],[],[],[],[],...
        prams.kappa,prams.viscCont,walls,wallsCoarse,updatePreco);
    accept = true;
    dtScale = 1;
    res = 0;
    % if doing second or higher order time stepping,
    % can't use the timeStepGL, yet
  end


  if options.tracers
    vel = tt.tracersVel(X,sigma,u,...
        prams.kappa,prams.viscCont,walls,eta,RS,Xtra);
    Xtra = Xtra + tt.dt*vel;
  end

  if ~accept
    time = time - dt;
  end
  % go back to old time

  if accept
    oc = curve;
    if options.correctShape
      X = oc.correctAreaAndLength(X,om.area,om.length);
    end

    if options.redistributeArcLength
      [X,u,sigma] = oc.redistributeParameterize(X,u,sigma);
    end
    % redistribute the points equally in arclength
    warning off
    terminate = om.outputInfo(X,sigma,u,eta,RS,Xwalls,Xtra,...
        time,iter,dtScale,res,iflag,Xrig,etaRig);
    % check if we have violated the error in area or length
    % also plot and save the current solution to dat file.
    % Print information to log file and console

    if terminate
      Xfinal = X;
      totTime = toc(tstart);
      message = ['Final Month is ' num2str(c(2),'%d')];
      om.writeMessage(message,'%s\n')
      message = ['Final Day is ' num2str(c(3),'%d')];
      om.writeMessage(message,'%s\n')
      message = ['Final Hour is ' num2str(c(4),'%d')];
      om.writeMessage(message,'%s\n')
      message = ['Final Minute is ' num2str(c(5),'%d')];
      om.writeMessage(message,'%s\n')
      message = ['Final Second is ' num2str(round(c(6)),'%d')];
      om.writeMessage(message,'%s\n')
      message = ' ';
      om.writeMessage(message,'%s\n')
      om.summary(matvecs,derivs,fmms,naccept,nreject,totTime);
      % write a final summary with number of matvecs, accepts/rejects,
      % and the CPU time
      return;
    end
    % if error in area or length is too large, stop simulation
    naccept = naccept + 1;
   
    if options.pressure
      op = poten(N);
      [press,stress1,stress2] = op.pressAndStress(...
          X,sigma,u,prams.kappa,prams.viscCont,...
          walls,pressTar,eta,RS,options.confined,...
          options.fmm,om);
    end % options.pressure
    % compute the pressure and stress due to the vesicles
    % and the solid walls
  else
    nreject = nreject + 1;
  end % if accept
  % save data if solution was accepted, compute pressure and stress

  for k = 1:options.order-1
    Xstore(:,:,k) = Xstore(:,:,k+1);
    sigStore(:,:,k) = sigStore(:,:,k+1);
    uStore(:,:,k) = uStore(:,:,k+1);
    etaStore(:,:,k) = etaStore(:,:,k+1);
    RSstore(:,:,k) = RSstore(:,:,k+1);
    XrigStore(:,:,k) = XrigStore(:,:,k+1);
    etaRigStore(:,:,k) = etaRigStore(:,:,k+1);
  end
  % Save the time steps still required if using second
  % order time stepping
  Xstore(:,:,options.order) = X;
  sigStore(:,:,options.order) = sigma;
  uStore(:,:,options.order) = u;
  etaStore(:,:,options.order) = eta;
  RSstore(:,:,options.order) = RS;
  XrigStore(:,:,options.order) = Xrig;
  etaRigStore(:,:,options.order) = etaRig;

  % update the positions, tension, and velocity field of
  % the vesicles, and the density function and rotlet and
  % stokeslets

end
% end of main 


totTime = toc(tstart);
% Total time doing time stepping
c = clock;
% get current time in case tic and toc give funny anwers
message = ['Final Month is ' num2str(c(2),'%d')];
om.writeMessage(message,'%s\n')
message = ['Final Day is ' num2str(c(3),'%d')];
om.writeMessage(message,'%s\n')
message = ['Final Hour is ' num2str(c(4),'%d')];
om.writeMessage(message,'%s\n')
message = ['Final Minute is ' num2str(c(5),'%d')];
om.writeMessage(message,'%s\n')
message = ['Final Second is ' num2str(round(c(6)),'%d')];
om.writeMessage(message,'%s\n')
message = ' ';
om.writeMessage(message,'%s\n')

if options.profile
  profile off;
  p = profile('info');
  filename = [options.logFile(1:end-4) 'Profile'];
  save([filename '.mat'],'p');
%  profview
  profsave(profile('info'),filename);
end
% Save profile

om.summary(matvecs,derivs,fmms,naccept,nreject,totTime);
% write a final summary with number of matvecs, accepts/rejects,
% and the CPU time

Xfinal = X;
% final configuation

message = 'SIMULATION SUCCESSFULLY COMPLETED';
om.writeMessage(message,'%s\n')
om.writeStars


