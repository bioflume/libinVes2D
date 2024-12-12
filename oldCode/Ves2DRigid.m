function Xfinal = Ves2DRigid(X,Xwalls,prams,options,Xtra,pressTar)

tt = tstep(options,prams);
Xm = X;

time = 0;

oc = curve;
while(time<prams.T)
  time = time + tt.dt;
  particle = capsules(Xm,[],[],0,0);
  cm = particle.center
  [Xm,Up,omg] = tt.timeStepRigid(Xm,tt.dt);
  Up
  omg
  [x,y] = oc.getXY(Xm);
  plot([x;x(1,:)],[y;y(1,:)],'k-x');
  axis equal;
  axis(options.axis);
  pause(0.1)

end


