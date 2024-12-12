set(0,'DefaultAxesFontSize',22)
options.confined = false;  % confined flow

if 0
  file = 'ex1_shear1VesData.bin';
  ax = [-6 6 -5 5];
end
if 0
  file = 'ex2_shear4VesData.bin';
  ax = [-6 6 -5 5];
end
if 0
  file = 'ex3_relaxationData.bin';
  ax = [-3 3 -4 4];
end
if 1
  file = 'flowextensionalN64nv2ts0.0125visc16384order1implicitinterpGLorder2nsdc0resolveCol0minSep0.bin';
  ax = [-2 2 -1 1];
end
if 0
  file = 'ex5_taylorGreenData.bin';
  ax = [0 pi 0 pi] + 0.1*[-1 1 -1 1];
end
if 0
  file = 'ex6_couetteData.bin';
  ax = [-21 21 -21 21];
end
if 0
  file = 'ex7_doubleCouetteData.bin';
  ax = [-21 21 -21 21];
end
if 0
  file = 'flowboxN64nv3Nbd256ts0.02visc100order1implicitinterpGLorder2nsdc0resolveCol0minSep0gCont80.bin';
  file = 'flowboxN64nv3Nbd256ts0.02visc100order1explicitinterpGLorder2nsdc0resolveCol1minSep1gCont80.bin';
  ax = [-1.5 1.5 -3.0 3.0]; 
end
if 0
  file = 'flowchokeN64nv6Nbd256ts0.000625visc100order1explicitinterpGLorder2nsdc0resolveCol1minSep1utest.bin';
  ax = [-10 10 -3.0 3.0]; 
end

[posx,posy,ten,eta,RS,wallx,wally,ea,el,time,n,nv] = loadFile(file);
% load positions, tension, errors, time, number of points,
% and number of vesicles

options.quiver = false;
% can add a quiver plot on the vesicles if desired

ntime = numel(time);
% number of time steps

options.savefig = false;
if options.savefig
  count = 1;
end
% for saving images to png files

for k=1:1:ntime
  if time(k) < 6
    continue;
  end
  vec1 = [posx(:,:,k);posx(1,:,k)];
  vec2 = [posy(:,:,k);posy(1,:,k)];
  figure(5); clf
  plot(vec1,vec2,'-x','linewidth',1)
  axis equal
  axis(ax);
  titleStr = ['t = ' num2str(time(k),'%4.2e') ...
      ' eA = ' num2str(ea(k),'%4.2e') ...
      ' eL = ' num2str(el(k),'%4.2e')];
  title(titleStr)

  if options.quiver && k~=ntime
    velx = (posx(:,:,k+1) - posx(:,:,k))/(time(k+1)-time(k));
    vely = (posy(:,:,k+1) - posy(:,:,k))/(time(k+1)-time(k));
    hold on
    quiver(posx(:,:,k),posy(:,:,k),velx,vely,'k')
    hold off
  end
  % Quiver plot of a finite difference applied to position
  
  if options.confined
    hold on
    vec1 = [wallx;wallx(1,:)];
    vec2 = [wally;wally(1,:)];
    plot(vec1,vec2,'k','linewidth',3);
    hold off
  end
%  set(gca,'visible','off')

  if options.savefig
    filename = ['./frames/image', sprintf('%04d',count),'.png'];
    count = count+1;
    figure(1);
    print(gcf,'-dpng','-r300',filename);
  end

  %pause(0.01)
  pause
end


