function [posx,posy,ten,eta,RS,wallx,wally,ea,el,time,n,nv] = ...
    loadFile(file)
fid = fopen(file,'r');
val = fread(fid,'double');
fclose(fid);

n = val(1);
nv = val(2);
nbd = val(3);
nvbd = val(4);
walls = val(7:7+2*nbd*nvbd-1);
val = val(7+2*nbd*nvbd:end);


ntime = numel(val)/(3*n*nv+3+(2*nbd+3)*nvbd);
if ntime ~= ceil(ntime);
  disp('PROBLEM WITH VALUES FOR n AND nv');
end

wallx = zeros(nbd,nvbd);
wally = zeros(nbd,nvbd);
for k = 1:nvbd
  istart = (k-1)*nbd+1;
  iend = k*nbd;
  wallx(:,k) = walls(istart:iend);
  wally(:,k) = walls((istart:iend)+nbd*nvbd);
end
%wallx = walls(1:nbd*nvbd);
%wally = walls(nbd*nvbd+1:2*nbd*nvbd);
posx = zeros(n,nv,ntime);
posy = zeros(n,nv,ntime);
ten = zeros(n,nv,ntime);
eta = zeros(2*nbd,ntime);
RS = zeros(3,ntime);
time = zeros(ntime,1);
ea = zeros(ntime,1);
el = zeros(ntime,1);

istart = 1;
for m = 1:ntime
  for k=1:nv
    iend = istart + n - 1;
    posx(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end
  % load x positions

  for k=1:nv
    iend = istart + n - 1;
    posy(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end
  % load y positions

  for k=1:nv
    iend = istart + n - 1;
    ten(:,k,m) = val(istart:iend);
    istart = iend + 1;
  end
  % load tensions
  if nbd
  iend = istart + 2*nbd - 1;
  eta(:,m) = val(istart:iend);
  istart = iend + 1;
  % load eta 
  
  iend = istart + 3 -1;
  RS(:,m) = val(istart:iend);
  istart = iend + 1;
  % load RS
  end
  ea(m) = val(istart);
  el(m) = val(istart+1);
  time(m) = val(istart+2);
  istart = istart + 3;

end



