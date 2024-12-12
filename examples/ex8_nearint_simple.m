function  ex8_nearint_simple() 

clear all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization

% Physics parameters
N=64;
prams.N = N;               % points per vesicle
prams.nv = 1;               % number of vesicles
prams.T = 10;               % time horizon
prams.m = 200;              % number of time steps
prams.kappa = 1e-1;         % bending coefficient

options.farField = 'shear'; % background velocity
options.order = 2;          % time stepping order
options.vesves = 'implicit';
% Discretization of vesicle-vesicle interactions.
% Either 'explicit' or 'implicit'
options.near = true;        % near-singular integration
options.stockeval = true;        % near-singular integration with new code
options.fmm = false;
options.logFile = 'output/ex8_nearint.log';
% Name of log file for saving messages
options.dataFile = 'output/ex8_nearint.bin';
% Name of binary data file for storing vesicle information
options.axis = [6 6 -5 5];
% Axis for the plot

[options,prams] = initVes2D(options,prams);
% Set options and parameters that the user doesn't
% Also add src to path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Near singular part


t=(1:N)'*2*pi/N;
Xsou=[cos(t);0.3*sin(t)];
f=[sin(t);2*cos(t)];
Xsou2=zeros(2*N,2);
f2=zeros(2*N,2);
Xsou2(:,1)=[cos(t);0.3*sin(t)];
Xsou2(:,2)=[3.+cos(t);0.3*sin(t)];
f2(:,1)=[sin(t);2*cos(t)];
f2(:,2)=[1.+sin(t);2*cos(t)];


ntar=500;
numtar=10;
t=(1:ntar)'*2*pi/ntar;
Xtar=zeros(2*ntar,numtar);
for i=1:numtar
    Xtar(:,i)=(1.-(i)/(2.*numtar))*[cos(t);0.3*sin(t)];
end

vesicleSou = capsules(Xsou,[],[],0.1,1);
vesicleSou2 = capsules(Xsou2,[],[],0.1,1);
vesicleTar = capsules(Xtar,[],[],0.1,1);
vesSou=vesicleSou;
vesTar=vesicleTar;
ntar=vesTar.N;
tEq=false;
fcur=f;
op = poten(N);
kernel = @op.exactStokesDL;
tic;
% selfMat = op.stokesSLmatrix(vesSou);
[~, NearO] = vesSou.getZone(vesTar,2);
rhs=zeros(2*N,1);
fcur=zeros(2*N,1);
% poiseulle flow
 rhs(1:N,1)=Xsou(N+1:2*N).*(1-Xsou(N+1:2*N));
 rhs(N+1:2*N,1)=0.0;

% cubic flow
% rhs(1:N,1)=Xsou(N+1:2*N).^3;
% rhs(N+1:2*N,1)=Xsou(1:N).^3;

% fcur=selfMat\rhs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% New code
addnewpath();
s.Z = @(t) cos(t)+1i*0.3*sin(t); 
s.Zp = @(t) -sin(t)+1i*0.3*cos(t);
s.Zpp = @(t) -cos(t)-1i*0.3*sin(t);
% R = @(t) abs(s.Z(t));
% s.inside = @(z) abs(z)<=R(angle(z));
s.a = -0.001;  % interior pt far from bdry
s = quadr(s,N);

fcur=zeros(N,1);
fcur(:)=1+2i;


% tau=f(1:N)+1i*f(N+1:2*N);
% p.x=zeros(ntar*numtar,1);
% for i=1:numtar
%     p.x((i-1)*ntar+1:i*ntar)=Xtar(1:ntar,i)+1i*Xtar(ntar+1:2*ntar,i);
% end

tic;
% [uc,ug] = StokesScloseeval(p, vesicleSou, f, 'e');
% a.x=vesicleTar.x;
ucn=nearIntStockeval(vesSou,fcur, @StokesDcloseeval, ...
    vesTar,tEq);
ucni=ucn(1:ntar,:)+1i*ucn(ntar+1:end,:);

plot(Xsou(1:N),Xsou(N+1:end),'-k', 'LineWidth',3);
hold on;
for i=1:numtar
    plot(Xtar(1:ntar,i),Xtar(ntar+1:end,i));
    hold on;
end
% figure;
% scatter(Xtar(1:ntar,1),Xtar(ntar+1:end,1),10.1,abs(ucni(:,1)+(1+2i)),'filled');
figure;
plot_log_image([1,numtar],[1,ntar],abs(ucni+(1+2i)));
title('Error in double layer evaluation inside the ellipse');
xlabel('# of ellipse');
ylabel('Point on ellepse');

% imagesc(real(vesTar.x(:,1)),imag(vesTar.x(:,1)),abs(ucni(:,1)-(1+2i)))
% [uc,ug] = StokesScloseeval(a, vesicleSou, f, 'e');
fprintf('Timing new code: %g\n',toc)
% fc=selfMat*f;
% fprintf('comp. selfMat: %g %g\n',max(real2cmplx(fc)-uc));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compearing
%resold=zeros(ntar*numtar,1);
%for i=1:numtar
%    resold((i-1)*ntar+1:i*ntar)=LP(1:ntar,i)+1i*LP(ntar+1:2*ntar,i);
%end

for i=1:numtar

%     curpoleNew=ucn(2*(i-1)*ntar+1:2*i*ntar)';    
%    fprintf('Max error in target vecicle for new code #%d = %g\n',i, ...
%        max(abs(curpoleNew(ntar+1:2*ntar)     )  ) ... 
%     );
    true_pole=-(1+2i);
    fprintf('Values in tagrets for new code #%d: max=%g,min=%g\n',i, ...
     max(abs(ucni(:,i)-true_pole )), min(abs(ucni(:,i)-true_pole ))   );
 
 
 end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions

function uc = nearIntStockeval(vesicleSou,f, ...
    kernel,vesicleTar,tEquals)
    
    nvsou=vesicleSou.nv;
    nvtar=vesicleTar.nv;
    N=vesicleTar.N;
    uc=zeros(2*N,nvtar);
    if tEquals
         t.x=zeros(N*(nvtar-1),1);
    else   
         t.x=zeros(N*nvtar,1);
         idxs=1:nvtar;
    end
    for i=1:nvsou
            if tEquals
                idxs=[(1:i-1) (i+1:nvsou)];
            end
            t.x(:)=vesicleTar.x(:,idxs);
            SouCur.x=vesicleSou.x(:,i);
            SouCur.xp=vesicleSou.xp(:,i);
            SouCur.sp=vesicleSou.sa(:,i);
            SouCur.nx=vesicleSou.nx(:,i);
            SouCur.tang=vesicleSou.tang(:,i);
            SouCur.t=vesicleSou.t;
            SouCur.x=vesicleSou.x(:,i);
            SouCur.cw=vesicleSou.cw(:,i);
            SouCur.w=vesicleSou.w(:,i);
            SouCur.a=vesicleSou.a(i);
            fcur=f(:,i);
            uccur = kernel(t.x, SouCur, fcur, 'i');
            for j=1:numel(idxs)
                uc(1:N,idxs(j))=uc(1:N,idxs(j))+...
                   real( uccur((j-1)*N+1:j*N) );
                uc(N+1:2*N,idxs(j))=uc(N+1:2*N,idxs(j))+...
                   imag( uccur((j-1)*N+1:j*N) );
            end
    end     

end

function Xc=real2cmplx(X)
n = size(X,1)/2;
Xc=X(1:n)+1i*X(n+1:2*n);
end

function X=cmplx2real(Xc)
n = size(Xc,1);
X=zeros(2*n,1);
X(1:n)=real(Xc);
X(n+1:2*n)=imag(Xc);
end

function addnewpath()
P = path; ii = find(pwd == filesep); ii = ii(end);
subPath = pwd; subPath = [subPath(1:ii) ['src' filesep 'stokeseval_code']];
if isempty(strfind(P, subPath)),addpath(subPath);end
end

function s = quadr(s, N)  % set up periodic trapezoid quadrature on a segment
% Note sign change in normal vs periodicdirpipe.m.  Barnett 4/21/13
t = (1:N)'/N*2*pi; s.x = s.Z(t);
s.xp = s.Zp(t);
s.sp = abs(s.Zp(t)); 
s.nx = -1i*s.Zp(t)./s.sp;
s.tang = s.Zp(t)./s.sp;
s.cur = -real(conj(s.Zpp(t)).*s.nx)./s.sp.^2; 
s.w = 2*pi/N*s.sp; % speed weights
s.t = t; 
s.cw = 1i*s.nx.*s.w;  % complex weights (incl complex speed)
end

function pot=stokes_dl_smooth_eval(trg,src,den,nor)

nTrg = size(trg,1);
pot  = zeros(nTrg,2);

% inefficient (for demonstration only)
for iT=1:nTrg
    rx = trg(iT,1)-src(:,1);
    ry = trg(iT,2)-src(:,2);
    r  = [rx ry];
    rho2 = rx.*rx+ry.*ry;
    c   = dot(r,nor,2).*dot(r,den,2)./rho2./rho2; 
    r   = [c.*rx c.*ry];
    pot(iT,:) = sum(r)/pi;
end

end

function D = stokes_dl_singular_matrix(src_trg,nor,speed,tan,curvature)
%src_trg is both source and target
nSrc   = size(src_trg,1);
weight = speed*(2*pi/nSrc)/pi;%quadrature weight times the normalizing factor
weight = repmat(weight,2,2)'; 
D      = zeros(2*nSrc);

for iT=1:nSrc
    rx = src_trg(iT,1) - src_trg(:,1);
    ry = src_trg(iT,2) - src_trg(:,2);
    rho2 = rx.*rx+ry.*ry;
   
%     c   = sum(rx.*nor(:,1)+ry.*nor(:,2),2)./rho2./rho2; 
    c   = (rx.*nor(:,1)+ry.*nor(:,2))./rho2./rho2; 
    I   = [iT, iT+nSrc];
    D(I,:) = [c.*rx.*rx c.*rx.*ry;c.*ry.*rx c.*ry.*ry]'; %dyadic product
    D(I,I) = -curvature(iT)/2*(tan(iT,:)'*tan(iT,:)); %singular limit
    D(I,:) = weight.*D(I,:);    
end

end


function plot_log_image(x,y,z,Contours)

if nargin<4, Contours=10.^(-16:2);end
minmax = [min(Contours) max(Contours)];
imagesc(x(:),y(:),log10(z),log10(minmax));
caxis(log10(minmax));
ylims = log10([min(z(:)) max(z(:))]);
colorbar('FontSize',12,'YTick',log10(Contours(1:4:end)),'YTickLabel',Contours(1:4:end),'ylim',ylims);
colormap(jet);
end