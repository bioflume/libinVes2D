function near_singular()

% construct boundary
n    = 256;
bdry = get_curve(n);

% evaluation points
ng = 128;
[xg,yg] = meshgrid(linspace(-1.5,1.5,ng),linspace(-.5,.5,ng));
trg = [xg(:) yg(:)];

% density (constant)
cnt = [1,2];
den = repmat(cnt,n,1);
den = den.*[bdry.speed bdry.speed]*2*pi/n; %trapezoidal quadrature 

% double layer evaluation
pot = stokes_dl_smooth_eval(trg,bdry.x,den,bdry.nor);
err = pot./repmat(cnt,ng*ng,1)+1;
en  = reshape(sqrt(dot(err,err,2)),ng,[]);

plot_log_image(xg,yg,en);
% colormap(goodbw)
title('relative error in double layer evaluation');
hold on;
plot(bdry.x(:,1),bdry.x(:,2),'-k', 'LineWidth',3);
axis equal tight;

% solving for density
U  = repmat(cnt,n,1);
qw = bdry.speed*2*pi/n;

D  = stokes_dl_singular_matrix(bdry.x,bdry.nor,bdry.speed,bdry.tan,bdry.curvature);
N  = bdry.nor(:)*(bdry.nor(:).*[qw;qw])';
% figure;
% [xmat,ymat] = meshgrid(linspace(-5,5,n),linspace(-5,5,n));
% imagesc(D(1:n,1:n));
% figure;
% imagesc(D(n+1:end,n+1:end));
% figure;
% plot(N);
D  = D + N; %rank modification to make D invertible
den = (-.5*eye(2*n)+D)\U(:);
figure;
plot(den+U(:));
title('Error in density inversion');
end

function curve=get_curve(n)
t = (0:n-1)'*2*pi/n;
a = .3;r=1;
curve.x = r*[cos(t) a*sin(t)];
curve.tan = r*[-sin(t) a*cos(t)];
curve.speed = sqrt(dot(curve.tan,curve.tan,2));
curve.tan = curve.tan./[curve.speed curve.speed];
curve.nor = [curve.tan(:,2) -curve.tan(:,1)];

xss = r*[-cos(t) -a*sin(t)];
curve.curvature = -dot(curve.nor,xss,2)./curve.speed.^2;
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
