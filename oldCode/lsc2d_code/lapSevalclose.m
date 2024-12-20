function [u ux uy do] = lapSevalclose(x, s, tau, side) 
% LAPSEVALCLOSE - Laplace potential and field for SLP on curve (global quadr)
%
% u = lapSevalclose(x,s,tau,side) returns potentials at targets x due to a
%  single-layer potential with density tau living on a smooth curve structure
%  s, sampled with the global periodic trapezoid rule.
%  "side" controls whether targets are inside (default) or outside
%  (a mixture isn't allowed).
%  A scheme generalizing the Helsing/Ioakimidis globally compensated quadrature
%  is used (getting complex values in Kress style then using Cauchy formula).
%
% [u ux uy] = lapSevalclose(x,s,tau,side) also returns field (ie first
%  derivatives of potential)
%
% To set up a curve structure s with needed elements, given only a list of
%  complex nodes s.x, use code quadr.m in this directory.
%
% Inputs:
% x = M-by-1 list of targets in complex plane
% s = curve struct containing N-by-1 vector s.x of source nodes (as complex
%     numbers), and all other fields in s generated by quadr(), and
%     s.a one interior point far from bdry (mean(s.x) used if not provided).
% tau = single-layer density values at nodes
% side = 'i','e' to indicate targets are interior or exterior to the curve.
%
% Outputs:
% u     = potential values at targets x (M-by-1)
% ux,uy = (optional) first partials of potential at targets
% do    = (optional) diagnostic outputs:
%   do.vb : complex values on bdry (ie, v^+ for side='e' or v^- for side='i')
%
% complexity O(N^2) for evaluation of v^+ or v^-, plus O(NM) for globally
% compensated quadrature to targets (same as DLP)
%
% Also see:  QUADR, and enclosed test code testCSLPselfmatrix
%
% (c) Barnett 10/18/13; use cauchycompeval 10/23/13.
% Non-zero mean density 11/4/13 & finally v_infty fixed 11/21/13
% s.xp inserted as Bobbie did; using external quadr(), 10/8/14

if nargin<4, side='i'; end              % default
wantder = nargout>1;
if strcmp(x,'test') 
  testCSLPselfmatrix;
  return; 
end % SLP matrices (re,im) alone test
if ~isfield(s,'a') 
  s.a = mean(s.x);
  % warning('I''m guessing s.a for you. I sure hope it''s inside the
  % curve!');
end

vb = CSLPselfmatrix(s,side) * tau; 
% Step 1: eval v^+ or v^- = cmplx SLP(tau)
if side=='e'
  sawlog = s.t/1i + log(s.a - s.x(:));  
  % sawtooth with jump cancelled by log
  for i=1:numel(sawlog)-1, p = imag(sawlog(i+1)-sawlog(i)); 
  % remove phase jumps
    sawlog(i+1) = sawlog(i+1) - 2i*pi*round(p/(2*pi)); 
  end
  totchgp = sum(s.w.*tau)/(2*pi);
  % total charge due to SLP, over 2pi
  vb = vb + totchgp * sawlog;
  % is now r
  cw = 1i*s.nx.*s.w;
  % complex speed weights for native quadr...
  vinf = sum(1./(s.x-s.a).*vb.*cw) / (2i*pi); 
  % interior Cauchy gets v_infty
  %cauchycompeval(s.a,s,vb,'i');
  % sneaky interior Cauchy to get v_infty
  vb = vb - vinf;
  % kill off v_infty so that v is in exterior Hardy space
end
do.vb = vb;  % save for diagnostics (if totchg=0 it's useful)
%figure; plot([cumsum(s.w); cumsum(s.w)+sum(s.w)],[real(vb) imag(vb); real(vb) imag(vb)],'+-'); % test bdry data, & it's periodic

% Step 2: compensated close-evaluation of u = Re(v) & its 1st derivs..
if wantder 
  [v vp] = cauchycompeval(x,s,vb,side);
  if side=='e'
    vp = vp + totchgp./(s.a - x(:).');
  end 
  % add log deriv
  ux = real(vp)'; uy = -imag(vp)';
else 
  v = cauchycompeval(x,s,vb,side);
end
u = real(v)';
if side=='e'
  % add real part of log and of v_infty back in...
  u = u - totchgp*log(abs(s.a - x(:))) + real(vinf); 
  % Re v_infty = 0 anyway
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = CSLPselfmatrix(s,side) 
% complex SLP Kress-split Nystrom matrix s = src seg, even # nodes.
% side = 'i'/'e' int/ext case. Barnett 10/18/13.  only correct for
% zero-mean densities.

N = numel(s.x);
d = repmat(s.x, [1 N]) - repmat(s.x.', [N 1]);
% C-# displacements mat, t-s
x = exp(1i*s.t);
% unit circle nodes relative to which angles measured.
S = -log(d) + log(repmat(x,[1 N]) - repmat(x.',[N 1]));
% NB circ angles
S(diagind(S)) = log(1i*x./s.xp);               
% complex diagonal limit O(N^2) hack to remove 2pi phase jumps in S
% (assumes enough nodes for smooth):
for i=1:numel(S)-1 
  p = imag(S(i+1)-S(i));  
  % phase jump between pixels
  S(i+1) = S(i+1) - 2i*pi*round(p/(2*pi)); 
end 
% (NB access matrix as 1d array!)

%figure; imagesc(real(S)); figure; imagesc(imag(S)); stop % check S mat smooth!
  
m = 1:N/2-1;
if side=='e' 
  Rjn = ifft([0 1./m 1/N 0*m]); 
  % imag sign dep on side
else 
  Rjn = ifft([0 0*m 1/N 1./m(end:-1:1)]); 
end % cmplx Kress Rj(N/2)/4pi
%m = 1:N/2-1; 
% Rjn = ifft([0 0*m 1/N 1./m(end:-1:1)]); 
% Rjn = [Rjn(1) Rjn(end:-1:2)]; % flips order
S = S/N + circulant(Rjn); 
% incl 1/2pi SLP prefac. drops const part in imag
S = S .* repmat(s.sp.',[N 1]);  
% include speed (2pi/N weights already in)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function i = diagind(A)
% function i = diagind(A)
%
% return diagonal indices of a square matrix, useful for changing a
% diagonal in O(N) effort, rather than O(N^2) if add a matrix to A
% using matlab diag()
%
% barnett 2/6/08
N = size(A,1);
if size(A,2)~=N
  disp('input must be square!');
end
i = sub2ind(size(A), 1:N, 1:N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = circulant(x)
% function A = circulant(x)
%
% return square circulant matrix with first row x
% barnett 2/5/08
x = x(:);
A = toeplitz([x(1); x(end:-1:2)], x);

% ---------- TEST ROUTINES ------------

function testCSLPselfmatrix % test spectral conv of complex SLP self-int Nystrom
% This tests SLPselfmatrix and ISLPselfmatrix. Barnett 10/18/13
a = .3; w = 5;         % smooth wobbly radial shape params
R = @(t) 1 + a*cos(w*t); Rp = @(t) -w*a*sin(w*t); Rpp = @(t) -w*w*a*cos(w*t);
s.Z = @(t) R(t).*exp(1i*t); s.Zp = @(t) (Rp(t) + 1i*R(t)).*exp(1i*t);
s.Zpp = @(t) (Rpp(t) + 2i*Rp(t) - R(t)).*exp(1i*t);
fprintf('\n'); for N=40:40:500  % convergence study of Re v, Im v...
  s = quadr(s,N); sig = cos(3*s.t + 1); % zero-mean real SLP density func
  v = CSLPselfmatrix(s,'i') * sig; % holomorphic potential bdry value
  fprintf('N=%d:\tu(s=0) = %.15g  Im[v(s=0)-v(s=pi)] = %.15g\n',N,real(v(end)),imag(v(end)-v(end/2))); % note overall const for Im not high-order convergence
end % NB needs N=320 for 13 digits in Re, but 480 for Im part (why slower?)
%figure; plot(s.t, [real(v) imag(v)], '+-'); title('Re and Im of v=S\sigma');

% ---------- UNUSED USEFUL ROUTINES -------------

function S = SLPselfmatrix(s) % single-layer Kress-split self-int Nyst matrix
% s = src seg, even # nodes. Barnett 10/18/13 based upon mpspack/@layerpot/S.m
N = numel(s.x);
d = repmat(s.x, [1 N]) - repmat(s.x.', [N 1]);    % C-# displacements mat
S = -log(abs(d)) + circulant(0.5*log(4*sin([0;s.t(1:end-1)]/2).^2)); % peri log
S(diagind(S)) = -log(s.sp);                       % diagonal limit
m = 1:N/2-1; Rjn = ifft([0 1./m 2/N 1./m(end:-1:1)])/2; % Kress Rj(N/2)/4pi
S = S/N + circulant(Rjn); % includes SLP prefac 1/2pi. Kress peri log matrix L
S = S .* repmat(s.sp.',[N 1]);  % include speed factors (not 2pi/N weights)

