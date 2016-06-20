function [lambdaN,lambdaL,vn,vl] = adjoint(N,epsilon,An,uN,uL)
%*************************************************************************%
%                                                                         %
%  This function computes the loc and nl ADJOINT of the coupling problem  %
%                                                                         %
%  Author: Marta D'Elia                                                   %
%                                                                         %
%  Modified: 06-02-2016                                                   %
%                                                                         %
%  NOTE 1: h is THE SAME for nonlocal and local problems                  %
%                                                                         %
%  NOTE 2: nonlocal discretization - piecewise linear Disc Galerkin       %
%          local discretization - piecewise linear Cont Galerkin          %
%                                                                         %
%  NOTE 3: uN and uL are the extended solutions                           %
%                                                                         %
%*************************************************************************%

%%% domains ---------------------------------------------------------------
% nonlocal (0-\epsilon, 1+\epsilon)
b  = 1;
a  = 0;
h  = (b-a)/N;
if(mod(epsilon,h)==0)
    n=round(epsilon/h);
else
    n = ceil(epsilon/h);
end
xN = linspace(a-(n-1)*h,b+(n-1)*h,N+1+2*(n-1))';
xN = [a-epsilon; xN; b+epsilon];
% local (0.75, 1.75)
aL = 0.75;
bL = 1.75;
hL = (bL-aL)/N;
xL = linspace(aL,bL,N+1);

%%% forcing term ----------------------------------------------------------
% vector of points on \Omega_b
xb = (linspace(aL,b+(n-1)*h,(b+(n-1)*h-aL)/h+1))';
xb = [xb;b+epsilon];
Nb = length(xb);
% selection of the local solution on \Omega_b
iL = floor((b+epsilon-aL)/h)+1;
% addition of a point on b+epsilon (nonlocal grid)
ul_end = (uL(iL+1)-uL(iL))./(xL(iL+1)-xL(iL)).*(b+epsilon-xL(iL)) + (uL(iL));
ul     = [uL(1:iL);ul_end];
% extension of the local vector to the discontinuous FE space
ulDG = ul(1);
for k=2:length(ul)-1
    ulDG=[ulDG;ul(k);ul(k)];
end
ulDG = [ulDG;ul(end)];
% selection of the nonlocal solution on \Omega_b
iN = floor((aL)/h) + 1;
un = uN(2*iN-1:end);

% nonlocal forcing term
wN    = zeros((Nb-1)*2,1);
for i = 1:Nb-1
    a         = xb(i); 
    b         = xb(i+1);
    wN(2*i-1) = wN(2*i-1)   + quadgk(@(v)gaussF_nl(i,  v,xb,ulDG,un,b-a,0),a, b,'RelTol',1e-10);
    wN(2*i)   = wN(2*i)     + quadgk(@(v)gaussF_nl(i+1,v,xb,ulDG,un,b-a,1),a, b,'RelTol',1e-10);
    
end
wn = zeros(size(An,1),1);
wn(end-2*(Nb-1)+1:end) = wN;

% local forcing term
wL    = zeros(Nb,1);
for i = 1:Nb-1
    a       = xb(i) ;
    b       = xb(i+1);
    wL(i)   = wL(i)   + quadgk(@(v)gaussF_l(i,  v,xL,xb,ul,un,h,0),a,b,'RelTol',1e-10);
    wL(i+1) = wL(i+1) + quadgk(@(v)gaussF_l(i+1,v,xL,xb,ul,un,h,1),a,b,'RelTol',1e-10);
end

wl = zeros(N+1,1);    
wl(1:Nb) = wL;

%%% solving ---------------------------------------------------------------
% local 
Al  = zeros(N+1, N+1); 
for i = 1:N
      Al(i,i)   =  2/hL;  
    if(i~=1)
      Al(i,i-1) = -1/hL; 
    end
    if(i~=N)
      Al(i,i+1) = -1/hL;
    end
end
Al(:,1)     = zeros(N+1,1);
Al(1,1)     = 1;
Al(N+1,:)  = zeros(1,N+1);
Al(:,end)   = zeros(N+1,1);
Al(end,end) = 1;
pl          = Al\wl;
lambdaL     = pl(1);
vl          = pl;
vl(1)       = 0;
% nonlocal
An(:,1:2*n)         = zeros(size(An(:,1:2*n)));
An(1:2*n,:)         = zeros(size(An(1:2*n,:)));
An(1:2*n,1:2*n)     =  eye(size(An(1:2*n,1:2*n)));
An(:,end-2*n+1:end) = zeros(size(An(:,end-2*n+1:end)));
Mass                = zeros(2*n,2*n);
An(end-2*n+1:end,end-2*n+1:end) = eye(2*n,2*n);
pn                  = An\wn;
lambdaN             = pn(end-2*n+1:end);
vn                  = pn;
vn(end-2*n+1:end)   = zeros(2*n,1);


%----- FUNCTIONS -----%


function value = phi(y_m_xi)
   %%% test function, piece-wise linear 
   value = 1 + (y_m_xi) .* (2*(-y_m_xi>=0)-1);

function value = gaussF_nl(i,v,x,uL,uN,h_current, on)%h_uniform,aL)
    %%%  nonlocal source integrals 
    ph_v = 0;
    if(abs(v-x(i))<=h_current)
        ph_v = phi((v-x(i))./h_current);
    end    
    k     = i-on;
    ul    = ((uL(2*k)-uL(2*k-1))./(x(k+1)-x(k)).*(v'-x(k)) + (uL(2*k-1)))';
    un    = ((uN(2*k)-uN(2*k-1))./(x(k+1)-x(k)).*(v'-x(k)) + (uN(2*k-1)))';
    value = (un-ul).*ph_v;
    
function value = gaussF_l(j,v,x,xb,uL,uN,h,on)
    %%%  local source integrals 
    ph_v  = 0;
    if(abs(v-x(j))<=h)
        ph_v = phi((v-x(j))./h);
    end
    k     = j - on;
    un    = ((uN(2*k)-uN(2*k-1))./(xb(k+1)-xb(k)).*(v'-xb(k)) + (uN(2*k-1)))';
    ul    = ((uL(k+1)-uL(k))./(xb(k+1)-xb(k)).*(v'-xb(k)) + (uL(k)))';
    value = (ul-un).*ph_v; 
   
function value = gaussMass(i,j,v,x,h_current)
      phi1 = 0;
      phi2 = 0;
   if(abs(v-x(i))<=h_current)
      phi1 = phi((v-x(i))./h_current);
   end
   if(abs(v-x(j))<=h_current)
      phi2 = phi((v-x(j))./h_current);
   end
   value = phi1.*phi2;
