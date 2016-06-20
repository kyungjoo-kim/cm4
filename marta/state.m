function [uN,uL,A_full,u_full,xN_plot,xL,errN,errL,errNT,l,wN,A] = state(N,epsilon,thetaN,thetaL,print,accuracy,A_full,test)
%*************************************************************************%
%                                                                         %
%  This function computes the local and nonlocal states given the values   %
%  the volume constraint and the boundary condition                        %
%                                                                         %
%  Author: Marta D'Elia                                                   %
%                                                                         %
%  Modified: 06-01-2016                                                   %
%                                                                         %
%  NOTE 1: h is THE SAME for nonlocal and local problems                  %
%                                                                         %
%  NOTE 2: nonlocal discretization - piecewise linear Disc Galerkin       %
%          local discretization - piecewise linear Cont Galerkin          %
%                                                                         %
%*************************************************************************%

%%% domains ---------------------------------------------------------------
% nonlocal (0-\epsilon, 1+\epsilon)
bN  = 1;
aN  = 0;
h  = (bN-aN)/N;
if(mod(epsilon,h)==0)
    shift=h;
    n=round(epsilon/h);
else
    shift = mod(epsilon,h);
    n = ceil(epsilon/h);
end
xN = linspace(aN-(n-1)*h,bN+(n-1)*h,N+1+2*(n-1))';
xN = [aN-epsilon; xN; bN+epsilon];
% local (0.75, 1.75) 
aL = 0.75;
bL = 1.75;
xL = linspace(aL,bL,N+1);
hL = (bL-aL)/N;

%%% local matrices and solving --------------------------------------------
AL  = zeros(N+1, N+1);  
for i = 1:N+1
      AL(i,i)   =  2/hL;  
    if(i~=1)
      AL(i,i-1) = -1/hL; 
    end
    if(i~=N+1)
      AL(i,i+1) = -1/hL;
    end
end
AL(1,:) = 0;
AL(1,1) = 1;
AL(end,:) = 0;
AL(end,end) = 1;
% right boundary conditions
switch test
    case 0
        theta_right = 0; 
    case 1
        theta_right=1.75;
    case 2
        theta_right=1.75^2;
    case 3
        theta_right=1.75^3;
    case 4
        theta_right=1.75^2-1.75^4;
end

% forcing term
wL = zeros(N+1,1);
for i = 1:N
    a       = xL(i); 
    b       = xL(i+1);
    wL(i)   = wL(i)     + quadgk(@(v)fL(i,v,xL,hL,epsilon,test),a,b,'RelTol',1e-10);
    wL(i+1) = wL(i+1)   + quadgk(@(v)fL(i+1,v,xL,hL,epsilon,test),a,b,'RelTol',1e-10);
end
wL(1) = thetaL;
wL(end) = theta_right;
uL  = AL\wL;

%%% nonlocal matrices and solving -----------------------------------------
if(A_full==0)   
    A      = zeros(2*N+4*n,2*N+4*n); 
    l      = zeros(2*N+4*n,1);       
    for i  = 1:N+2*n
     for k = max(1-i,-n):min(n,N+2*n-i)
       if(i==1)
          a  = xN(1); 
          b  = xN(2);
          hx = shift;
       elseif (i==N+2*n)
          a  = xN(end-1); 
          b  = xN(end);
          hx = shift;
       else
          hx = h;
          a  = xN(i)   + (k==n) *(h-shift);
          b  = xN(i+1) - (k==-n)*(h-shift);
       end
       if(i+k==1 || i+k==N+2*n)
              hy=shift; 
       else
              hy = h; 
       end
       I=i;
       dx_cv =  (xN(i+k) - xN(I))/hy; % shift for Change of Variables (CV)
       a = (a-xN(I))/hx;
       b = (b-xN(I))/hx;
        if(k==0)
             flags=[1 1 1 1]; 
             for ib = [2*i-1, 2*i]
                for jb = [2*i-1,2*i]
                        A(ib,jb) = A(ib,jb)...
                        +hx*quadgk(@(v)gauss(mod(jb,2), mod(ib,2),max(0,hx/hy*v-epsilon/hy-dx_cv),hx/hy*v-dx_cv                        , v,hx/hy,epsilon,flags,dx_cv),a,b,'RelTol',1e-10)...
                        +hx*quadgk(@(v)gauss(mod(jb,2), mod(ib,2),  hx/hy*v-dx_cv,min((xN(i+k+1)-xN(i+k))/hy, hx/hy*v+epsilon/hy-dx_cv) , v,hx/hy,epsilon,flags,dx_cv),a,b,'RelTol',1e-10);
                end
             end
        else
         for ib = [2*i-1, 2*i]
           for jb = [2*(i+k)-1,2*(i+k)] 
               Ib = floor((ib+1)/2);
               Jb = floor((jb+1)/2);
               flags = [i+k==Jb, Jb==I, i+k==Ib, Ib==I];           
                  t = hx*quadgk(@(v)gauss(mod(jb,2), mod(ib,2),  max(0, hx/hy*v-epsilon/hy-dx_cv),  min((xN(i+k+1)-xN(i+k))/hy, hx/hy*v+epsilon/hy-dx_cv), v,hx/hy,epsilon,flags,dx_cv),a,b,'RelTol',1e-10);            
                  A(ib,jb) = A(ib,jb) + t;
                  A(jb,ib) = A(jb,ib) + t; 
           end
         end
         for ib = [2*(i+k)-1, 2*(i+k)]
           for jb = [2*(i+k)-1, 2*(i+k)]
               Ib = floor((ib+1)/2);
               Jb = floor((jb+1)/2);
               flags = [i+k==Jb, Jb==I, i+k==Ib, Ib==I];
                    A(ib,jb) = A(ib,jb) + hx*quadgk(@(v)gauss(mod(jb,2), mod(ib,2),  max(0, hx/hy*v-epsilon/hy-dx_cv),  min((xN(i+k+1)-xN(i+k))/hy, hx/hy*v+epsilon/hy-dx_cv), v,hx/hy,epsilon,flags,dx_cv),a,b,'RelTol',1e-10);                           
           end
         end          
        for ib = [2*i-1, 2*i]
           for jb = [2*i-1, 2*i]
               Ib = floor((ib+1)/2);
               Jb = floor((jb+1)/2);
                flags = [i+k==Jb, Jb==I, i+k==Ib, Ib==I];
                  A(ib,jb) = A(ib,jb) + hx*quadgk(@(v)gauss(mod(jb,2), mod(ib,2),   max(0, hx/hy*v-epsilon/hy-dx_cv),  min((xN(i+k+1)-xN(i+k))/hy, hx/hy*v+epsilon/hy-dx_cv), v,hx/hy,epsilon,flags,dx_cv),a,b,'RelTol',1e-10);
           end
         end
       end
     end
    end
    A_full = A;
else
    A      = A_full;
end
A(1:2*n,:)                     = zeros(size(A(1:2*n,:)));
A(1:2*n,1:2*n)                 = eye(2*n,2*n);
A(end-2*n+1:end,:)             = zeros(size(A(end-2*n+1:end,:)));
A(end-2*n+1:end,end-2*n+1:end) = eye(2*n,2*n);

% forcing term
xN_int = linspace(aN,bN,N+1);
wN     = zeros(N*2,1);
for i  = 1:N
    a         = xN_int(i); 
    b         = xN_int(i+1);
    wN(2*i-1) = wN(2*i-1) + quadgk(@(v)fN(i,  v,xN_int,h,epsilon,test),a,b,'RelTol',1e-10);
    wN(2*i)   = wN(2*i)   + quadgk(@(v)fN(i+1,v,xN_int,h,epsilon,test),a,b,'RelTol',1e-10);
end 

% left volume constraint pointwise
% exact control for testing state.m: see at the bottom of the file
switch test
    case 0
        t_left = xN(1:n+1).*0;
    case 1
        t_left = xN(1:n+1);
    case 2
        t_left = xN(1:n+1).^2;
    case 3
        t_left = xN(1:n+1).^3;
    case 4
        t_left = xN(1:n+1).^2 - xN(1:n+1).^4;
end
% left volume constraint for DG DOFs
theta_left = t_left(1);
for k=2:length(t_left)-1
    theta_left=[theta_left;t_left(k);t_left(k)];
end
theta_left=[theta_left;t_left(end)];

% solving
rhs                = zeros(2*N+4*n,1);
rhs(1:2*n)         = theta_left;
rhs(2*n+1:end-2*n) = wN;
rhs(end-2*n+1:end) = thetaN;
u_full             = A\rhs; % un in the domain and interaction domain
uN                 = u_full(2*n+1:end); % un in the domain only

%%% plots -----------------------------------------------------------------
xN_plot = 0;
xN = linspace(aN-(n-1)*h,bN+(n-1)*h,N+1+2*(n-1))';
xN = [aN-epsilon; xN; bN+epsilon];
xN_plot(1) = xN(1);
for i=2:length(xN)-1
xN_plot=[xN_plot; xN(i);xN(i)];
end
xN_plot=[xN_plot;xN(end)];
if(print)
    plot(xL',uL,'g--','Linewidth',4);hold on;
    plot(xN_plot,u_full,'*b','Linewidth',4), grid
    pause
end

%%% accuracy --------------------------------------------------------------
errN=0;
errL=0;
errNT=0;
if(accuracy)
    % states
    value = 0;
    for i = 1:length(xN)-1
        a     = xN(i); 
        b     = xN(i+1);
        value = value + quadgk(@(v)EN(i,v,xN,u_full,test),a, b,'RelTol',1e-10);
    end
    errN  = sqrt(value);
    value = 0;
    for i = 1:length(xL)-1
        a     = xL(i); 
        b     = xL(i+1);
        value = value + quadgk(@(v)EL(i,v,xL,uL,test),a, b,'RelTol',1e-10);
    end
    errL = sqrt(value);
    % theta
    value = 0;
    xC = xN(end-n:end);
    uC = u_full(end-2*n+1:end);
    for i = 1:length(xC)-1
        a     = xC(i); 
        b     = xC(i+1);
        value = value + quadgk(@(v)EN(i,v,xC,uC,test),a, b,'RelTol',1e-10);
    end
    errNT = sqrt(value);
end


%----- FUNCTIONS -----%


function value = phi(y_m_xi)
    %%% test function, piece-wise linear
    value = 1 + (y_m_xi) .* (2*(-y_m_xi>=0)-1);

function value = gauss(phiJ,phiI,a,b,xj,hratio,epsilon,flags,dx_cv)
   %%% matrix integrals, gaussian quadgk rules 
   y1    =  (b-a)/2  *    sqrt((3 - 2*sqrt(6/5))/7) + (a+b)/2;
   y2    = -(b-a)/2  *    sqrt((3 - 2*sqrt(6/5))/7) + (a+b)/2;
   y3    =  (b-a)/2  *    sqrt((3 + 2*sqrt(6/5))/7) + (a+b)/2;
   y4    = -(b-a)/2  *    sqrt((3 + 2*sqrt(6/5))/7) + (a+b)/2;
   value =  (b-a)/2 .*  ((18+sqrt(30))/36 * sub_f(y1,xj,dx_cv,phiJ,phiI,flags,hratio) + ... 
                         (18+sqrt(30))/36 * sub_f(y2,xj,dx_cv,phiJ,phiI,flags,hratio) + ...
                         (18-sqrt(30))/36 * sub_f(y3,xj,dx_cv,phiJ,phiI,flags,hratio) + ...
                         (18-sqrt(30))/36 * sub_f(y4,xj,dx_cv,phiJ,phiI,flags,hratio)   )/epsilon^2;

function value = sub_f(v,xj,dx_cv,phiJ,phiI,flags,hratio)
        %%% sub integrand - matrix
        phi1_y  = 0;
        phi1_xj = 0;
        phi2_y  = 0;
        phi2_xj = 0;
    if(flags(1))
        phi1_y  = (phiJ*(1-v)+(1-phiJ)*v);
    end
    if(flags(2))
        phi1_xj = (phiJ*(1-xj)+(1-phiJ)*xj);
    end
        phi1    = phi1_y-phi1_xj;
    if(flags(3))
       phi2_y   = (phiI*(1-v)+(1-phiI)*v);  
    end
    if(flags(4))
        phi2_xj = (phiI*(1-xj)+(1-phiI)*xj);
    end
       phi2     = phi2_y - phi2_xj;
       value    = (phi1.*phi2)./(abs(v-hratio*xj+dx_cv)+1.e-9);

function value = fL(i,v,x,h,epsilon,test)
    %%% local source integrals 
    switch test
        case 0
            f = (v<0.5-epsilon).*0 + ...
            (v>=0.5-epsilon & v<0.5).*(-(2/epsilon^2).* (0.5*epsilon^2-epsilon+3/8+(2*epsilon-3/2-log(epsilon)).*v + (3/2+log(epsilon)).*v.^2-(v.^2-v).*log(0.5-v) )) + ...
            (v>0.5 & v<0.5+epsilon).*(-(2/epsilon^2).*(0.5*epsilon^2-epsilon-3/8+(2*epsilon+3/2+log(epsilon)).*v -(3/2+log(epsilon)).*v.^2+(v.^2-v).*log(v-0.5) )) + ...
            (v>=0.5+epsilon).*(-2);
        case 1
            f = 0.*v;
        case 2
            f = v.*0-2;
        case 3
            f = -6.*v;
        case 4
            f = 12.*v.^2-2;
    end
    ph_v  = 0;
    if(abs(v-x(i))<=h)
    ph_v  = phi((v-x(i))./h);
    end
    value = f.*ph_v;
    
function value = fN(i,v,x,h,epsilon,test)
    %%%  nonlocal source integrals
    switch test
        case 0
            f = (v<0.5-epsilon).*0 + ...
            (v>=0.5-epsilon & v<0.5).*(-(2/epsilon^2).* (0.5*epsilon^2-epsilon+3/8+(2*epsilon-3/2-log(epsilon)).*v + (3/2+log(epsilon)).*v.^2-(v.^2-v).*log(0.5-v) )) + ...
            (v>0.5 & v<0.5+epsilon).*(-(2/epsilon^2).*  (0.5*epsilon^2-epsilon-3/8+(2*epsilon+3/2+log(epsilon)).*v - (3/2+log(epsilon)).*v.^2+(v.^2-v).*log(v-0.5) )) + ...
            (v>=0.5+epsilon).*(-2);
        case 1
            f = 0.*v;
        case 2
            f = v.*0-2;
        case 3
            f = -6.*v;
        case 4
            f = 12.*v.^2-2 + epsilon^2;
    end
    ph_v  = 0;
    if(abs(v-x(i))<=h)
    ph_v  = phi((v-x(i))./h);
    end    
    value = f.*ph_v;
   
function value = EN(k,v,x,uN,test)
    %%% nonlocal error
    un = ((uN(2*k)-uN(2*k-1))./(x(k+1)-x(k)).*(v'-x(k)) + (uN(2*k-1)))';
    switch test
        case 0
            value = 0;
        case 1
            value = (un-v).^2;
        case 2
            value = (un-v.^2).^2;
        case 3
            value = (un-v.^3).^2;
        case 4
            value = (un-v.^2+v.^4).^2;
    end           


function value = EL(k,v,x,uL,test)
    %%% local error
    ul = ((uL(k+1)-uL(k))./(x(k+1)-x(k)).*(v'-x(k)) + (uL(k)))';
    switch test
        case 0
            value = 0;
        case 1
            value = (ul-v).^2;
        case 2
            value = (ul-v.^2).^2;
        case 3
            value = (ul-v.^3).^2;
        case 4
            value = (ul-v.^2+v.^4).^2;
    end


%__________________________________________________________________________
%
% exact control for testing state.m, ignore for the coupling
% % thetaN_right = x^2/x^3, for testing the state
% x_right = xN(end-n:end).^2;
% x_right = xN(end-n:end).^2 - xN(end-n:end).^4;
% x_right = xN(end-n:end)-1;
% thetaN_right = x_right(1);
% for k=2:length(x_right)-1
%     thetaN_right=[thetaN_right;x_right(k);x_right(k)];
% end
% thetaN_right=[thetaN_right;x_right(end)];
% thetaN = thetaN_right;
%__________________________________________________________________________

