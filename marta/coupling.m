function [thetaL,thetaN,uN,uL,errN,errL,errNT,errLT] = coupling(N,epsilon,test)
%*************************************************************************%
%                                                                         %
%  This function solves the optimization-based local-nonlocal coupling    %
%  problem using the built-in Matlab function fminunc.m                   %
%                                                                         %
%  Author: Marta D'Elia                                                   %
%                                                                         %
%  Modified: 01-06-2016                                                   %
%                                                                         %
%  NOTE 1: h is THE SAME for nonlocal and local problems                  %
%                                                                         %
%  NOTE 2: nonlocal discretization - piecewise linear Disc Galerkin       %
%          local discretization - piecewise linear Cont Galerkin          %
%                                                                         %
%*************************************************************************%
close all
accuracy = 1;
if(test==0), accuracy = 0; end

%%% domains ---------------------------------------------------------------
% nonlocal (0-\epsilon, 1+\epsilon)
b = 1;
a = 0;
h = (b-a)/N
if(mod(epsilon,h)==0) % # of discretization steps within \epsilon
    n = round(epsilon/h);
else
    n = ceil(epsilon/h);
end
xN = linspace(a-(n-1)*h,b+(n-1)*h,N+1+2*(n-1))'; % uniform spacing
xN = [a-epsilon; xN; b+epsilon]; % adding end points
% local (0.75, 1.75)
aL = 0.75;
bL = 1.75;
xL = linspace(aL,bL,N+1); 

%%% initializing the control ----------------------------------------------
thetaL = 1.0;
xc     = xN(end-n:end);
thetaN = zeros((length(xc)-1)*2,1); % this can be any function
thetaN = thetaN + 1.5;
theta  = [thetaL;thetaN]; % control vector
t0     = theta;

%%% solving --------------------------------------------------------------- 
% computing the nonlocal matrix once an for all
[uN,uL,A_full,uN_plot,xN_plot]= state(N,epsilon,thetaN,thetaL,0,0,0,test);
figure;
plot(xL,uL,'k-',xN_plot,uN_plot,'r-','Linewidth',4)
tic;
options = optimset('Display','iter','TolFun',1.e-10,'TolX',1.e-10,'MaxIter',1000,'LargeScale','off','MaxFunEvals',1000,'GradObj','on'); 
[theta,fval,info,output] = fminunc(@(t)functional(t,N,epsilon,h,xL,aL,b,n,A_full,1,test),t0,options);
output
thetaL = theta(1);
thetaN = theta(2:end);
[uN,uL,A_full,uN_plot,xN_plot,xL,errN,errL,errNT]=state(N,epsilon, ...
                                                  thetaN,thetaL,0,accuracy,A_full,test);
toc;
figure;
plot(xL,uL,'k-',xN_plot,uN_plot,'r-','Linewidth',4)



%----- FUNCTIONS -----%

function [Jvalue,gvalue] = functional(t,N,epsilon,h,xL,aL,b,n,A_full,grad,test)
    %%% computes functional and gradient 
    % states
    [uN,uL] = state(N,epsilon,t(2:end),t(1),0,0,A_full,test);
    % overlap domain
    xb = (linspace(aL,b+(n-1)*h,(b+(n-1)*h-aL)/h+1))';
    xb = [xb;b+epsilon];
    Nb = length(xb);
    % selection of the DOFs of ul on the overlap
    iL     = floor((b+epsilon-aL)/h)+1;
    ul_end = (uL(iL+1)-uL(iL))./(xL(iL+1)-xL(iL)).*(b+epsilon-xL(iL)) + (uL(iL));
    ul     = [uL(1:iL);ul_end];
    % projection of ul onto the discontinuous FE space
    ulDG = ul(1);
    for k=2:length(ul)-1, ulDG = [ulDG;ul(k);ul(k)]; end
    ulDG = [ulDG;ul(end)];
    % selection of the DOFs of un on the overlap
    iN = floor((aL)/h) + 1;
    un = uN(2*iN-1:end);
    % computation of J
    Jvalue = 0;
    for i = 1:Nb-1
        a      = xb(i); 
        b      = xb(i+1);
        Jvalue = Jvalue + quadgk(@(v)gaussJ(i,v,xb,ulDG,un),a, b,'RelTol',1e-10);
    end
    Jvalue = 0.5*Jvalue;
    
    % adjoint and gradient
    gvalue = t;
    if(grad)
        [lambdaN,lambdaL] = adjoint(N,epsilon,A_full,uN,uL);
        gvalue = [lambdaL;lambdaN];
        % gvalue = gradient(N,epsilon,lambdaL,lambdaN);
    end

function value = gaussJ(k,v,x,uL,uN)  
    %%% evaluation of the L2 norm |un-ul|^2
    ul = ((uL(2*k)-uL(2*k-1))./(x(k+1)-x(k)).*(v'-x(k)) + (uL(2*k-1)))';
    un = ((uN(2*k)-uN(2*k-1))./(x(k+1)-x(k)).*(v'-x(k)) + (uN(2*k-1)))';
    value = (un-ul).^2;
     