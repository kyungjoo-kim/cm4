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
% x - grid points (assumed to be uniform)
% h - gird size
% epsilon - nonlocal operator neighborhood size
% n - points in epsilon
% theta_left - values at left boundary region
% theta_right - values at right boundary region
% print - plot result
% accuracy - compute error
% test - problem case
function [A,b,u,x_plot,err] = nonlocal_problem(x,h,epsilon,n,theta_left,theta_right,reuse,print,accuracy,test)
    
%%% nonlocal matrices and solving -----------------------------------------
    if (mod(epsilon,h)==0)
        shift=h;
    else
        shift = mod(epsilon,h);
    end  

    N = length(x)-2*n-1; %% active elements
    m = 2*(N+2*n); %% matrix size (double dof and epsilon in both sides)
    A = zeros(m, m);

    persistent A_stored;
    if (reuse && ~isempty(A_stored))
        A = A_stored;
    else
        for i  = 1:N+2*n
            for k = max(1-i,-n):min(n,N+2*n-i)
                if(i==1)
                    a  = x(1); 
                    b  = x(2);
                    hx = shift;
                elseif (i==N+2*n)
                    a  = x(end-1); 
                    b  = x(end);
                    hx = shift;
                else
                    hx = h;
                    a  = x(i)   + (k==n) *(h-shift);
                    b  = x(i+1) - (k==-n)*(h-shift);
                end
                if(i+k==1 || i+k==N+2*n)
                    hy=shift; 
                else
                    hy = h; 
                end
                I=i;
                dx_cv =  (x(i+k) - x(I))/hy; % shift for Change of Variables (CV)
                a = (a-x(I))/hx;
                b = (b-x(I))/hx;
                if(k==0)
                    flags=[1 1 1 1]; 
                    for ib = [2*i-1, 2*i]
                        for jb = [2*i-1,2*i]
                            A(ib,jb) = A(ib,jb)...
                                +hx*quadgk(@(v)gauss(mod(jb,2), mod(ib,2),max(0,hx/hy*v-epsilon/hy-dx_cv),hx/hy*v-dx_cv                        , v,hx/hy,epsilon,flags,dx_cv),a,b,'RelTol',1e-10)...
                                +hx*quadgk(@(v)gauss(mod(jb,2), mod(ib,2),  hx/hy*v-dx_cv,min((x(i+k+1)-x(i+k))/hy, hx/hy*v+epsilon/hy-dx_cv) , v,hx/hy,epsilon,flags,dx_cv),a,b,'RelTol',1e-10);
                        end
                    end
                else
                    for ib = [2*i-1, 2*i]
                        for jb = [2*(i+k)-1,2*(i+k)] 
                            Ib = floor((ib+1)/2);
                            Jb = floor((jb+1)/2);
                            flags = [i+k==Jb, Jb==I, i+k==Ib, Ib==I];           
                            t = hx*quadgk(@(v)gauss(mod(jb,2), mod(ib,2),  max(0, hx/hy*v-epsilon/hy-dx_cv),  min((x(i+k+1)-x(i+k))/hy, hx/hy*v+epsilon/hy-dx_cv), v,hx/hy,epsilon,flags,dx_cv),a,b,'RelTol',1e-10);            
                            A(ib,jb) = A(ib,jb) + t;
                            A(jb,ib) = A(jb,ib) + t; 
                        end
                    end
                    for ib = [2*(i+k)-1, 2*(i+k)]
                        for jb = [2*(i+k)-1, 2*(i+k)]
                            Ib = floor((ib+1)/2);
                            Jb = floor((jb+1)/2);
                            flags = [i+k==Jb, Jb==I, i+k==Ib, Ib==I];
                            A(ib,jb) = A(ib,jb) + hx*quadgk(@(v)gauss(mod(jb,2), mod(ib,2),  max(0, hx/hy*v-epsilon/hy-dx_cv),  min((x(i+k+1)-x(i+k))/hy, hx/hy*v+epsilon/hy-dx_cv), v,hx/hy,epsilon,flags,dx_cv),a,b,'RelTol',1e-10);                           
                        end
                    end          
                    for ib = [2*i-1, 2*i]
                        for jb = [2*i-1, 2*i]
                            Ib = floor((ib+1)/2);
                            Jb = floor((jb+1)/2);
                            flags = [i+k==Jb, Jb==I, i+k==Ib, Ib==I];
                            A(ib,jb) = A(ib,jb) + hx*quadgk(@(v)gauss(mod(jb,2), mod(ib,2),   max(0, hx/hy*v-epsilon/hy-dx_cv),  min((x(i+k+1)-x(i+k))/hy, hx/hy*v+epsilon/hy-dx_cv), v,hx/hy,epsilon,flags,dx_cv),a,b,'RelTol',1e-10);
                        end
                    end
                end
            end
        end
        A_stored = A;
    end
    
    %% epsilon region
    dof = 2*n; 
    
    A(1:dof, :   ) = 0;             A(end-dof+1:end,:)             = 0;
    A(1:dof,1:dof) = eye(dof);      A(end-dof+1:end,end-dof+1:end) = eye(dof);

    % forcing term
    xi = x(n+1:end-n);
    bi = zeros(N*2,1);
    for i=1:N
        bi(2*i-1) = bi(2*i-1) + quadgk(@(v)source_integral(i,  v,xi,h,epsilon,test),xi(i),xi(i+1),'RelTol',1e-10);
        bi(2*i)   = bi(2*i)   + quadgk(@(v)source_integral(i+1,v,xi,h,epsilon,test),xi(i),xi(i+1),'RelTol',1e-10);
    end 
    
    % solving
    b  = [theta_left; bi; theta_right];
    u  = A\b; % un in the domain and interaction domain
    ui = u(2*n+1:end); % ui in the domain only

    %%% plots -----------------------------------------------------------------
    x_plot = zeros(length(x(1:end-1))*2,1);
    x_plot(1:2:end) = x(1:end-1);
    x_plot(2:2:end) = x(2:end);
    % if (print)
    %     plot(x_plot,u,'*b','Linewidth',4), grid
    %     pause;
    % end

    %%% accuracy --------------------------------------------------------------
    err = 0;
    if (accuracy && test)
        % states
        value = 0;
        for i=1:N+2*n-1
            value = value + quadgk(@(v)nonlocal_error(i,v,x,u,test),x(i),x(i+1),'RelTol',1e-10);
        end
        err = sqrt(value);
    end


%----- FUNCTIONS -----%
function value = nonlocal_error(k,v,x,uN,test)
    value = (evaluate_nonlocal_solution(k,v,x,uN) - exact_solution(v,test)).^2;

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