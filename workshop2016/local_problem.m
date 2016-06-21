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
% h - element size 
% theta_left - left boundary value
% theta_right - right boundary value
% print - plot solution
% accuracy - compute error when it is on
% test - problem case
function [A,b,u,err] = local_problem(x,h,epsilon,theta_left,theta_right,reuse,print,accuracy,test)

    hidden_flag = false;
    %%% matrix --------------------------------------------------------------
    N = length(x);
    A = zeros(N, N);

    persistent A_stored;
    if (reuse && ~isempty(A_stored))
        A = A_stored;
    else
        for i = 1:N
            A(i,i) = 2/h;  
            if (i~=1)
                A(i,i-1) = -1/h; 
            end
            if (i~=N)
                A(i,i+1) = -1/h;
            end
        end
        A_stored = A;
    end
    A(1,:) = 0;
    A(1,1) = 1;
    A(end,:) = 0;
    A(end,end) = 1;

    % forcing term
    b = zeros(N,1);
    epsilon = epsilon*(test==0);
    for i=1:N-1
        b(i)   = b(i)     + quadgk(@(v)source_integral(i,  v,x,h,epsilon,test,hidden_flag),x(i),x(i+1),'RelTol',1e-10);
        b(i+1) = b(i+1)   + quadgk(@(v)source_integral(i+1,v,x,h,epsilon,test,hidden_flag),x(i),x(i+1),'RelTol',1e-10);
    end

    % boundary modification
    b(1)   = theta_left;
    b(end) = theta_right;

    % solve problem
    u = A\b;

    %%% accuracy --------------------------------------------------------------
    err = 0;
    if (accuracy && test)
        value = 0;
        for i=1:N-1
            value = value + quadgk(@(v)local_error(i,v,x,u,test),x(i),x(i+1),'RelTol',1e-10);
        end
        err = sqrt(value);
    end
    
    %%% plots -----------------------------------------------------------------
    % if (print)
    %     plot(x,u,'g--','Linewidth',4);
    %     pause;
    % end

%----- FUNCTIONS -----%
function value = local_error(k,v,x,uL,test)
    value = (evaluate_local_solution(k,v,x,uL) - exact_solution(v,test)).^2;

