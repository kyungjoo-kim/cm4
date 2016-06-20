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
function [errL,AL,bL,uL] = run_local_problem(N,epsilon,test)
    close all;
    more off;

    %%% domains ---------------------------------------------------------------

    % non local dummy
    nonlocal_domain = 0;

    if (test == 0)
        %% for integration of source term it needs nonzero epsilon
    else
        %% otherwise, epsilon is dummy and should be zero
        epsilon = 0;
    end
    
    % local domain
    local_domain = [ 0 1.75 ];
    [xN,xL,n,h] = domain(N,epsilon,nonlocal_domain,local_domain);

    %%% initializing the control ----------------------------------------------
    accuracy = true;
    print    = false;
    reuse    = false;

    theta_local_left  = exact_solution(xL(1),  test);
    theta_local_right = exact_solution(xL(end),test);

    tic;
    [AL,bL,uL,errL] = local_problem(xL,h,epsilon,...
                                    theta_local_left,theta_local_right,...
                                    reuse,print,accuracy,test);
    toc;
    
    figure;
    plot(xL,uL,'k-*');

    SL = sparse(AL);

    figure;
    spy(SL);

    printf('nnz(A) = %d\n', nnz(SL));    