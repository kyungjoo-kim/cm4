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
function [errN,AN,bN,uN] = run_nonlocal_problem(N,epsilon,test)
    close all;
    more off;

    %%% domains ---------------------------------------------------------------

    % local dommy
    local_domain = 0;
    
    nonlocal_domain = [ 0 1.75];
    [xN,xL,n,h] = domain(N,epsilon,nonlocal_domain,local_domain);

    if (epsilon < h)
        warning('nonlocal neighborhood (epsilon = %f) is smaller than  mesh size (h = %f)\n', epsilon, h);
    end
    
    %%% initializing the control ----------------------------------------------
    accuracy = true;
    print    = false;
    reuse    = false;

    % dofs are doubled for discontinuous galerkin
    theta_nonlocal_left  = zeros(2*n,1);
    theta_nonlocal_right = zeros(2*n,1);
    
    theta_tmp = exact_solution(xN(1:n+1),test);
    theta_nonlocal_left(1:2:end) = theta_tmp(1:end-1);
    theta_nonlocal_left(2:2:end) = theta_tmp(2:end);
    
    theta_tmp = exact_solution(xN(end-n:end),test);
    theta_nonlocal_right(1:2:end) = theta_tmp(1:end-1);
    theta_nonlocal_right(2:2:end) = theta_tmp(2:end);

    tic;
    [AN,bN,uN,xN_plot,errN] = nonlocal_problem(xN,h,epsilon,n,...
                                       theta_nonlocal_left,theta_nonlocal_right,...
                                       reuse,print,accuracy,test);
    toc;
    
    figure;
    plot(xN_plot,uN,'r-+');

    SN = sparse(AN);

    figure;
    spy(SN);

    fprintf('nnz(A) = %d\n', nnz(SN));    