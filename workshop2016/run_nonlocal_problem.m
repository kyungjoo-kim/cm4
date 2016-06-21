% The code solves a nonlocal 1D poisson problem. 
%
% [errN, AN, bN, uN] = run_nonlocal_problem(N, epsilon, test)
% * Input
%   - N: # of elements
%   - epsilon: interaction radius
%   - test: problem id
%     0 -> a source function with discontinuity
%     1 -> problem with a manufactured solution u = x
%     2 -> problem with a manufactured solution u = x^2
%     3 -> problem with a manufactured solution u = x^3
%     4 -> problem with a manufactured solution u = x^2 - x^4
% * Output 
%   - errN: nonlocal error with respect to manufactured solutions
%   - AN: nonlocal stiffness matrix
%   - bN: load vector
%   - uN: solution vector
%
%  Author: Marta D'Elia
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
    plot(xN_plot,uN,'r*');

    SN = sparse(AN);

    figure;
    spy(SN);

    fprintf('nnz(A) = %d\n', nnz(SN));    