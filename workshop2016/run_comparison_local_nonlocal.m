% The code demonstrates that a nonlocal operator converges to a
% local operator with a decreasing epsilon.
%
% run_comparison_local_nonlocal(N, epsilon, niter)
% * Input
%   - N: # of elements
%   - epsilon: interaction radius
%   - niter: # of iteration to bisect epsilon
%
%     for iter=1:niter
%       solve local and nonlocal problem 
%       epsilon = epsilon / 2
%     end
%
% * Output 
%   - plot the result for each iteration
%
%  Author: Marta D'Elia
function run_comparison_local_nonlocal(N,epsilon,niter)
    close all;
    more off;
    %%% domains ---------------------------------------------------------------

    for iter=1:niter
        % local domain
        local_domain    = [ 0 1];
        nonlocal_domain = [ 0 1];
        [xN,xL,n,h] = domain(N,epsilon,nonlocal_domain,local_domain);
        
        if (epsilon < h)
            warning('nonlocal neighborhood (epsilon = %f) is smaller than  mesh size (h = %f)\n', epsilon, h);
        end
        
        %%% initializing the control ----------------------------------------------
        accuracy = true;
        print    = false;
        reuse    = false;
        test     = 5;
        
        theta_local_left  = exact_solution(xL(1),  test);
        theta_local_right = exact_solution(xL(end),test);
        
        fprintf('local problem is running\n');
        tic;
        [AL,bL,uL,errL] = local_problem(xL,h,epsilon,...
                                        theta_local_left,theta_local_right,...
                                        reuse,print,accuracy,test);
        toc;
        
        % dofs are doubled for discontinuous galerkin
        theta_nonlocal_left  = zeros(2*n,1);
        theta_nonlocal_right = zeros(2*n,1);
        
        theta_tmp = exact_solution(xN(1:n+1),test);
        theta_nonlocal_left(1:2:end) = theta_tmp(1:end-1);
        theta_nonlocal_left(2:2:end) = theta_tmp(2:end);
        
        theta_tmp = exact_solution(xN(end-n:end),test);
        theta_nonlocal_right(1:2:end) = theta_tmp(1:end-1);
        theta_nonlocal_right(2:2:end) = theta_tmp(2:end);
        
        fprintf('nonlocal problem is running\n');    
        tic;
        [AN,bN,uN,xN_plot,errN] = nonlocal_problem(xN,h,epsilon,n,...
                                                   theta_nonlocal_left,theta_nonlocal_right,...
                                                   reuse,print,accuracy,test);
        toc;
        
        plot(xL,uL,'k-+',xN_plot,uN,'r-*');
        pause;
        
        % difference at cell centers
        diff = zeros(N,1);
        for i=1:N
            pt = xL(i) + h/2;
            
            u_local    = evaluate_local_solution(i,pt,xL,uL);
            u_nonlocal = evaluate_nonlocal_solution(i+n,pt,xN,uN);
            
            diff(i) = abs(u_nonlocal - u_local);
        end
        fprintf('norm(diff(nonlocal,local)) = %f\n', norm(diff));

        epsilon = epsilon/2;
    end