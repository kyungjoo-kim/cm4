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
% for local operator epsilon should be zero
function value = source_integral(i,v,x,h,epsilon,test,flag)
    if nargin > 7
        error('source_integral: wrong number of input');
    end
    if (nargin == 7)
    else
        flag = true;
    end
    
    switch test
      case 0
        f = (v<0.5-epsilon).*0 + ...
            (v>=0.5-epsilon & v<0.5).*(-(2/epsilon^2).* (0.5*epsilon^2-epsilon+3/8+(2*epsilon-3/2-log(epsilon)).*v + (3/2+log(epsilon)).*v.^2-(v.^2-v).*log(0.5-v) )) + ...
            (v>0.5 & v<0.5+epsilon).*(-(2/epsilon^2).*(0.5*epsilon^2-epsilon-3/8+(2*epsilon+3/2+log(epsilon)).*v -(3/2+log(epsilon)).*v.^2+(v.^2-v).*log(v-0.5) )) + ...
            (v>=0.5+epsilon).*(-2);
      case 1
        f = v.*0;
      case 2
        f = v.*0-2;
      case 3
        f = -6.*v;
      case 4
        f = 12.*v.^2-2 + epsilon^2;
      case 5
        %% arbitrary source term (p>=3) for comparing local and
        %% nonlocal operators
        f = (v-.5).^3;
      case 6
        if (flag == true)
            f = 2*sin(80*v/pi)./(v.^2); % -6.*v;
        else
            f = 0;
        end
    end
    ph_v  = 0;
    if (abs(v-x(i)) <= h)
        ph_v  = phi((v-x(i))./h);
    end
    value = f.*ph_v;
    