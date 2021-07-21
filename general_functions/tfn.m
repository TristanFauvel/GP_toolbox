function [value, dvalue_dx, dvalue_dfx] = tfn ( x, fx )

%*****************************************************************************80
%
%% TFN calculates the T-function of Owen.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    20 January 2008
%
%  Author:
%
%    Original FORTRAN77 version by JC Young, Christoph Minder.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    MA Porter, DJ Winstanley,
%    Remark AS R30:
%    A Remark on Algorithm AS76:
%    An Integral Useful in Calculating Noncentral T and Bivariate
%    Normal Probabilities,
%    Applied Statistics,
%    Volume 28, Number 1, 1979, page 113.
%
%    JC Young, Christoph Minder,
%    Algorithm AS 76:
%    An Algorithm Useful in Calculating Non-Central T and
%    Bivariate Normal Distributions,
%    Applied Statistics,
%    Volume 23, Number 3, 1974, pages 455-457.
%
%  Parameters:
%
%    Input, real X, FX, the parameters of the function.
%
%    Output, real VALUE, the value of the T-function.
%
ng = 5;

r = [ ...
    0.1477621, ...
    0.1346334, ...
    0.1095432, ...
    0.0747257, ...
    0.0333357 ];
tp = 0.159155;
tv1 = 1.0E-35;
tv2 = 15.0;
tv3 = 15.0;
tv4 = 1.0E-05;
u = [ ...
    0.0744372, ...
    0.2166977, ...
    0.3397048, ...
    0.4325317, ...
    0.4869533 ];
%
%  Test for X near zero.
%

if numel(x) ~= numel(fx)
    error('x and fx must have the same size')
end
n = numel(x);
value = NaN(size(x));
value(abs(x)<tv1) = tp*atan(fx(abs(x)<tv1));
value(tv2 < abs(x)) = 0.0; %  Test for large values of abs (X).
value(abs(fx)< tv1) = 0.0;%  Test for FX near zero.

id = 1:n;
id = id(isnan(value));

xs = NaN(1,n);
x1 = NaN(1,n);
x2 = NaN(1,n);
fxs = NaN(1,n);
for i= id%  Test whether abs ( FX ) is so large that it must be truncated.
    %
    xs(i) = - 0.5 * x(i).* x(i);
    x2(i) = fx(i);
    fxs(i) = fx(i).* fx(i);
    %
    %  Computation of truncation point by Newton iteration.
    %

        if (tv3 <= log (1.0 + fxs(i) ) - xs(i)*fxs(i) )

            x1(i) = 0.5 * fx(i);
            fxs(i) = 0.25 * fxs(i);

            while ( true )

                rt = fxs(i) + 1.0;

                x2(i) = x1(i) + ( xs(i).*fxs(i) + tv3 - log(rt)) ...
                    / ( 2.0 * x1(i) .* ( 1.0 ./ rt - xs(i) ) );

                fxs(i) = x2(i) .* x2(i);

                if ( abs ( x2(i) - x1(i) ) < tv4 )
                    break
                end

                x1(i) = x2(i);            
            end        
        end
    %
    %  Gaussian quadrature.
    %
    rt = 0.0;
    
    for j = 1 : ng
        
        r1 = 1.0 + fxs(i) .* ( 0.5 + u(j) )^2;
        r2 = 1.0 + fxs(i) .* ( 0.5 - u(j) )^2;
        
        rt = rt + r(j) * ( exp ( xs(i).* r1 ) ./ r1 + exp ( xs(i).* r2 ) ./ r2 );
        
    end
    
    value(i) = rt .* x2(i) .* tp;
end

if nargout > 1
     dvalue_dx = -(normcdf(x.*fx) - 0.5).*exp(-0.5*x.^2)/sqrt(2*pi);
     dvalue_dfx = exp(-0.5*x.^2.*(1+fx.^2))./(2*pi*(1+fx.^2));
%     dvalue_dx = -x.*(1+x.^2).*value;
%     dvalue_dfx = sqrt(2*pi)./x*(normcdf(x*x2) - 0.5);
%     
% %     dvalue_dx = -x.*(1+x.^2).*value;
% %     dvalue_dfx = sqrt(2*pi)./x*(normcdf(x*fx) - 0.5);
end
return
end

% figure(); 
% plot(x2); hold on;
% plot(x1); hold on;
% plot(fx); hold off;
