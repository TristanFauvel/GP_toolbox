
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GRAMACY & LEE (2012) FUNCTION
%
% Authors: Sonja Surjanovic, Simon Fraser University
%          Derek Bingham, Simon Fraser University
% Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
%
% Copyright 2013. Derek Bingham, Simon Fraser University.
%
% THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
% FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
% derivative works, such modified software should be clearly marked.
% Additionally, this program is free software; you can redistribute it
% and/or modify it under the terms of the GNU General Public License as
% published by the Free Software Foundation; version 2.0 of the License.
% Accordingly, this program is distributed in the hope that it will be
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty
% of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.
%
% For function details and reference information, see:
% http://www.sfu.ca/~ssurjano/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef grlee12
    properties
            D=1;
            xbounds = [0.5,2.5];
        name = 'Gramacy and Lee (2012)';
                opt = 'max';

    end
    methods
        function y = do_eval(obj, xx)
            if size(xx,1)~=obj.D
                error('Problem with input size')
            end
            term1 = sin(10*pi.*xx) ./ (2.*xx);
            term2 = (xx-1).^4;
            
            y = term1 + term2;
            
            
            % if nargout>1
            % dydx = 4*(x - 1)^3 - sin(10*pi*x)/(2*x^2) + (5*pi*cos(10*pi*x))/x;
            % dydx = -dydx;
            %
            % end
            
            if strcmp(obj.opt, 'max')
                y = -y;
            end
            y(xx > obj.xbounds(:,2) | xx <  obj.xbounds(:,1)) = NaN;
        end
    end
end
