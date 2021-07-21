%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SCHAFFER FUNCTION N. 4
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
% xx = [x1, x2]
% The function is usually evaluated on the square xi ∈ [-100, 100], for all i = 1, 2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef schaffer4
    properties
        D = 2;
        xbounds = [-100, 100;-100, 100];
        name = 'Schaffer n4';
        opt = 'max';
    end
    methods
        function y = do_eval(obj, xx)
            x1 = xx(1,:);
            x2 = xx(2,:);
            fact1 = (cos(sin(abs(x1.^2-x2.^2)))).^2 - 0.5;
            fact2 = (1 + 0.001*(x1.^2+x2.^2)).^2;
            
            y = 0.5 + fact1./fact2;
            
            if strcmp(obj.opt, 'max')
                y = -y;
            end
        end
    end
end
