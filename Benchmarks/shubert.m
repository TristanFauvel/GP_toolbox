%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SHUBERT FUNCTION
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
%The function is usually evaluated on the square xi ∈ [-10, 10], for all i = 1, 2, although this may be restricted to the square xi ∈ [-5.12, 5.12], for all i = 1, 2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef shubert
    properties
        D = 2
        xbounds = [0, 10; 0, 10];
        name = 'Schubert';
        opt = 'max';
    end
    methods
        function y = do_eval(obj,xx)
            if size(xx,1)~=obj.D
                error('Problem with input size')
            end
            x1 = xx(1,:);
            x2 = xx(2,:);
            sum1 = 0;
            sum2 = 0;
            
            for ii = 1:5
                new1 = ii * cos((ii+1)*x1+ii);
                new2 = ii * cos((ii+1)*x2+ii);
                sum1 = sum1 + new1;
                sum2 = sum2 + new2;
            end
            
            y = sum1 .* sum2;
            
            if strcmp(obj.opt, 'max')
                y = -y;
            end
            
        end
    end
end

