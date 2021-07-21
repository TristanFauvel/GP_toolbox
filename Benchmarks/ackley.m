%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ACKLEY FUNCTION
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
% INPUTS:
%
% xx = [x1, x2, ..., xd]
% a = constant (optional), with default value 20
% b = constant (optional), with default value 0.2
% c = constant (optional), with default value 2*pi
%The function is usually evaluated on the hypercube xi ∈ [-32.768, 32.768], for all i = 1, …, d, although it may also be restricted to a smaller domain.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef ackley
    properties
        D = 2;
        name = 'Ackley';
        c = 2*pi;
        b = 0.2;
        a = 20;
        xbounds
        opt = 'max';

    end
    methods
        function obj = ackley(D,a,b,c)
            if (nargin < 4)
                obj.c = 2*pi;
            end
            if (nargin < 3)
                obj.b = 0.2;
            end
            if (nargin < 2)
                obj.a = 20;
            end
            if (nargin < 1)
                obj.D = 2;
            end
            obj.xbounds = repmat([-32.768, 32.768], obj.D, 1);
        end

        function y = do_eval(obj, xx)
            if size(xx,1)~=obj.D
                error('Problem with input size')
            end

            sum1 = sum(xx.^2);
            sum2 = sum(cos(obj.c*xx));

            term1 = -obj.a * exp(-obj.b*sqrt(sum1/obj.D));
            term2 = -exp(sum2/obj.D);

            y = term1 + term2 + obj.a + exp(1);
            y(xx > obj.xbounds(:,2) | xx <  obj.xbounds(:,1)) = NaN;
            if strcmp(obj.opt, 'max')
                y = -y;
            end
        end
    end
end
