%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ROSENBROCK FUNCTION
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
%The function is usually evaluated on the hypercube xi ∈ [-5, 10], for all i = 1, …, d,
%although it may be restricted to the hypercube xi ∈ [-2.048, 2.048], for all i = 1, …, d.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
% xx = [x1, x2, ..., xd]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef rosen
    properties
        xbounds
        D
        name = 'Rosenbrock';
        opt = 'max';
    end
    methods
        function obj = rosen(D)
            if nargin == 0
                D = 2;
            end
            obj.D = D;
            obj.xbounds = repmat([-2.048, 2.048], D, 1);
          
        end
        function y = do_eval(obj, xx)
            if size(xx,1)~=obj.D
                error('Problem with input size')
            end
            sum = 0;
            for ii = 1:(obj.D-1)
                xi = xx(ii,:);
                xnext = xx(ii+1,:);
                new = 100*(xnext-xi.^2).^2 + (xi-1).^2;
                sum = sum + new;
            end
            y = sum;    
            
            if strcmp(obj.opt, 'max')
                y = -y;
            end
            y(xx > obj.xbounds(:,2) | xx <  obj.xbounds(:,1)) = NaN;
        end
    end
end
