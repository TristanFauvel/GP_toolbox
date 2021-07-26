%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SUM SQUARES FUNCTION
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
% xx = [x1, x2, ..., xd]
%The function is usually evaluated on the hypercube xi ∈ [-10, 10], for all i = 1, …, d, although this may be restricted to the hypercube xi ∈ [-5.12, 5.12], for all i = 1, …, d.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef sumsqu
    properties
        xbounds
        D
        name = 'Sum Squares';
        opt = 'max';
                mean
        var
        takelog
        rescaling

    end
    methods
        function obj = sumsqu(rescaling, D)
           if nargin ==0
                D = 2;
                rescaling = 0;
            elseif nargin == 1
                D = 2;
            end
            obj.rescaling = rescaling;
            obj.D = D;
            obj.xbounds = repmat([-10, 10], D, 1);
            if obj.rescaling
                load('benchmarks_rescaling.mat', 't');
                obj.var = t(t.Names == obj.name,:).Variance;
                obj.mean = t(t.Names == obj.name,:).Mean;
                obj.takelog = t(t.Names == obj.name,:).TakeLog;
            end
        end
        function y = do_eval(obj, xx)
            if size(xx,1)~=obj.D
                error('Problem with input size')
            end
            d = length(xx);
            sum = 0;
            for ii = 1:obj.D
                xi = xx(ii,:);
                sum = sum + ii*xi.^2;
            end            
            y = sum;     
            if obj.rescaling
                if obj.takelog
                    y = log(y);
                end
                y = (y- obj.mean)./sqrt(obj.var);
            end
            if strcmp(obj.opt, 'max')
                y = -y;
            end
            y(xx > obj.xbounds(:,2) | xx <  obj.xbounds(:,1)) = NaN;
        end
    end
end

