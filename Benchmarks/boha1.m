%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BOHACHEVSKY FUNCTION 1
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
%The functions are usually evaluated on the square xi ∈ [-100, 100], for all i = 1, 2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef boha1
    properties
        D = 2
        xbounds = [-100, 100;-100, 100];
        name = 'Bohachevsky';
                opt = 'max';
        mean
        var
        takelog
        rescaling

    end
    methods
        function obj = boha1(rescaling)
            if nargin<1
                obj.rescaling = 0;
            else
                obj.rescaling =rescaling;
            end
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
            x1 = xx(1,:);
            x2 = xx(2,:);

            term1 = x1.^2;
            term2 = 2*x2.^2;
            term3 = -0.3 * cos(3*pi*x1);
            term4 = -0.4 * cos(4*pi*x2);

            y = term1 + term2 + term3 + term4 + 0.7;
            y(xx > obj.xbounds(:,2) | xx <  obj.xbounds(:,1)) = NaN;
            
           if obj.rescaling
                if obj.takelog
                    y = log(y);
                end
                y = (y- obj.mean)./sqrt(obj.var);
            end
            
            if strcmp(obj.opt, 'max')
                y = -y;
            end
        end
    end
end
