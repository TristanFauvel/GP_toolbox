%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LANGERMANN FUNCTION
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
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS:
%
% xx = [x1, x2, ..., xd]
% m  = constant (optional), with default value 5
% c  = m-dimensional vector (optional), with default value [1, 2, 5, 2, 3]
%      (when m=5)
% A  = (mxd)-dimensional matrix (optional), with default value
%      [3, 5; 5, 2; 2, 1; 1, 4; 7, 9] (when m=5 and d=2)
%The function is usually evaluated on the hypercube xi ∈ [0, 10], for all i = 1, …, d.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef langer
    properties
        xbounds
        D
        name = 'Langer';
        m = 5;
        c = [1, 2, 5, 2, 3];
        A = [3, 5; 5, 2; 2, 1; 1, 4; 7, 9];
        opt = 'max';
                mean
        var
        takelog
        rescaling

    end
    methods
        function obj = langer(rescaling, D, m, c, A)
            if nargin == 0
                rescaling = 0;
                D = 2;
            elseif nargin ==1
                D = 2;           
            elseif nargin == 3
                obj.m = m;
            elseif nargin == 4
                obj.m = m;
                obj.c = c;
            elseif nargin == 5
                obj.m = m;
                obj.c = c;
                obj.A = A;
            end      
            obj.rescaling = rescaling;
            obj.D = D;
            obj.xbounds = repmat([0,10], D, 1);
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
            outer = 0;
            for ii = 1:obj.m
                inner = 0;
                for jj = 1:obj.D
                    xj = xx(jj,:);
                    Aij = obj.A(ii,jj);
                    inner = inner + (xj-Aij).^2;
                end
                new = obj.c(ii) * exp(-inner/pi) .* cos(pi*inner);
                outer = outer + new;
            end
            y = outer; 
            
            
            if obj.rescaling
                if obj.takelog
                     if any(y<=0)
                        error('Log of negative value')
                    end
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

