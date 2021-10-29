%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GOLDSTEIN-PRICE FUNCTION
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef goldpr
    properties
        D = 2;
        xbounds = [-2, 2; -2, 2];
        name = 'Goldstein-Price';
        opt = 'max';
        mean
        var
        takelog
        rescaling
        
    end
    methods
        function obj = goldpr(rescaling)
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
            
            fact1a = (x1 + x2 + 1).^2;
            fact1b = 19 - 14*x1 + 3*x1.^2 - 14*x2 + 6*x1.*x2 + 3*x2.^2;
            fact1 = 1 + fact1a.*fact1b;
            
            fact2a = (2*x1 - 3*x2).^2;
            fact2b = 18 - 32*x1 + 12*x1.^2 + 48*x2 - 36*x1.*x2 + 27*x2.^2;
            fact2 = 30 + fact2a.*fact2b;
            
            y = fact1.*fact2;
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
