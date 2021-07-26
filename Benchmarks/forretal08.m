%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FORRESTER ET AL. (2008) FUNCTION
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


classdef forretal08
    properties
        D = 1
        xbounds = [0,1];
        name = 'Forrester et al (2008)';
        opt = 'max';
        mean
        var
        takelog
        rescaling
        
    end
    methods
        function obj = forretal(rescaling)
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
            fact1 = (6*xx - 2).^2;
            fact2 = sin(12*xx - 4);
            
            
            y = fact1 .* fact2;
            
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

