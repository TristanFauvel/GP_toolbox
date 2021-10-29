%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% TRID FUNCTION
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
%The function is usually evaluated on the hypercube xi ∈ [-d2, d2], for all i = 1, …, d.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef trid
    properties
        xbounds
        D
        name = 'Trid';
        opt = 'max';
                mean
        var
        takelog
        rescaling

    end
    methods
        function obj = trid(rescaling, D)
            if nargin ==0
                D = 2;
                rescaling = 0;
            elseif nargin == 1
                D = 2;
            end
            obj.rescaling = rescaling;
            obj.D = D;
            obj.xbounds = repmat([-D^2, D^2], D, 1);
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
            sum1 = (xx(1,:)-1).^2;
            sum2 = 0;
            
            for ii = 2:obj.D
                xi = xx(ii,:);
                xold = xx(ii-1,:);
                sum1 = sum1 + (xi-1).^2;
                sum2 = sum2 + xi.*xold;
            end
            
            y = sum1 - sum2;
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
        end
    end
end
