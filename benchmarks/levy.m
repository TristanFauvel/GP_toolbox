%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEVY FUNCTION
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef levy
    properties
        name = 'Levy';
        xbounds
        D
        opt = 'max'; 
        mean
        var
        takelog
        rescaling

    end
    methods
        function obj = levy(rescaling, D)
            if nargin <1
                rescaling = 0;
                D =2;
            elseif nargin <2
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
            n=size(xx,2);
            w = NaN(obj.D, n);
            for ii = 1:obj.D
                w(ii,:) = 1 + (xx(ii,:) - 1)/4;
            end
            
            term1 = (sin(pi*w(1,:))).^2;
            term3 = (w(obj.D,:)-1).^2 .* (1+(sin(2*pi*w(obj.D,:))).^2);
            
            sum = 0;
            for ii = 1:(obj.D-1)
                wi = w(ii,:);
                new = (wi-1).^2 .* (1+10*(sin(pi*wi+1)).^2);
                sum = sum + new;
            end
            
            y = term1 + sum + term3;
            
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