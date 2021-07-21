%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PERM FUNCTION 0, d, beta
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
% xx = [x1, x2]
% b = constant (optional), with default value 10
%The function is usually evaluated on the hypercube xi ∈ [-d, d], for all i = 1, …, d.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef perm0db
    properties
        xbounds
        D
        name = 'Perm 0,d,$\beta$';
        b = 10;
        opt = 'max';
    end
    methods
        function obj = perm0db(D,b)            
            if nargin ==0
                D = 2;
            elseif nargin ==1
                obj.b = b;
            end
            obj.D = D;
            obj.xbounds = repmat([-D, D], D, 1);
        end
        function y = do_eval(obj, xx)
            if size(xx,1)~=obj.D
                error('Problem with input size')
            end
            outer = 0;
            
            for ii = 1:obj.D
                inner = 0;
                for jj = 1:obj.D
                    xj = xx(jj,:);
                    inner = inner + (jj+obj.b)*(xj.^ii-(1/jj)^ii);
                end
                outer = outer + inner.^2;
            end
            
            y = outer;    
            
            if strcmp(obj.opt, 'max')
                y = -y;
            end
            y(xx > obj.xbounds(:,2) | xx <  obj.xbounds(:,1)) = NaN;
        end
    end
end

