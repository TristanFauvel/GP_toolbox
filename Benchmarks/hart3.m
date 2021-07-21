%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%

% HARTMANN 3-DIMENSIONAL FUNCTION

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

% xx = [x1, x2, x3]

%The function is usually evaluated on the hypercube xi ∈ (0, 1), for all i = 1, 2, 3.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef hart3
    properties
        D = 3;
        xbounds = [0,1;0,1;0,1];
        name = 'Hartmann 3-D';
                opt = 'max';

    end
    methods
        function y = do_eval(obj, xx)
            if size(xx,1)~=obj.D
                error('Problem with input size')
            end
            alpha = [1.0, 1.2, 3.0, 3.2]';
            
            A = [3.0, 10, 30;
                
            0.1, 10, 35;
            
            3.0, 10, 30;
            
            0.1, 10, 35];
        
        P = 10^(-4) * [3689, 1170, 2673;
            
        4699, 4387, 7470;
        
        1091, 8732, 5547;
        
        381, 5743, 8828];
    
    
    outer = 0;
    
    for ii = 1:4
        
        inner = 0;
        
        for jj = 1:3
            
            xj = xx(jj,:);
            
            Aij = A(ii, jj);
            
            Pij = P(ii, jj);
            
            inner = inner + Aij*(xj-Pij).^2;
            
        end
        
        new = alpha(ii) * exp(-inner);
        
        outer = outer + new;
        
    end
    
    
    y = -outer;
    if strcmp(obj.opt, 'max')
                y = -y;
            end
    
    y(xx > obj.xbounds(:,2) | xx <  obj.xbounds(:,1)) = NaN;
        end
    end
end
