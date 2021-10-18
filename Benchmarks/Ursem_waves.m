classdef Ursem_waves
    properties
        D = 2
        xbounds = [-0.9, 1.2; -1.2, 1.2];
        name = 'Ursem Waves';
        opt = 'max';
        mean
        var
        takelog
        rescaling
        
    end
    methods
        function obj = Ursem_waves(rescaling)
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
            y=-(0.3*x1).^3+(x2.^2-4.5*x2.^2).*x1.*x2+4.7*cos(3*x1-(x2.^2).*(2+x1)).*sin(2.5*pi*x1);
            
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
% clear
% clc
% warning off
%
% x1min=-0.9;
% x1max=1.2;
% x2min=-1.2;
% x2max=1.2;
% R=1500; % steps resolution
% x1=x1min:(x1max-x1min)/R:x1max;
% x2=x2min:(x2max-x2min)/R:x2max;
%
% for j=1:length(x1)
%
%     for i=1:length(x2)
%         f(i)=-(0.3*x1(j)).^3+(x2(i).^2-4.5*x2(i).^2)*x1(j)*x2(i)+4.7*cos(3*x1(j)-(x2(i).^2)*(2+x1(j)))*sin(2.5*pi*x1(j));
%     end
%
%     f_tot(j,:)=f;
%
% end
%
% figure(1)
% meshc(x1,x2,f_tot);colorbar;set(gca,'FontSize',12);
% xlabel('x_2','FontName','Times','FontSize',20,'FontAngle','italic');
% set(get(gca,'xlabel'),'rotation',25,'VerticalAlignment','bottom');
% ylabel('x_1','FontName','Times','FontSize',20,'FontAngle','italic');
% set(get(gca,'ylabel'),'rotation',-25,'VerticalAlignment','bottom');
% zlabel('f(X)','FontName','Times','FontSize',20,'FontAngle','italic');
% title('3D View','FontName','Times','FontSize',24,'FontWeight','bold');
%
% figure(2)
% mesh(x1,x2,f_tot);view(0,90);colorbar;set(gca,'FontSize',12);
% xlabel('x_2','FontName','Times','FontSize',20,'FontAngle','italic');
% ylabel('x_1','FontName','Times','FontSize',20,'FontAngle','italic');
% zlabel('f(X)','FontName','Times','FontSize',20,'FontAngle','italic');
% title('X-Y Plane View','FontName','Times','FontSize',24,'FontWeight','bold');
%
% figure(3)
% mesh(x1,x2,f_tot);view(90,0);colorbar;set(gca,'FontSize',12);
% xlabel('x_2','FontName','Times','FontSize',20,'FontAngle','italic');
% ylabel('x_1','FontName','Times','FontSize',20,'FontAngle','italic');
% zlabel('f(X)','FontName','Times','FontSize',20,'FontAngle','italic');
% title('X-Z Plane View','FontName','Times','FontSize',24,'FontWeight','bold');
%
% figure(4)
% mesh(x1,x2,f_tot);view(0,0);colorbar;set(gca,'FontSize',12);
% xlabel('x_2','FontName','Times','FontSize',20,'FontAngle','italic');
% ylabel('x_1','FontName','Times','FontSize',20,'FontAngle','italic');
% zlabel('f(X)','FontName','Times','FontSize',20,'FontAngle','italic');
% title('Y-Z Plane View','FontName','Times','FontSize',24,'FontWeight','bold');