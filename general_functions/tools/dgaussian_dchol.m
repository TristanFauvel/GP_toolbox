
dGdtheta = zeros(1, 3);
detS  = det(S);

dmahadL= -2*invS*x*x'*invS*L;

for i = 1:3
    if i ==1
        dLdtheta_i = [1,0;0,0];
        ddetSdtheta_i = 2*theta(1)*theta(3)^2;
    elseif i ==2
        dLdtheta_i = [0,0;1,0];
        ddetSdtheta_i = 0;        
    elseif i==3
        dLdtheta_i = [0,0;0,1];
        ddetSdtheta_i = 2*theta(1)^2*theta(3);        
    end
    dmahadtheta_i = dmahadL*dLdtheta_i;
    
    ddetSdtheta(i)=ddetSdtheta_i;
    dmahadtheta(i) = dmahadL*dLdtheta_i;
    dGtheta(i) = -1/detS*G*ddetSdtheta(i)-0.5*G*dmahadtheta(i);
end