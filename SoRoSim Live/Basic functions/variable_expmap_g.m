function g=variable_expmap_g(Gamma)

k      = Gamma(1:3);
theta  = norm(k);

Gammahat  = dinamico_hat(Gamma);

if (theta<=1e-6)
    g   = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]+Gammahat;
else
    Gammahatp2 = Gammahat*Gammahat;
    Gammahatp3 = Gammahatp2*Gammahat;
    
    tp2        = theta*theta;
    tp3        = tp2*theta;
    
    sintheta   = sin(theta);
    costheta   = cos(theta);
    
    g   = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]+Gammahat+...
          (1-costheta)/(tp2)*Gammahatp2+...
          ((theta-sintheta)/(tp3))*Gammahatp3;
end
