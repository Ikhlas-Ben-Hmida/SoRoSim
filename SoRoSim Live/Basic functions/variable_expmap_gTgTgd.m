function [g,Tg,Tgd]=variable_expmap_gTgTgd(Gamma,Gammad)

k      = Gamma(1:3);
theta  = norm(k);
kd     = Gammad(1:3);
thetad = (kd'*k)/theta;

Gammahat  = dinamico_hat(Gamma);
adjGamma  = dinamico_adj(Gamma);
adjGammad = dinamico_adj(Gammad);

if (theta<=1e-6)
    g   = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]+Gammahat;
    Tg  = [1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1]+0.5*adjGamma;
    Tgd = 0.5*adjGammad;
else
    Gammahatp2 = Gammahat*Gammahat;
    Gammahatp3 = Gammahatp2*Gammahat;
    
    adjGammap2 = adjGamma*adjGamma;
    adjGammap3 = adjGammap2*adjGamma;
    adjGammap4 = adjGammap3*adjGamma;
    
    adjGammad2 = adjGammad*adjGamma+adjGamma*adjGammad;
    adjGammad3 = adjGammad2*adjGamma+adjGammap2*adjGammad;
    adjGammad4 = adjGammad3*adjGamma+adjGammap3*adjGammad;
    
    tp2        = theta*theta;
    tp3        = tp2*theta;
    tp4        = tp3*theta;
    tp5        = tp4*theta;
    tp6        = tp5*theta;
    
    sintheta   = sin(theta);
    costheta   = cos(theta);
    
    t1 = theta*sintheta;
    t2 = theta*costheta;
    t3 = -8+(8-tp2)*costheta+5*t1;
    t4 = -8*theta+(15-tp2)*sintheta-7*t2;
    t5 = (4-4*costheta-t1)/(2*tp2);
    t6 = (4*theta-5*sintheta+t2)/(2*tp3);
    t7 = (2-2*costheta-t1)/(2*tp4);
    t8 = (2*theta-3*sintheta+t2)/(2*tp5);
    
    g   = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]+Gammahat+...
          (1-costheta)/(tp2)*Gammahatp2+...
          ((theta-sintheta)/(tp3))*Gammahatp3;
    Tg  = [1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1]+t5*adjGamma+...
          t6*adjGammap2+...
          t7*adjGammap3+...
          t8*adjGammap4;
    Tgd = (thetad*(t3)/(2*tp3))*adjGamma+t5*adjGammad+...
          (thetad*(t4)/(2*tp4))*adjGammap2+t6*adjGammad2+...
          (thetad*(t3)/(2*tp5))*adjGammap3+t7*adjGammad3+...
          (thetad*(t4)/(2*tp6))*adjGammap4+t8*adjGammad4;
end
