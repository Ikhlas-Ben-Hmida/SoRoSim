function Texpg = variable_Texpmap(h,theta,Gamma)

adjGamma       = dinamico_adj(Gamma);

if (theta<=1e-7)
    Texpg      = h*[1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1]+((h^2)/2)*adjGamma;
else
    
    tp2 = theta*theta;
    tp3 = tp2*theta;
    tp4 = tp3*theta;
    tp5 = tp4*theta;
    
    adjGammap2 = adjGamma*adjGamma;
    adjGammap3 = adjGammap2*adjGamma;
    adjGammap4 = adjGammap3*adjGamma;
    
    sin_htheta = sin(h*theta);
    cos_htheta = cos(h*theta);
    
    t1 = h*theta*sin_htheta;
    t2 = h*theta*cos_htheta;
    
    Texpg      = h*[1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1]+((4-4*cos_htheta-t1)/(2*tp2))*adjGamma+((4*h*theta-5*sin_htheta+t2)/(2*tp3))*adjGammap2+((2-2*cos_htheta-t1)/(2*tp4))*adjGammap3+((2*h*theta-3*sin_htheta+t2)/(2*tp5))*adjGammap4;
end

% eof