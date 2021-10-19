function dotTexpGamma = variable_dotTexpmap(h,theta,dottheta,Gamma,dotGamma)

adj_dotGamma          = dinamico_adj(dotGamma);

if (theta<=1e-7)
    dotTexpGamma = ((h^2)/2)*adj_dotGamma;
else
    adj_Gamma     = dinamico_adj(Gamma);
    adj_Gammap2   = adj_Gamma*adj_Gamma;
    adj_Gammap3   = adj_Gammap2*adj_Gamma;
    adj_Gammap4   = adj_Gammap3*adj_Gamma;
    
    adj_dot2Gamma = adj_dotGamma*adj_Gamma+adj_Gamma*adj_dotGamma;
    adj_dot3Gamma = adj_dot2Gamma*adj_Gamma+adj_Gammap2*adj_dotGamma;
    adj_dot4Gamma = adj_dot3Gamma*adj_Gamma+adj_Gammap3*adj_dotGamma;
    
    tp2 = theta*theta;
    tp3 = tp2*theta;
    tp4 = tp3*theta;
    tp5 = tp4*theta;
    
    sin_htheta = sin(h*theta);
    cos_htheta = cos(h*theta);
    
    t1 = h*theta*sin_htheta;
    t2 = h*theta*cos_htheta;
    t3 = -8+(8-h^2*tp2)*cos_htheta+5*t1;
    t4 = -8*h*theta+(15-h^2*tp2)*sin_htheta-7*t2;
    
    
    dotTexpGamma = (dottheta*(t3)/(2*tp3))*adj_Gamma+...
                   ((4-4*cos_htheta-t1)/(2*tp2))*adj_dotGamma+...
                   (dottheta*(t4)/(2*tp4))*adj_Gammap2+...
                   ((4*h*theta-5*sin_htheta+t2)/(2*tp3))*adj_dot2Gamma+...
                   (dottheta*(t3)/(2*tp5))*adj_Gammap3+...
                   ((2-2*cos_htheta-t1)/(2*tp4))*adj_dot3Gamma+...
                   (dottheta*(t4)/(2*theta^6))*adj_Gammap4+...
                   ((2*h*theta-3*sin_htheta+t2)/(2*tp5))*adj_dot4Gamma;
end

% eof