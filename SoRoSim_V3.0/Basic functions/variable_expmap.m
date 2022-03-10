function g = variable_expmap(theta,Gamma)

Gammahat   = dinamico_hat(Gamma);

if (theta<=1e-7)
    g      = diag([1 1 1 1])+Gammahat;
else
    Gammahatp2 = Gammahat*Gammahat;
    Gammahatp3 = Gammahatp2*Gammahat;
    tp2        = theta*theta;
    tp3        = tp2*theta;
    g          = diag([1 1 1 1])+Gammahat+...
                 ((1-cos(theta))/(tp2))*Gammahatp2+...
                 ((theta-sin(theta))/(tp3))*Gammahatp3;
end

% eofs