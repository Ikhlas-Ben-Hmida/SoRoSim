function g = joint_expmap(xi)
x = 1;

theta      = norm(xi(1:3));
xihat      = dinamico_hat(xi);
if theta == 0
    g      = diag([1 1 1 1])+x*xihat;
else
    g      = diag([1 1 1 1])+x*xihat+...
            ((1-cos(x*theta))/(theta^2))*xihat^2+...
            ((x*theta-sin(x*theta))/(theta^3))*xihat^3;
end


% eof



