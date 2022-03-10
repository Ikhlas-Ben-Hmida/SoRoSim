function se3 = dinamico_hat(screw)

se3          = [dinamico_tilde(screw(1:3)) screw(4:6);zeros(1,4)];
% se3(1:3,1:3)=dinamico_tilde(screw(1:3));
% se3(1:3,4)  =screw(4:6);

% eof