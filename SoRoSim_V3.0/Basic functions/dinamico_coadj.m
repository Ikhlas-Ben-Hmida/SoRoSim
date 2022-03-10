function coadj = dinamico_coadj(screw)
t1             = dinamico_tilde(screw(1:3));
coadj          = [t1 dinamico_tilde(screw(4:6));zeros(3,3) t1];
% coadj(1:3,1:3)=dinamico_tilde(screw(1:3));
% coadj(1:3,4:6)=dinamico_tilde(screw(4:6));
% coadj(4:6,4:6)=coadj(1:3,1:3);

% eof