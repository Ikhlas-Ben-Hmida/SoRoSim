function coAdj = dinamico_coAdjoint(rototras)

coAdj          =[rototras(1:3,1:3) dinamico_tilde(rototras(1:3,4))*rototras(1:3,1:3);zeros(3,3) rototras(1:3,1:3)];
% coAdj(1:3,1:3)=rototras(1:3,1:3);
% coAdj(1:3,4:6)=dinamico_tilde(rototras(1:3,4))*rototras(1:3,1:3);
% coAdj(4:6,4:6)=rototras(1:3,1:3);

% eof