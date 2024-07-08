%Use this function to evaluate the kinetic, gravitational potential and
%elastic potential energy of the beam. Solve dynamics before runing this
%code
% load('FreeFall.mat')
load('DynamicsSolution.mat')

ndof   = T1.ndof;
nip    = T1.CVTwists{1}(2).nip;
lpf    = T1.VLinks(1).lp{1};
Ws     = T1.CVTwists{1}(2).Ws;
Xs     = T1.CVTwists{1}(2).Xs;
tt     = 0:0.01:t(end);
nt     = length(tt);
T      = zeros(nt,1);
Vg     = zeros(nt,1);
Ve     = zeros(nt,1);
lp     = T1.VLinks(1).lp{1};

Ms       = T1.CVTwists{1}(2).Ms;
Es       = T1.CVTwists{1}(2).Es;
Bdof     = T1.CVTwists{1}(2).Bdof;
Bodr     = T1.CVTwists{1}(2).Bodr;
xi_star  = T1.CVTwists{1}(2).xi_star;

for i=1:nt

    qqdt  = interp1(t,qqd,tt(i));
    eta   = T1.ScrewVelocity(qqdt(1:ndof)',qqdt(ndof+1:2*ndof)');
    g     = T1.FwdKinematics(qqdt(1:ndof)');
    i_sig = 3;

    for ii=2:nip-1
        
        eta_here = eta((i_sig-1)*6+1:i_sig*6);
        g_here   = g((i_sig-1)*4+1:i_sig*4,:);
        M_here   = Ms((ii-1)*6+1:ii*6,:);
        E_here   = Es((ii-1)*6+1:ii*6,:);
        B_here   = zeros(6,T1.ndof); 
        x        = Xs(ii)*lpf;

        for jj = 1:6
            for k = 1:Bdof(jj)*Bodr(jj)+Bdof(jj)
                kk              = sum(Bdof(1:jj-1).*Bodr(1:jj-1))+sum(Bdof(1:jj-1))+k;
                B_here(jj,kk)   = x^(k-1);
            end
        end
        
        xi_starhere = xi_star(6*(ii-1)+1:6*(ii),1);
        xi_here     = B_here*qqdt(1:ndof)'+xi_starhere;

        T(i)  = T(i)+Ws(ii)*lp*0.5*eta_here'*M_here*eta_here;
        Vg(i) = Vg(i)+Ws(ii)*lp*M_here(6,6)*-T1.G(6)*g_here(3,4);
        Ve(i) = Ve(i)+Ws(ii)*lp*0.5*(xi_here-xi_starhere)'*E_here*(xi_here-xi_starhere);
        i_sig = i_sig+1;
    end

end
close all
plot(tt,T,'r','LineWidth',1.5)
hold on
plot(tt,Vg,'g','LineWidth',1.5)
plot(tt,Ve,'b','LineWidth',1.5)
E=T+Vg+Ve;
plot(tt,E,'k','LineWidth',1.5)
legend('Kinetic energy','Gravitational potential energy','Elastic potential energy','Total energy')
xlabel('Time in s')
ylabel('Energy in J')
set(gca,'FontSize',18)