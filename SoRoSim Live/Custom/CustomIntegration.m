function [nip,Xs,Ws]=CustomIntegration(Tr,i,j,t,q)
% Edit this if you want to change the integration scheme during the
% analysis. Make sure to enable the property CI (make it true (1)) and the property 
% CIFn should be the function handle for custom integration

r02 = log(50)^1/2/q(2); %distance at which value of basis becomes 0.02 (2%)
X1 = t/5-r02;
X2 = t/5+r02;
if (X1<=0&&X2>=1)||X1>=1||X2<=0
    [Xs,Ws,nip]=GaussQuadrature(5);
elseif X1>0&&X2<1
    
    X0 = 0;
    %from 0 to X1
    [Xs1,Ws1,nip1]=GaussQuadrature(3);
    Ws1 = Ws1*X1;
    Xs1 = X0+Xs1*X1;
    X0  = X1;
    %from X1 to X2
    [Xs2,Ws2,nip2]=GaussQuadrature(5);
    Ws2 = Ws2*(X2-X1);
    Xs2 = X0+Xs2*(X2-X1);
    X0  = X2;
    %from X2 to 1
    [Xs3,Ws3,nip3]=GaussQuadrature(3);
    Ws3 = Ws3*(1-X2);
    Xs3 = X0+Xs3*(1-X2);
    
    nip = nip1+nip2+nip3-2;
    Xs  = [Xs1;Xs2(2:end);Xs3(2:end)];
    Ws  = [Ws1;Ws2(2:end);Ws3(2:end)];
    
elseif X1>0
     X0 = 0;
    %from 0 to X1
    [Xs1,Ws1,nip1]=GaussQuadrature(3);
    Ws1 = Ws1*X1;
    Xs1 = X0+Xs1*X1;
    X0  = X1;
    %from X1 to 1
    [Xs2,Ws2,nip2]=GaussQuadrature(5);
    Ws2 = Ws2*(1-X1);
    Xs2 = X0+Xs2*(1-X1);
    
    nip = nip1+nip2-1;
    Xs  = [Xs1;Xs2(2:end)];
    Ws  = [Ws1;Ws2(2:end)];
else %X2<1
     X0 = 0;

    %from 0 to X2
    [Xs2,Ws2,nip2]=GaussQuadrature(5);
    Ws2 = Ws2*X2;
    Xs2 = X0+Xs2*X2;
    X0  = X2;
    %from X2 to 1
    [Xs3,Ws3,nip3]=GaussQuadrature(3);
    Ws3 = Ws3*(1-X2);
    Xs3 = X0+Xs3*(1-X2);
    
    nip = nip2+nip3-1;
    Xs  = [Xs2;Xs3(2:end)];
    Ws  = [Ws2;Ws3(2:end)];
end
end