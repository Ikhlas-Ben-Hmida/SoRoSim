%Class that assigns DoFs and Bases to link pieces and can calculate the
%twist given joint angles 
%Last modified by Anup Teejo Mathew 14.06.2022

classdef SorosimTwist < handle
    
    properties
        Type         %Base type (Monomial, Lagrange Polynomial, Linear Interpolation, Gaussian, Custom, Non-linear Gaussian)
        SubClass     %Now only for FEM Like basis (linear,quadratic,cubic)
        Bdof         %(6x1) array specifying the allowable DoFs of a soft piece. 1 if allowed 0 if not.
        Bodr         %(6x1) array specifying the order of allowed DoF (0: constant, 1: linear, 2: quadratic,...)
        dof          %degress of freedom of each base
        
        Bh           %Function handle for base
        B            %(6xdof) Base matrix calculated at lumped joints or ((6xnGauss)xdof) base matrices computed at every significant points of a soft division
        B_Z1         %Base calculated at 4th order first Zanna point (Xs+Z1*(delta(Xs)))
        B_Z2         %Base calculated at 4th order second Zanna point (Xs+Z2*(delta(Xs)))
        B_Z          %Base calculated at 2nd order Zanna point 
        
        xi_starfn    %Reference strain vector as a function of X
        xi_star      %(6x1) reference strain vector at the lumped joint or ((6xnGauss)x4) reference strain vectors computed at Gauss quadrature and Zannah collocation points
        
        Link         %Link associated with this twist only for soft link
        div          %Division associated with this twist
        nip           %number of integration point including boundaries
        Xs           %integration points 
        Ws           %weights of integration point
        Ms           %Inertia matrix of cross-section (6nip x 6) matrix
        Es           %Stiffness matrix (6nip x 6) matrix
        Gs           %Damping matrix (6nip x 6) matrix
        
        Xadd         %additional integration points (nx1) vector 
        CI           %logical 0 by default 1 if custom integration is enabled
        CIFn         %function handle for custom integration
    end
    
    methods
        
        function T = SorosimTwist(varargin)
            
            if nargin==3
                i = varargin{1};
                j = varargin{2};
                Link = varargin{3};
                
                T.Link = Link;
                T.div = j;
                
                Z1     = 1/2-sqrt(3)/6;          %Zanna quadrature coefficient 4th order
                Z2     = 1/2+sqrt(3)/6;          %Zanna quadrature coefficient 4th order
                Z      = 1/2;                    %Zanna quadrature coefficient 2nd order
               
                basedef(i,j);
                load('Base_properties.mat','Type','Bdof','Bodr','xi_stars','SubClass')
                
                T.Type = Type;
                T.Bdof = Bdof;
                T.Bodr = Bodr;
                
                if any(strcmp(Type,{'Monomial','Legendre Polynomial','Chebychev','Fourier'}))
                    if strcmp(Type,'Fourier')
                        nGauss = max(Bodr)*4;
                    else
                        nGauss = max(Bodr)+1;
                    end
                    nGauss(nGauss<5)=5;
                    prompt           = ['Number of Gaussian Quadrature points (min ',num2str(nGauss),'):'];
                    dlgtitle         = 'Gaussian Points ';
                    definput         = {num2str(nGauss)};
                    opts.Interpreter = 'tex';
                    answer2 = inputdlg(prompt,dlgtitle,[1 70],definput,opts);
                    nGauss = str2num(answer2{1});
                    [Xs,Ws,nip]=GaussQuadrature(nGauss);
                end
                    
                switch Type
                    case 'Monomial'
                        
                        dof  = sum(Bdof.*(Bodr+1));
                        
                        B    = zeros(nip*6,dof);
                        B_Z1 = zeros(nip*6,dof);
                        B_Z2 = zeros(nip*6,dof);
                        B_Z  = zeros(nip*6,dof);

                        X = Xs(1);
                        B(1:6,:) = B_Monomial(X,Bdof,Bodr);
                        
                        for ii=2:nip
                            X = Xs(ii);
                            B(6*(ii-1)+1:6*ii,:) = B_Monomial(X,Bdof,Bodr);
                            X = Xs(ii-1)+Z1*(Xs(ii)-Xs(ii-1));
                            B_Z1(6*(ii-2)+1:6*(ii-1),:) = B_Monomial(X,Bdof,Bodr);
                            X = Xs(ii-1)+Z2*(Xs(ii)-Xs(ii-1));
                            B_Z2(6*(ii-2)+1:6*(ii-1),:) = B_Monomial(X,Bdof,Bodr);
                            X = Xs(ii-1)+Z*(Xs(ii)-Xs(ii-1));
                            B_Z(6*(ii-2)+1:6*(ii-1),:)  = B_Monomial(X,Bdof,Bodr);                                
                        end
                        
                        file = 'B_Monomial';
                        Bh   = str2func(['@(X,Bdof,Bodr)',file,'(X,Bdof,Bodr)']);
                        
                    case 'Legendre Polynomial'
                        
                        dof      = sum(Bdof.*(Bodr+1));
                        
                        B        = zeros(nip*6,dof);
                        B_Z1     = zeros(nip*6,dof);
                        B_Z2     = zeros(nip*6,dof);
                        B_Z      = zeros(nip*6,dof);

                        X = Xs(1);
                        B(1:6,:) = B_LegendrePolynomial(X,Bdof,Bodr);
                        
                        for ii=2:nip
                            X = Xs(ii);
                            B(6*(ii-1)+1:6*ii,:) = B_LegendrePolynomial(X,Bdof,Bodr);
                            X = Xs(ii-1)+Z1*(Xs(ii)-Xs(ii-1));
                            B_Z1(6*(ii-2)+1:6*(ii-1),:) = B_LegendrePolynomial(X,Bdof,Bodr);
                            X = Xs(ii-1)+Z2*(Xs(ii)-Xs(ii-1));
                            B_Z2(6*(ii-2)+1:6*(ii-1),:) = B_LegendrePolynomial(X,Bdof,Bodr);
                            X = Xs(ii-1)+Z*(Xs(ii)-Xs(ii-1));
                            B_Z(6*(ii-2)+1:6*(ii-1),:)  = B_LegendrePolynomial(X,Bdof,Bodr);                                
                        end
                        
                        file = 'B_LegendrePolynomial';
                        Bh   = str2func(['@(X,Bdof,Bodr)',file,'(X,Bdof,Bodr)']);
                        
                    case 'Chebychev'
                        
                        dof      = sum(Bdof.*(Bodr+1));
                        
                        B        = zeros(nip*6,dof);
                        B_Z1     = zeros(nip*6,dof);
                        B_Z2     = zeros(nip*6,dof);
                        B_Z      = zeros(nip*6,dof);

                        X = Xs(1);
                        B(1:6,:) = B_Chebychev(X,Bdof,Bodr);
                        
                        for ii=2:nip
                            X = Xs(ii);
                            B(6*(ii-1)+1:6*ii,:) = B_LegendrePolynomial(X,Bdof,Bodr);
                            X = Xs(ii-1)+Z1*(Xs(ii)-Xs(ii-1));
                            B_Z1(6*(ii-2)+1:6*(ii-1),:) = B_LegendrePolynomial(X,Bdof,Bodr);
                            X = Xs(ii-1)+Z2*(Xs(ii)-Xs(ii-1));
                            B_Z2(6*(ii-2)+1:6*(ii-1),:) = B_LegendrePolynomial(X,Bdof,Bodr);
                            X = Xs(ii-1)+Z*(Xs(ii)-Xs(ii-1));
                            B_Z(6*(ii-2)+1:6*(ii-1),:)  = B_LegendrePolynomial(X,Bdof,Bodr);                                
                        end
                        
                        file = 'B_Chebychev';
                        Bh   = str2func(['@(X,Bdof,Bodr)',file,'(X,Bdof,Bodr)']);
                        
                    case 'Fourier'
                        
                        dof  = sum(Bdof.*(2*Bodr+1));
                        
                        B    = zeros(nip*6,dof);
                        B_Z1 = zeros(nip*6,dof);
                        B_Z2 = zeros(nip*6,dof);
                        B_Z  = zeros(nip*6,dof);

                        X = Xs(1);
                        B(1:6,:) = B_Fourier(X,Bdof,Bodr);
                        
                        for ii=2:nip
                            X = Xs(ii);
                            B(6*(ii-1)+1:6*ii,:) = B_Fourier(X,Bdof,Bodr);
                            X = Xs(ii-1)+Z1*(Xs(ii)-Xs(ii-1));
                            B_Z1(6*(ii-2)+1:6*(ii-1),:) = B_Fourier(X,Bdof,Bodr);
                            X = Xs(ii-1)+Z2*(Xs(ii)-Xs(ii-1));
                            B_Z2(6*(ii-2)+1:6*(ii-1),:) = B_Fourier(X,Bdof,Bodr);
                            X = Xs(ii-1)+Z*(Xs(ii)-Xs(ii-1));
                            B_Z(6*(ii-2)+1:6*(ii-1),:)  = B_Fourier(X,Bdof,Bodr);                                
                        end
                        
                        file = 'B_Fourier';
                        Bh   = str2func(['@(X,Bdof,Bodr)',file,'(X,Bdof,Bodr)']);
                        
                    case 'FEM Like'
                        
                        n_ele = max(Bodr); % All numbers must be same max is used considering user error
                        nGauss = 5; %minimum
                        if strcmp(SubClass,'Linear')
                            dof = sum(Bdof*(n_ele+1));
                            if nGauss/n_ele-floor(nGauss/n_ele)==0
                                nGausse_min = floor(nGauss/n_ele);
                            else
                                nGausse_min = floor(nGauss/n_ele)+1;
                            end
                            nGausse_min(nGausse_min<2)=2;
                        elseif strcmp(SubClass,'Quadratic')
                            dof = sum(Bdof*(2*n_ele+1));
                            if nGauss/n_ele-floor(nGauss/n_ele)==0
                                nGausse_min = floor(nGauss/n_ele);
                            else
                                nGausse_min = floor(nGauss/n_ele)+1;
                            end
                            nGausse_min(nGausse_min<3)=3;
                        else % cubic
                            dof = sum(Bdof*(3*n_ele+1));
                            if nGauss/n_ele-floor(nGauss/n_ele)==0
                                nGausse_min = floor(nGauss/n_ele);
                            else
                                nGausse_min = floor(nGauss/n_ele)+1;
                            end
                            nGausse_min(nGausse_min<4)=4;
                        end
                        
                        prompt           = ['Gauss Quadrature points per element (min ', num2str(nGausse_min),'):'];
                        dlgtitle         = 'Gaussian Points ';
                        definput         = {num2str(nGausse_min)};
                        opts.Interpreter = 'tex';
                        answer2 = inputdlg(prompt,dlgtitle,[1 70],definput,opts);
                        nGausse = str2num(answer2{1});
                        
                        [Xse,Wse,npe]=GaussQuadrature(nGausse);
                        
                        nip = n_ele*(npe-1)+1;
                        Ws = zeros(nip,1);
                        Xs = zeros(nip,1);
                        Ws(1:end-1)=repmat(Wse(1:end-1),n_ele,1)/n_ele;
                        
                        X0 = 0;
                        for i=1:n_ele
                            Xs((i-1)*(npe-1)+1:i*(npe-1)) = X0+Xse(1:npe-1)/n_ele;
                            X0 = X0+1/n_ele;
                        end
                        Xs(end)=1;
                        
                        B        = zeros(nip*6,dof);
                        B_Z1     = zeros(nip*6,dof);
                        B_Z2     = zeros(nip*6,dof);
                        B_Z      = zeros(nip*6,dof);
                        
                        X = Xs(1);
                        B(1:6,:) = B_FEMLike(X,Bdof,Bodr,SubClass);
                        for ii=2:nip
                            X = Xs(ii);
                            B(6*(ii-1)+1:6*ii,:) = B_FEMLike(X,Bdof,Bodr,SubClass);
                            X = Xs(ii-1)+Z1*(Xs(ii)-Xs(ii-1));
                            B_Z1(6*(ii-2)+1:6*(ii-1),:) = B_FEMLike(X,Bdof,Bodr,SubClass);
                            X = Xs(ii-1)+Z2*(Xs(ii)-Xs(ii-1));
                            B_Z2(6*(ii-2)+1:6*(ii-1),:) = B_FEMLike(X,Bdof,Bodr,SubClass);
                            X = Xs(ii-1)+Z*(Xs(ii)-Xs(ii-1));
                            B_Z(6*(ii-2)+1:6*(ii-1),:)  = B_FEMLike(X,Bdof,Bodr,SubClass);                                
                        end
                        file = 'B_FEMLike';
                        Bh   = str2func(['@(X,Bdof,Bodr,SubClass)',file,'(X,Bdof,Bodr,SubClass)']);
                        T.SubClass = SubClass;
                    case 'Hermite Spline'
                        
                        n_ele = max(Bodr); % All numbers must be same max is used considering user error
                        nGauss = 5; %minimum
                        
                        dof      = sum(Bdof*(2*n_ele+2));
                        if nGauss/n_ele-floor(nGauss/n_ele)==0
                            nGausse_min = floor(nGauss/n_ele);
                        else
                            nGausse_min = floor(nGauss/n_ele)+1;
                        end
                        nGausse_min(nGausse_min<4)=4;
                        
                        prompt           = ['Gauss Quadrature points per element (min ', num2str(nGausse_min),'):'];
                        dlgtitle         = 'Gaussian Points ';
                        definput         = {num2str(nGausse_min)};
                        opts.Interpreter = 'tex';
                        answer2 = inputdlg(prompt,dlgtitle,[1 70],definput,opts);
                        nGausse = str2num(answer2{1});
                        
                        [Xse,Wse,npe]=GaussQuadrature(nGausse);
                        
                        nip = n_ele*(npe-1)+1;
                        Ws = zeros(nip,1);
                        Xs = zeros(nip,1);
                        Ws(1:end-1)=repmat(Wse(1:end-1),n_ele,1)/n_ele;
                        
                        X0 = 0;
                        for i=1:n_ele
                            Xs((i-1)*(npe-1)+1:i*(npe-1)) = X0+Xse(1:npe-1)/n_ele;
                            X0 = X0+1/n_ele;
                        end
                        Xs(end)=1;
                        
                        B    = zeros(nip*6,dof);
                        B_Z1 = zeros(nip*6,dof);
                        B_Z2 = zeros(nip*6,dof);
                        B_Z  = zeros(nip*6,dof);

                        ii = 1;
                        X = Xs(ii);
                        B(1:6,:) = B_HermiteSpline(X,Bdof,Bodr);

                        for ii=2:nip
                            X = Xs(ii);
                            B(6*(ii-1)+1:6*ii,:) = B_HermiteSpline(X,Bdof,Bodr);
                            X = Xs(ii-1)+Z1*(Xs(ii)-Xs(ii-1));
                            B_Z1(6*(ii-2)+1:6*(ii-1),:) = B_HermiteSpline(X,Bdof,Bodr);
                            X = Xs(ii-1)+Z2*(Xs(ii)-Xs(ii-1));
                            B_Z2(6*(ii-2)+1:6*(ii-1),:) = B_HermiteSpline(X,Bdof,Bodr);
                            X = Xs(ii-1)+Z*(Xs(ii)-Xs(ii-1));
                            B_Z(6*(ii-2)+1:6*(ii-1),:)  = B_HermiteSpline(X,Bdof,Bodr);                                
                        end
                        file = 'B_HermiteSpline';
                        Bh   = str2func(['@(X,Bdof,Bodr)',file,'(X,Bdof,Bodr)']);    
                    case 'Gaussian Element'
                        
                        n_ele = max(Bodr); % All numbers must be same max is used considering user error
                        nGauss = 5; %minimum
                        
                        dof      = sum(Bdof*(n_ele+1));
                        if nGauss/n_ele-floor(nGauss/n_ele)==0
                            nGausse_min = floor(nGauss/n_ele);
                        else
                            nGausse_min = floor(nGauss/n_ele)+1;
                        end
                        nGausse_min(nGausse_min<3)=3;
                        
                        prompt           = ['Gauss Quadrature points per element (min ', num2str(nGausse_min),'):'];
                        dlgtitle         = 'Gaussian Points ';
                        definput         = {num2str(nGausse_min)};
                        opts.Interpreter = 'tex';
                        answer2 = inputdlg(prompt,dlgtitle,[1 70],definput,opts);
                        nGausse = str2num(answer2{1});
                        
                        [Xse,Wse,npe]=GaussQuadrature(nGausse);
                        
                        nip = n_ele*(npe-1)+1;
                        Ws = zeros(nip,1);
                        Xs = zeros(nip,1);
                        Ws(1:end-1)=repmat(Wse(1:end-1),n_ele,1)/n_ele;
                        
                        X0 = 0;
                        for i=1:n_ele
                            Xs((i-1)*(npe-1)+1:i*(npe-1)) = X0+Xse(1:npe-1)/n_ele;
                            X0 = X0+1/n_ele;
                        end
                        Xs(end)=1;
                        
                        B        = zeros(nip*6,dof);
                        B_Z1     = zeros(nip*6,dof);
                        B_Z2     = zeros(nip*6,dof);
                        B_Z      = zeros(nip*6,dof);

                        ii = 1;
                        X = Xs(ii);
                        B(1:6,:) = B_Gaussian(X,Bdof,Bodr);
                        for ii=2:nip
                            X = Xs(ii);
                            B(6*(ii-1)+1:6*ii,:) = B_Gaussian(X,Bdof,Bodr);
                            X = Xs(ii-1)+Z1*(Xs(ii)-Xs(ii-1));
                            B_Z1(6*(ii-2)+1:6*(ii-1),:) = B_Gaussian(X,Bdof,Bodr);
                            X = Xs(ii-1)+Z2*(Xs(ii)-Xs(ii-1));
                            B_Z2(6*(ii-2)+1:6*(ii-1),:) = B_Gaussian(X,Bdof,Bodr);
                            X = Xs(ii-1)+Z*(Xs(ii)-Xs(ii-1));
                            B_Z(6*(ii-2)+1:6*(ii-1),:)  = B_Gaussian(X,Bdof,Bodr);                                
                        end
                        file = 'B_Gaussian';
                        Bh   = str2func(['@(X,Bdof,Bodr)',file,'(X,Bdof,Bodr)']);
                    case 'Inverse Multi-quadratic'
                        
                        n_ele = max(Bodr); % All numbers must be same max is used considering user error
                        nGauss = 5; %minimum
                        
                        dof      = sum(Bdof*(n_ele+1));
                        if nGauss/n_ele-floor(nGauss/n_ele)==0
                            nGausse_min = floor(nGauss/n_ele);
                        else
                            nGausse_min = floor(nGauss/n_ele)+1;
                        end
                        nGausse_min(nGausse_min<3)=3;
                        
                        prompt           = ['Gauss Quadrature points per element (min ', num2str(nGausse_min),'):'];
                        dlgtitle         = 'Gaussian Points ';
                        definput         = {num2str(nGausse_min)};
                        opts.Interpreter = 'tex';
                        answer2 = inputdlg(prompt,dlgtitle,[1 70],definput,opts);
                        nGausse = str2num(answer2{1});
                        
                        [Xse,Wse,npe]=GaussQuadrature(nGausse);
                        
                        nip = n_ele*(npe-1)+1;
                        Ws = zeros(nip,1);
                        Xs = zeros(nip,1);
                        Ws(1:end-1)=repmat(Wse(1:end-1),n_ele,1)/n_ele;
                        
                        X0 = 0;
                        for i=1:n_ele
                            Xs((i-1)*(npe-1)+1:i*(npe-1)) = X0+Xse(1:npe-1)/n_ele;
                            X0 = X0+1/n_ele;
                        end
                        Xs(end)=1;
                        
                        B    = zeros(nip*6,dof);
                        B_Z1 = zeros(nip*6,dof);
                        B_Z2 = zeros(nip*6,dof);
                        B_Z  = zeros(nip*6,dof);

                        ii = 1;
                        X = Xs(ii);
                        B(1:6,:) = B_IMQ(X,Bdof,Bodr);
                        for ii=2:nip
                            X = Xs(ii);
                            B(6*(ii-1)+1:6*ii,:) = B_IMQ(X,Bdof,Bodr);
                            X = Xs(ii-1)+Z1*(Xs(ii)-Xs(ii-1));
                            B_Z1(6*(ii-2)+1:6*(ii-1),:) = B_IMQ(X,Bdof,Bodr);
                            X = Xs(ii-1)+Z2*(Xs(ii)-Xs(ii-1));
                            B_Z2(6*(ii-2)+1:6*(ii-1),:) = B_IMQ(X,Bdof,Bodr);
                            X = Xs(ii-1)+Z*(Xs(ii)-Xs(ii-1));
                            B_Z(6*(ii-2)+1:6*(ii-1),:)  = B_IMQ(X,Bdof,Bodr);                                
                        end
                        file = 'B_IMQ';
                        Bh   = str2func(['@(X,Bdof,Bodr)',file,'(X,Bdof,Bodr)']);
                    case 'Custom Independent'
                        %select function, find dof
                        file = 'B_CustomIndependent';
                        Bh   = str2func(file);
                        B    = Bh(0);
                        dof  = size(B,2);
                        B    = zeros(nip*6,dof);
                        B_Z1 = zeros(nip*6,dof);
                        B_Z2 = zeros(nip*6,dof);
                        B_Z  = zeros(nip*6,dof);
                        
                        ii = 1;
                        X = Xs(ii);
                        B(1:6,:) = Bh(X);
                        for ii=2:nip
                            X = Xs(ii);
                            B(6*(ii-1)+1:6*ii,:) = Bh(X);
                            X = Xs(ii-1)+Z1*(Xs(ii)-Xs(ii-1));
                            B_Z1(6*(ii-2)+1:6*(ii-1),:) = Bh(X);
                            X = Xs(ii-1)+Z2*(Xs(ii)-Xs(ii-1));
                            B_Z2(6*(ii-2)+1:6*(ii-1),:) = Bh(X);
                            X = Xs(ii-1)+Z*(Xs(ii)-Xs(ii-1));
                            B_Z(6*(ii-2)+1:6*(ii-1),:)  = Bh(X);                                
                        end
                end
                
                syms X;
                xi_stars1 = str2sym(xi_stars{1});
                xi_stars2 = str2sym(xi_stars{2});
                xi_stars3 = str2sym(xi_stars{3});
                xi_stars4 = str2sym(xi_stars{4});
                xi_stars5 = str2sym(xi_stars{5});
                xi_stars6 = str2sym(xi_stars{6});
                xi_stars  = [xi_stars1;xi_stars2;xi_stars3;xi_stars4;xi_stars5;xi_stars6];

                if any(has(xi_stars,X))
                    xi_starfn = matlabFunction(xi_stars);
                else 
                    xi_starfn = str2func(['@(X) [' num2str(double(xi_stars')) ']''']);
                end

                T.xi_starfn = xi_starfn;

                xi_star        = zeros(6*nip,4); % precomputation at all gauss and zannah gauss points
                xi_star(1:6,1) = xi_starfn(Xs(1)); %1
                for ii=2:nip
                    xi_star((ii-1)*6+1:ii*6,1)     = xi_starfn(Xs(ii)); %2 to nGauss
                    xi_star((ii-2)*6+1:(ii-1)*6,2) = xi_starfn(Xs(ii-1)+Z1*(Xs(ii)-Xs(ii-1))); %1 to nGauss-1
                    xi_star((ii-2)*6+1:(ii-1)*6,3) = xi_starfn(Xs(ii-1)+Z2*(Xs(ii)-Xs(ii-1))); %1 to nGauss-1
                    xi_star((ii-2)*6+1:(ii-1)*6,4) = xi_starfn(Xs(ii-1)+Z*(Xs(ii)-Xs(ii-1))); %1 to nGauss-1
                end
                
                [Ms,Es,Gs] = MEG(Link,j,Xs);

                T.xi_star = xi_star;
                T.nip = nip;
                T.Xs = Xs;
                T.Ws = Ws;
                T.Ms = Ms;
                T.Es = Es;
                T.Gs = Gs;
                T.B    = B;
                T.B_Z1 = B_Z1;
                T.B_Z2 = B_Z2;
                T.B_Z  = B_Z;
                T.Bh   = Bh;
                T.dof  = dof;
                T.CI = false;    
            elseif nargin==1
                
                T.B = varargin{1};

            elseif nargin==0
                
                xi_starfn   = str2func('@(X) [0 0 0 1 0 0]''');
                T.B         = [];
                T.xi_starfn = xi_starfn;
                T.Bodr      = zeros(6,1);
                T.Bdof      = zeros(6,1);
                T.dof       = 0;
                
            end
            
        end
    end
        
    methods
        %% Set and property update functions
        function UpdateBh(T) % change integration scheme
            if isempty(T.dof)
                return
            end
            switch T.Type
                case 'Monomial'
                    file = 'B_Monomial';
                    T.Bh   = str2func(['@(X,Bdof,Bodr)',file,'(X,Bdof,Bodr)']);
                case 'Legendre Polynomial'
                    file = 'B_LegendrePolynomial';
                    T.Bh   = str2func(['@(X,Bdof,Bodr)',file,'(X,Bdof,Bodr)']);
                case 'Chebychev'
                    file = 'B_Chebychev';
                    T.Bh   = str2func(['@(X,Bdof,Bodr)',file,'(X,Bdof,Bodr)']);
                case 'Fourier'
                    file = 'B_Fourier';
                    T.Bh   = str2func(['@(X,Bdof,Bodr)',file,'(X,Bdof,Bodr)']);
                case 'FEM Like'
                    file = 'B_FEMLike';
                    T.Bh = str2func(['@(X,Bdof,Bodr,SubClass)',file,'(X,Bdof,Bodr,SubClass)']);
                case 'Hermite Spline'
                    file = 'B_Gaussian';
                    T.Bh   = str2func(['@(X,Bdof,Bodr)',file,'(X,Bdof,Bodr)']);
                case 'Gaussian Element'
                    file = 'B_Gaussian';
                    T.Bh   = str2func(['@(X,Bdof,Bodr)',file,'(X,Bdof,Bodr)']);
                case 'Inverse Multi-quadratic'
                    file = 'B_IMQ';
                    T.Bh   = str2func(['@(X,Bdof,Bodr)',file,'(X,Bdof,Bodr)']);
                case 'Custom Independent'
                    file = 'B_CustomIndependent';
                    T.Bh   = str2func(file);
            end
        end
        
%%%%%%%%%%%%%%%%%%%%%

        function UpdateIntegration(T) %updates np, Xs, Ws, and dof
            if isempty(T.dof)
                return
            end
            switch T.Type
                case {'Monomial','Legendre Polynomial','Chebychev','Fourier','FEM Like','Hermite Spline','Gaussian Element','Inverse Multi-quadratic'}
                    if any(strcmp(T.Type,{'Monomial','Legendre Polynomial','Chebychev','Fourier'})) %Fourier maybe separate
                        
                        if strcmp(T.Type,'Fourier')
                            nGauss = max(T.Bodr)*4;
                        else
                            nGauss = max(T.Bodr)+1;
                        end
                        nGauss(nGauss<5)=5;
                        prompt           = ['Number of Gaussian Quadrature points (min ',num2str(nGauss),'):'];
                        dlgtitle         = 'Gaussian Points ';
                        definput         = {num2str(nGauss)};
                        opts.Interpreter = 'tex';
                        answer2 = inputdlg(prompt,dlgtitle,[1 70],definput,opts);
                        nGauss = str2num(answer2{1});
                        [T.Xs,T.Ws,T.nip]=GaussQuadrature(nGauss);
                    else % 'FEM Like','Hermite Spline','Gaussian Element','Inverse Multi-quadratic'
                        n_ele = max(T.Bodr); % All numbers must be same max is used considering user error
                        nGauss = 5; %minimum
                        if strcmp(T.Type,'FEM Like')
                            if strcmp(T.SubClass,'Linear')
                                if nGauss/n_ele-floor(nGauss/n_ele)==0
                                    nGausse_min = floor(nGauss/n_ele);
                                else
                                    nGausse_min = floor(nGauss/n_ele)+1;
                                end
                                nGausse_min(nGausse_min<2)=2;
                            elseif strcmp(T.SubClass,'Quadratic')
                                if nGauss/n_ele-floor(nGauss/n_ele)==0
                                    nGausse_min = floor(nGauss/n_ele);
                                else
                                    nGausse_min = floor(nGauss/n_ele)+1;
                                end
                                nGausse_min(nGausse_min<3)=3;
                            else % cubic
                                if nGauss/n_ele-floor(nGauss/n_ele)==0
                                    nGausse_min = floor(nGauss/n_ele);
                                else
                                    nGausse_min = floor(nGauss/n_ele)+1;
                                end
                                nGausse_min(nGausse_min<4)=4;
                            end
                        elseif strcmp(T.Type,'Hermite Spline')
                            if nGauss/n_ele-floor(nGauss/n_ele)==0
                                nGausse_min = floor(nGauss/n_ele);
                            else
                                nGausse_min = floor(nGauss/n_ele)+1;
                            end
                            nGausse_min(nGausse_min<4)=4;
                        else% 'Gaussian Element','Inverse Multi-quadratic'
                            if nGauss/n_ele-floor(nGauss/n_ele)==0
                                nGausse_min = floor(nGauss/n_ele);
                            else
                                nGausse_min = floor(nGauss/n_ele)+1;
                            end
                            nGausse_min(nGausse_min<3)=3;
                        end
                        
                        prompt           = ['Gauss Quadrature points per element (min ', num2str(nGausse_min),'):'];
                        dlgtitle         = 'Gaussian Points ';
                        definput         = {num2str(nGausse_min)};
                        opts.Interpreter = 'tex';
                        answer2 = inputdlg(prompt,dlgtitle,[1 70],definput,opts);
                        nGausse = str2num(answer2{1});
                        
                        [Xse,Wse,npe]=GaussQuadrature(nGausse);
                        T.nip = n_ele*(npe-1)+1;
                        T.Ws = zeros(T.nip,1);
                        T.Xs = zeros(T.nip,1);
                        T.Ws(1:end-1)=repmat(Wse(1:end-1),n_ele,1)/n_ele;
                        
                        X0 = 0;
                        for i=1:n_ele
                            T.Xs((i-1)*(npe-1)+1:i*(npe-1)) = X0+Xse(1:npe-1)/n_ele;
                            X0 = X0+1/n_ele;
                        end
                        T.Xs(end)=1;
                    end
            end
        end
%%%%%%%%%%%%%%%%%%%        
        function Updatedof(T)
            if isempty(T.dof)
                return
            end
            if any(strcmp(T.Type,{'Monomial','Legendre Polynomial','Chebychev'}))
                T.dof  = sum(T.Bdof.*(T.Bodr+1));
            elseif strcmp(T.Type,'Fourier')
                T.dof  = sum(T.Bdof.*(2*T.Bodr+1));
            elseif strcmp(T.Type,'FEM Like')
                n_ele = max(T.Bodr);
                n_ele(n_ele==0)=1;
                if strcmp(T.SubClass,'Linear')
                    T.dof = sum(T.Bdof*(n_ele+1));
                elseif strcmp(T.SubClass,'Quadratic')
                    T.dof = sum(T.Bdof*(2*n_ele+1));
                else % cubic
                    T.dof = sum(T.Bdof*(3*n_ele+1));
                end
            elseif strcmp(T.Type,'Hermite Spline')
                n_ele = max(T.Bodr);
                n_ele(n_ele==0)=1;
                T.dof  = sum(T.Bdof*(2*n_ele+2));  
            elseif any(strcmp(T.Type,{'Gaussian Element','Inverse Multi-quadratic'})) 
                n_ele = max(T.Bodr);
                n_ele(n_ele==0)=1;
                T.dof  = sum(T.Bdof*(n_ele+1));    
            elseif strcmp(T.Type,'Custom Independent')
                B0     = T.Bh(0);
                T.dof  = size(B0,2); 
            end
            
        end
%%%%%%%%%%%%%%%%%%%%%%%%%        
        function UpdateBs(T)
            if isempty(T.dof)
                return
            end
            Z1     = 1/2-sqrt(3)/6;          %Zanna quadrature coefficient
            Z2     = 1/2+sqrt(3)/6;          %Zanna quadrature coefficient
            Z      = 1/2;                    %Zanna quadrature coefficient 2nd order
            if any(strcmp(T.Type,{'Monomial','Legendre Polynomial','Chebychev','Fourier','FEM Like','Hermite Spline','Gaussian Element','Inverse Multi-quadratic','Custom Independent'}))

                T.B    = zeros(T.nip*6,T.dof);
                T.B_Z1 = zeros(T.nip*6,T.dof);
                T.B_Z2 = zeros(T.nip*6,T.dof);
                T.B_Z  = zeros(T.nip*6,T.dof);

                ii = 1;
                X = T.Xs(ii);
                if strcmp(T.Type,'FEM Like')
                    T.B(1:6,:) = T.Bh(X,T.Bdof,T.Bodr,T.SubClass);
                    for ii=2:T.nip
                        X = T.Xs(ii);
                        T.B(6*(ii-1)+1:6*ii,:) = T.Bh(X,T.Bdof,T.Bodr,T.SubClass);
                        X = T.Xs(ii-1)+Z1*(T.Xs(ii)-T.Xs(ii-1));
                        T.B_Z1(6*(ii-2)+1:6*(ii-1),:) = T.Bh(X,T.Bdof,T.Bodr,T.SubClass);
                        X = T.Xs(ii-1)+Z2*(T.Xs(ii)-T.Xs(ii-1));
                        T.B_Z2(6*(ii-2)+1:6*(ii-1),:) = T.Bh(X,T.Bdof,T.Bodr,T.SubClass);
                        X = T.Xs(ii-1)+Z*(T.Xs(ii)-T.Xs(ii-1));
                        T.B_Z(6*(ii-2)+1:6*(ii-1),:)  = T.Bh(X,T.Bdof,T.Bodr,T.SubClass);                                
                    end
                elseif strcmp(T.Type,'Custom Independent')
                    T.B(1:6,:) = T.Bh(X);
                    for ii=2:T.nip
                        X = T.Xs(ii);
                        T.B(6*(ii-1)+1:6*ii,:) = T.Bh(X);
                        X = T.Xs(ii-1)+Z1*(T.Xs(ii)-T.Xs(ii-1));
                        T.B_Z1(6*(ii-2)+1:6*(ii-1),:) = T.Bh(X);
                        X = T.Xs(ii-1)+Z2*(T.Xs(ii)-T.Xs(ii-1));
                        T.B_Z2(6*(ii-2)+1:6*(ii-1),:) = T.Bh(X);
                        X = T.Xs(ii-1)+Z*(T.Xs(ii)-T.Xs(ii-1));
                        T.B_Z(6*(ii-2)+1:6*(ii-1),:)  = T.Bh(X);                                
                    end
                else
                    T.B(1:6,:) = T.Bh(X,T.Bdof,T.Bodr);
                    for ii=2:T.nip
                        X = T.Xs(ii);
                        T.B(6*(ii-1)+1:6*ii,:) = T.Bh(X,T.Bdof,T.Bodr);
                        X = T.Xs(ii-1)+Z1*(T.Xs(ii)-T.Xs(ii-1));
                        T.B_Z1(6*(ii-2)+1:6*(ii-1),:) = T.Bh(X,T.Bdof,T.Bodr);
                        X = T.Xs(ii-1)+Z2*(T.Xs(ii)-T.Xs(ii-1));
                        T.B_Z2(6*(ii-2)+1:6*(ii-1),:) = T.Bh(X,T.Bdof,T.Bodr);
                        X = T.Xs(ii-1)+Z*(T.Xs(ii)-T.Xs(ii-1));
                        T.B_Z(6*(ii-2)+1:6*(ii-1),:)  = T.Bh(X,T.Bdof,T.Bodr);                                
                    end
                end
            end  
        end
%%%%%%%%%%%%%%%%%%%        
        function UpdateMEG(T)
            if isempty(T.dof)
                return
            end
            [T.Ms,T.Es,T.Gs] = MEG(T.Link,T.div,T.Xs);
        end
%%%%%%%%%%%%%%%%%%        
        function Updatexi_star(T)

            if isempty(T.dof)
                return
            end
            Z1     = 1/2-sqrt(3)/6;          %Zanna quadrature coefficient
            Z2     = 1/2+sqrt(3)/6;          %Zanna quadrature coefficient
            Z      = 1/2;                    %Zanna quadrature coefficient 2nd order

            T.xi_star        = zeros(6*T.nip,4);
            T.xi_star(1:6,1) = T.xi_starfn(T.Xs(1)); %1
            for ii=2:T.nip
                T.xi_star((ii-1)*6+1:ii*6,1)     = T.xi_starfn(T.Xs(ii)); %2 to nGauss
                T.xi_star((ii-2)*6+1:(ii-1)*6,2) = T.xi_starfn(T.Xs(ii-1)+Z1*(T.Xs(ii)-T.Xs(ii-1))); %1 to nGauss-1
                T.xi_star((ii-2)*6+1:(ii-1)*6,3) = T.xi_starfn(T.Xs(ii-1)+Z2*(T.Xs(ii)-T.Xs(ii-1))); %1 to nGauss-1
                T.xi_star((ii-2)*6+1:(ii-1)*6,4) = T.xi_starfn(T.Xs(ii-1)+Z*(T.Xs(ii)-T.Xs(ii-1))); %1 to nGauss-1
            end
        end
        function Add_more_X(T)
            for k=1:size(T.Xadd)
                X1 = T.Xadd(k);
                if ~any(T.Xs==X1) %if already present don't add
                    iv = find(T.Xs>X1);
                    iX1 = iv(1);
                    T.Xs = [T.Xs(1:(iX1-1));X1;T.Xs(iX1:end)];
                    T.Ws = [T.Ws(1:(iX1-1));0;T.Ws(iX1:end)];
                    T.nip = T.nip+1;
                end
            end
        end
%%%%%%%%%%%%%%%%%%%        
        function UpdateAll(T)
            if isempty(T.dof)
                return
            end
            T.Updatedof();
            T.UpdateIntegration();
            T.Add_more_X();
            T.UpdateBs();
            T.Updatexi_star();
            T.UpdateMEG();
        end
        function UpdatePreCompute(T)
            if ~isempty(T.Type) %This was run when Closed loop joint twist was loaded!!! if condition avoids the change required
                T.UpdateBs();
                T.Updatexi_star();
                T.UpdateMEG();
            end
        end
        
        
        
        function set.Bdof(T, value)
            T.Bdof = value;
            T.Updatedof();
            T.UpdateBs();
        end
        function set.Bodr(T, value)
            T.Bodr = value;
            T.UpdateAll();
        end
        function set.Type(T, value)
            T.Type = value;
            T.UpdateBh();
        end
        function set.xi_starfn(T, value)
            T.xi_starfn = value;
            T.Updatexi_star();
        end
        
        function set.Bh(T, value)
            T.Bh = value;
            T.UpdateAll();
        end
        function set.SubClass(T, value)
            T.SubClass = value;
            T.UpdateAll();
        end
        
        function set.Xadd(T,value)
            T.Xadd = value;
            T.Add_more_X();
            T.UpdatePreCompute();
        end
        function s = saveobj(T) 
            % save properties into a struct to avoid set method when loading
            % Pass that struct onto the SAVE command.
            s.Type = T.Type;
            s.SubClass = T.SubClass;
            s.Bdof = T.Bdof;
            s.Bodr = T.Bodr;
            s.dof = T.dof;
            s.Bh = T.Bh;
            s.B = T.B;
            s.B_Z1 = T.B_Z1;
            s.B_Z2 = T.B_Z2;
            s.B_Z = T.B_Z;
            s.xi_starfn = T.xi_starfn;
            s.xi_star = T.xi_star;
            s.Link = T.Link;
            s.div = T.div;
            s.nip = T.nip;
            s.Xs = T.Xs;
            s.Ws = T.Ws;
            s.Ms = T.Ms;
            s.Es = T.Es;
            s.Gs = T.Gs;
            s.Xadd = T.Xadd;
            s.CI = T.CI;
            s.CIFn = T.CIFn;
        end
  
    end
    
    methods(Static)
        
        function T = loadobj(s) 
            T = SorosimTwist;
            T.dof = []; %to overload set methods when loading
            T.Type = s.Type;
            T.SubClass = s.SubClass;
            T.Bdof = s.Bdof;
            T.Bodr = s.Bodr;
            T.Bh = s.Bh;
            T.B = s.B;
            T.B_Z1 = s.B_Z1;
            T.B_Z2 = s.B_Z2;
            T.B_Z = s.B_Z;
            T.xi_starfn = s.xi_starfn;
            T.xi_star = s.xi_star;
            T.Link = s.Link;
            T.div = s.div;
            T.nip = s.nip;
            T.Xs = s.Xs;
            T.Ws = s.Ws;
            T.Ms = s.Ms;
            T.Es = s.Es;
            T.Gs = s.Gs;
            T.Xadd = s.Xadd;
            T.CI = s.CI;
            T.CIFn = s.CIFn;
            T.dof = s.dof;
        end
    
    end
end
