%Class to define a SorosimLink
%A SorosimLink consists of a rigid joint and a body.
%Last modified by Anup Teejo Mathew 02.03.2022
classdef SorosimLink < handle
    
    properties
        
        %General Properties
        
        jointtype    %Type of joint used to connect the link (lumped DoF). (R) for Revolute,(P) for Prismatic, (H) for Helical, (U) for Universal, (C) for Cylindrical, (A) for Planar, (S) for Spherical, (F) for Free motion and (N) for Fixed
        linktype     %'s' for soft or 'r' for rigid
        CS           %Cross section: 'R' for rectangular 'C' for circular 'E' for elliptical
        npie         %Number of pieces. For a rigid link, npie = 1. For a soft link, npie=1+number of divisions 
        nGauss       %Number of Gaussian Quadrature points + 2
        Xs           %Gaussian Quadrature points inc 0 and 1
        Ws           %Gaussian Quadrature weights
        
        
        %Geometric Properties
        
        lp        %Length of each divisions of the link (soft link) [m]
        L         %Total length of the link [m]
        r         %Radius as a function of X1 (X1=X/L, X1 varies from 0 to 1) [m]
        h         %Height as a function of X1 (X1=X/L, X1 varies from 0 to 1) [m]
        w         %Width as a function of X1 (X1=X/L, X1 varies from 0 to 1) [m]
        a         %Semi-major axis as a function of X1 (X1=X/L, X1 varies from 0 to 1) [m]
        b         %Semi-minor axis as a function of X1 (X1=X/L, X1 varies from 0 to 1) [m]
        gi        %Transformation from joint to center of mass for ridig link to center of area for soft link
        gf        %Transformation to joint from center of mass for ridig link from center of area for soft link
        
        %Material
        
        E         %Young's modulus [Pa]
        Poi       %Poisson's ratio [-]
        G         %Shear modulus [Pa]
        Eta       %Material Damping [Pa.s]
        Rho       %Density [kg/m^3]
        Kj        %Joint Stiffness Matrix
        
        Ms        %Inertia matrix (cross sectional for soft link)
        Es        %Cross sectional stiffness matrix
        Gs        %Cross sectional damping matrix
        
        %Plot Properties
        
        n_l       %Number of cross sections per division. (default value: 10)
        n_r       %Number of radial points if the cross section is circular or ellipsoidal (default value: 18)
        color     %color of link (random by default)
        
    end
    
    
    %%
    methods
        function Li = SorosimLink
            
            quest = 'Select link type:';
            answ  = questdlg(quest,'Link type',...
                'Soft','Rigid','Soft');
            
            switch answ
                case 'Rigid'
                    
                    Li.linktype='r';
                    
                    badanswer = true;
                    
                    while badanswer
                        
                        prompt           = {'Joint type (Enter (R) for Revolute,(P) for Prismatic, (H) for Helical, (U) for Universal, (C) for Cylindrical, (A) for Planar, (S) for Spherical, (F) for Free motion and (N) for Fixed): ',...
                                            'Joint stiffness (matrix) [Nm/rad or N/m]:',...
                                            'Density [kg/m^3]:','Cross-section shape (C for circular, R for rectangular, E for elliptical):','y coordinate of origin [m] wrt. previous frame:','z coordinate of origin [m] wrt. previous frame:','Length of rigid link [m]:'};
                        dlgtitle         = 'Rigid Link Properties';
                        definput         = {'N','0','1000','C','0','0','0.3'};
                        opts.Interpreter = 'tex';
                        answer           = inputdlg(prompt,dlgtitle,[1 60],definput,opts);
                        
                        jointtype    = answer{1};
                        Kj           = str2num(answer{2});
                        Rho          = str2double(answer{3}); 
                        CS           = answer{4};
                        cy           = str2double(answer{5});
                        cz           = str2double(answer{6});
                        L            = str2double(answer{7});
                        
                        if ~any(strcmp(jointtype,{'R','P','H','U','C','A','S','F','N'}))
                            uiwait(msgbox('Choose joint type from the given options','Error','error'));
                        elseif ~any(strcmp(CS,{'C','R','E'}))
                            uiwait(msgbox('Choose cross-section shape from the given options','Error','error'));
                        elseif any(Kj,'all')
                            switch jointtype
                                case 'N'
                                    uiwait(msgbox('Cannot assign a stiffness for a fixed joint','Error','error'));
                                case 'R'
                                    if ~isequal(size(Kj),[1,1])
                                        uiwait(msgbox('Incorrect stiffness matrix dimension for a revolute joint','Error','error'));
                                    else
                                        badanswer = false;
                                    end
                                case 'P'
                                    if ~isequal(size(Kj),[1,1])
                                        uiwait(msgbox('Incorrect stiffness matrix dimension for a prismatic joint','Error','error'));
                                    else
                                        badanswer = false;
                                    end
                                case 'H'
                                    if ~isequal(size(Kj),[1,1])
                                        uiwait(msgbox('Incorrect stiffness matrix dimension for a helical joint','Error','error'));
                                    else
                                        badanswer = false;
                                    end
                                case 'U'
                                    if ~isequal(size(Kj),[2,2])
                                        uiwait(msgbox('Incorrect stiffness matrix dimension for a universal joint','Error','error'));
                                    else
                                        badanswer = false;
                                    end
                                case 'C'
                                    if ~isequal(size(Kj),[2,2])
                                        uiwait(msgbox('Incorrect stiffness matrix dimension for a cylindrical joint','Error','error'));
                                    else
                                        badanswer = false;
                                    end
                                case 'A'
                                    if ~isequal(size(Kj),[3,3])
                                        uiwait(msgbox('Incorrect stiffness matrix dimension for a planar joint','Error','error'));
                                    else
                                        badanswer = false;
                                    end
                                case 'S'
                                    if ~isequal(size(Kj),[3,3])
                                        uiwait(msgbox('Incorrect stiffness matrix dimension for a spherical joint','Error','error'));
                                    else
                                        badanswer = false;
                                    end
                                case 'F'
                                    if ~isequal(size(Kj),[6,6])
                                        uiwait(msgbox('Incorrect stiffness matrix dimension for a free joint','Error','error'));
                                    else
                                        badanswer = false;
                                    end
                            end
                        else
                            badanswer = false;
                        end
                        
                        if ~badanswer
                            
                            Li.jointtype = jointtype;
                            if jointtype=='N'
                                Li.Kj = [];
                            else
                                Li.Kj = Kj;
                            end
                            Li.npie      = 1;
                            Li.Rho       = Rho;
                            Li.CS        = CS;
                            Li.L         = L;
                        end
                    
                    end
                    
                    if CS=='C'
                        
                        prompt           = {'Initial radius [m]:','Final radius [m]:'};
                        dlgtitle         = 'Circular cross section properties';
                        definput         = {'0.03','0.03'};
                        opts.Interpreter = 'tex';
                        answer           = inputdlg(prompt,dlgtitle,[1 70],definput,opts);
                        
                        ri        = str2double(answer{1});
                        rf        = str2double(answer{2});
                        
                        syms X1
                        r   = ri+X1*(rf-ri);
                        
                        if has(r,X1)
                            r = matlabFunction(r);
                        else
                            r = str2func(['@(X1)' num2str(ri)]);
                        end
                        
                        Li.r = r;
                        
                        [Xs,Ws,nGauss] = GaussQuadrature(10);
                        V      = 0;
                        M_area = 0;
                        
                        for ii=2:nGauss-1
                            r_here = r(Xs(ii));
                            A_here = pi*r_here^2;
                            V      = V + Ws(ii)*L*A_here;
                            M_area = M_area+ Ws(ii)*L*A_here*Xs(ii)*L;
                        end
                        
                        cx   = M_area/V;
                        mass = V*Rho;
                        
                        Ix = 0;
                        Iy = 0;
                        for ii=2:nGauss-1
                            r_here = r(Xs(ii));
                            Iy     = Iy+Ws(ii)*L*((pi/4)*r_here^4+pi*r_here^2*(Xs(ii)*L-cx)^2);
                            Ix     = Ix+Ws(ii)*L*((pi/2)*r_here^4);
                        end
                        Iz = Iy;
   
                    elseif CS=='R'
                        
                        prompt           = {'Initial height [m]:','Initial width [m]:','Final height [m]:','Final width [m]:'};
                        dlgtitle         = 'Rectangular cross section properties';
                        definput         = {'0.02','0.02','0.02','0.02'};
                        opts.Interpreter = 'tex';
                        answer           = inputdlg(prompt,dlgtitle,[1 70],definput,opts);
                        
                        hi       = str2double(answer{1});
                        wi       = str2double(answer{2});
                        hf       = str2double(answer{3});
                        wf       = str2double(answer{4});
                        
                        syms X1
                        h   = hi+X1*(hf-hi);
                        w   = wi+X1*(wf-wi);
                        
                        if has(h,X1)
                            h = matlabFunction(h);
                        else
                            h = str2func(['@(X1)' num2str(hi)]);
                        end
                        
                        if has(w,X1)
                            w = matlabFunction(w);
                        else
                            w = str2func(['@(X1)' num2str(wi)]);
                        end
                        
                        Li.h = h;
                        Li.w = w;
                        
                        [Xs,Ws,nGauss] = GaussQuadrature(10);
                        V      = 0;
                        M_area = 0;
                        
                        for ii=2:nGauss-1
                            h_here = h(Xs(ii));
                            w_here = w(Xs(ii));
                            A_here = h_here*w_here;
                            V      = V + Ws(ii)*L*A_here;
                            M_area = M_area+ Ws(ii)*L*A_here*Xs(ii)*L;
                        end
                        
                        cx   = M_area/V;
                        mass = V*Rho;
                        
                        Ix = 0;
                        Iy = 0;
                        Iz = 0;
                        for ii=2:nGauss-1
                            h_here = h(Xs(ii));
                            w_here = w(Xs(ii));
                            Iy     = Iy+Ws(ii)*L*((Xs(ii)*L-cx)^2*h_here*w_here+h_here*w_here^3/12);
                            Iz     = Iz+Ws(ii)*L*((Xs(ii)*L-cx)^2*h_here*w_here+h_here^3*w_here/12);
                            Ix     = Iz+Ws(ii)*L*(h_here*w_here^3/12+h_here^3*w_here/12);
                        end
                        
                    elseif CS=='E'
                        
                        prompt           = {'Initial semi-major axis [m]:','Initial semi-minor axis [m]:','Final semi-major axis [m]:','Final semi-minor axis [m]:'};
                        dlgtitle         = 'Elliptical cross section properties';
                        definput         = {'0.04','0.02','0.04','0.02'};
                        opts.Interpreter = 'tex';
                        answer           = inputdlg(prompt,dlgtitle,[1 70],definput,opts);
                        
                        ai       = str2double(answer{1});
                        bi       = str2double(answer{2});
                        af       = str2double(answer{3});
                        bf       = str2double(answer{4});
                        
                        syms X1
                        a   = ai+X1*(af-ai);
                        b   = bi+X1*(bf-bi);
                        
                        if has(a,X1)
                            a = matlabFunction(a);
                        else
                            a = str2func(['@(X1)' num2str(ai)]);
                        end
                        
                        if has(b,X1)
                            b = matlabFunction(b);
                        else
                            b = str2func(['@(X1)' num2str(bi)]);
                        end
                        
                        Li.a = a;
                        Li.b = b;
                        
                        [Xs,Ws,nGauss] = GaussQuadrature(10);
                        V      = 0;
                        M_area = 0;
                        
                        for ii=2:nGauss-1
                            a_here = a(Xs(ii));
                            b_here = b(Xs(ii));
                            A_here = pi*a_here*b_here;
                            V      = V + Ws(ii)*L*A_here;
                            M_area = M_area+ Ws(ii)*L*A_here*Xs(ii)*L;
                        end
                        
                        cx   = M_area/V;
                        mass = V*Rho;
                        
                        Ix = 0;
                        Iy = 0;
                        Iz = 0;
                        for ii=2:nGauss-1
                            a_here = a(Xs(ii));
                            b_here = b(Xs(ii));
                            Iy     = Iy+Ws(ii)*L*((Xs(ii)*L-cx)^2*pi*a_here*b_here+pi*a_here*b_here^3/4);
                            Iz     = Iz+Ws(ii)*L*((Xs(ii)*L-cx)^2*pi*a_here*b_here+pi*a_here^3*b_here/4);
                            Ix     = Ix+Ws(ii)*L*(pi*a_here*b_here^3/4+pi*a_here^3*b_here/4);
                        end
                    end
                    
                    
                    
                    Ms   = double([Ix 0 0 0 0 0;0 Iy 0 0 0 0;0 0 Iz 0 0 0;...
                                   0 0 0 mass 0 0;0 0 0 0 mass 0;0 0 0 0 0 mass]);
                    gi = [eye(3),[cx;cy;cz];[0 0 0 1]];
                    gf = [eye(3),[L-cx;-cy;-cz];[0 0 0 1]];
                    
                    Li.gi = gi;
                    Li.gf = gf;
                    Li.Ms = Ms;
                    
                    %%
                case 'Soft'
                    
                    Li.linktype = 's';
                    badanswer = true;
                    
                    while badanswer
                        
                        prompt           = {'Joint type (Enter (R) for Revolute,(P) for Prismatic, (H) for Helical, (U) for Universal, (C) for Cylindrical, (A) for Planar, (S) for Spherical, (F) for Free motion and (N) for Fixed):',...
                                            'Joint stiffness (matrix) [Nm/rad or N/m]:',...    
                                            'Number of divisions:','Cross-section shape (C for circular, R for rectangular, E for elliptical):',...
                                            'Density (kg/m^3):','Youngs Modulus (N/m^2)','Poissons Ratio:','Material Damping (Pa.s)(for dynamic simulation):'};
                        dlgtitle         = 'Soft Link Properties';

                        definput         = {'N','0','1','C','1000','1e6','0.5','0.112e5'};
                        opts.Interpreter = 'tex';
                        answer           = inputdlg(prompt,dlgtitle,[1 60],definput,opts);

                        jointtype     = answer{1};
                        Kj            = str2num(answer{2});
                        ndiv          = str2double(answer{3});
                        CS            = answer{4};
                        Rho           = str2double(answer{5});
                        E             = str2double(answer{6});
                        Poi           = str2double(answer{7});
                        Eta           = str2double(answer{8});
                        
                        if ~any(strcmp(jointtype,{'R','P','H','U','C','A','S','F','N'}))
                            uiwait(msgbox('Choose joint type from the given options','Error','error'));
                        elseif ~any(strcmp(CS,{'C','R','E'}))
                            uiwait(msgbox('Choose cross-section shape from the given options','Error','error'));
                        elseif any(Kj,'all')
                            switch jointtype
                                case 'N'
                                    uiwait(msgbox('Cannot assign a stiffness for a fixed joint','Error','error'));
                                case 'R'
                                    if ~isequal(size(Kj),[1,1])
                                        uiwait(msgbox('Incorrect stiffness matrix dimension for a revolute joint','Error','error'));
                                    else
                                        badanswer = false;
                                    end
                                case 'P'
                                    if ~isequal(size(Kj),[1,1])
                                        uiwait(msgbox('Incorrect stiffness matrix dimension for a prismatic joint','Error','error'));
                                    else
                                        badanswer = false;
                                    end
                                case 'H'
                                    if ~isequal(size(Kj),[1,1])
                                        uiwait(msgbox('Incorrect stiffness matrix dimension for a helical joint','Error','error'));
                                    else
                                        badanswer = false;
                                    end
                                case 'U'
                                    if ~isequal(size(Kj),[2,2])
                                        uiwait(msgbox('Incorrect stiffness matrix dimension for a universal joint','Error','error'));
                                    else
                                        badanswer = false;
                                    end
                                case 'C'
                                    if ~isequal(size(Kj),[2,2])
                                        uiwait(msgbox('Incorrect stiffness matrix dimension for a cylindrical joint','Error','error'));
                                    else
                                        badanswer = false;
                                    end
                                case 'A'
                                    if ~isequal(size(Kj),[3,3])
                                        uiwait(msgbox('Incorrect stiffness matrix dimension for a planar joint','Error','error'));
                                    else
                                        badanswer = false;
                                    end
                                case 'S'
                                    if ~isequal(size(Kj),[3,3])
                                        uiwait(msgbox('Incorrect stiffness matrix dimension for a spherical joint','Error','error'));
                                    else
                                        badanswer = false;
                                    end
                                case 'F'
                                    if ~isequal(size(Kj),[6,6])
                                        uiwait(msgbox('Incorrect stiffness matrix dimension for a free joint','Error','error'));
                                    else
                                        badanswer = false;
                                    end
                            end
                        else
                            badanswer = false;
                        end
                        
                        if ~badanswer
                            
                            Li.jointtype = jointtype;
                            if jointtype=='N'
                                Li.Kj = [];
                            else
                                Li.Kj = Kj;
                            end
                            Li.npie      = ndiv+1;
                            Li.CS        = CS;
                            Li.Rho       = Rho;
                            Li.E         = E;
                            Li.Poi       = Poi;
                            Li.Eta       = Eta;
                            G            = E/(2*(1+Poi));
                            Li.G         = G;
                        end
                        
                    end
                    
                    lp   = cell(1,ndiv);
                    r    = cell(1,ndiv);
                    h    = cell(1,ndiv); w  = cell(1,ndiv);
                    a    = cell(1,ndiv); b  = cell(1,ndiv);
                    Ms   = cell(1,ndiv); Es = cell(1,ndiv); Gs = cell(1,ndiv);
                    gi   = cell(1,ndiv); gf = cell(1,ndiv);
                    
                    Xs = cell(1,ndiv); Ws = cell(1,ndiv); nGauss = cell(1,ndiv);
                    
                    for i=1:ndiv
                        
                        if CS=='C'
                            
                            prompt           = {'Length (m):','Number of Gaussian Quadrature points (min 5):','y coordinate of origin wrt. previous frame (m):','z coordinate of origin wrt. previous frame (m):'...
                                ,'Initial radius (m):','Final radius (m):'};
                            dlgtitle         = ['Enter geometric properties of division ',num2str(i)];
                            definput         = {'0.3','5','0','0','0.03','0.03'};
                            opts.Interpreter = 'tex';
                            answer2          = inputdlg(prompt,dlgtitle,[1 70],definput,opts);
                            
                            l_p   = str2double(answer2{1});
                            N_GQ  = str2num(answer2{2});
                            ri_p  = str2double(answer2{5});
                            rf_p  = str2double(answer2{6});
                                                     
                            syms X1
                            r_sym = ri_p+X1*(rf_p-ri_p);

                            if has(r_sym,X1)
                                r_fn = matlabFunction(r_sym);
                            else
                                r_fn = str2func(['@(X1)' num2str(ri_p)]);
                            end

                            r{i} = r_fn;
                            
                            [Xs_p,Ws_p,nGauss_p] = GaussQuadrature(N_GQ);
                            
                            r_nGauss = zeros(nGauss_p,1);
                            for ii=1:nGauss_p
                                r_nGauss(ii) = r_fn(Xs_p(ii));
                            end
                            Iy = (pi/4)*r_nGauss.^4;
                            Iz = Iy;
                            Ix = Iy+Iz;
                            A  = pi*r_nGauss.^2;
                            
                        elseif CS=='R'
                            
                            prompt           = {'Length (m):','Number of Gaussian Quadrature points (min 5):','y coordinate of origin wrt. previous frame (m):','z coordinate of origin wrt. previous frame (m):'...
                                ,'Initial height (m):','Initial width (m):','Final height (m):','Final width (m):'};
                            dlgtitle         = ['Enter geometric properties of division ',num2str(i)];
                            definput         = {'0.3','5','0','0','0.02','0.02','0.02','0.02'};
                            opts.Interpreter = 'tex';
                            answer2 = inputdlg(prompt,dlgtitle,[1 70],definput,opts);
                            
                            l_p    = str2double(answer2{1});
                            N_GQ   = str2num(answer2{2});
                            hi_p   = str2double(answer2{5});
                            wi_p   = str2double(answer2{6});
                            hf_p   = str2double(answer2{7});
                            wf_p   = str2double(answer2{8});
                            
                            syms X1
                            h_sym = hi_p+X1*(hf_p-hi_p);
                            w_sym = wi_p+X1*(wf_p-wi_p);

                            if has(h_sym,X1)
                                h_fn = matlabFunction(h_sym);
                            else
                                h_fn = str2func(['@(X1)' num2str(hi_p)]);
                            end

                            if has(w_sym,X1)
                                w_fn = matlabFunction(w_sym);
                            else
                                w_fn = str2func(['@(X1)' num2str(wi_p)]);
                            end
                            
                            h{i} = h_fn;
                            w{i} = w_fn;
                            
                            [Xs_p,Ws_p,nGauss_p] = GaussQuadrature(N_GQ);
                            
                            h_nGauss = zeros(nGauss_p,1);
                            w_nGauss = zeros(nGauss_p,1);
                            for ii=1:nGauss_p
                                h_nGauss(ii) = h_fn(Xs_p(ii));
                                w_nGauss(ii) = w_fn(Xs_p(ii));
                            end
                            
                            Iy     = (1/12)*h_nGauss.*(w_nGauss.^3);
                            Iz     = (1/12)*(h_nGauss.^3).*w_nGauss;
                            Ix     = Iy+Iz;
                            A      = h_nGauss.*w_nGauss;
                            
                         elseif CS=='E'
                            
                            prompt           = {'Length (m):','Number of Gaussian Quadrature points (min 5):','y coordinate of origin wrt. previous frame (m):','z coordinate of origin wrt. previous frame (m):'...
                                ,'Initial semi-major axis (m):','Initial semi-minor axis (m):','Final semi-major axis (m):','Final semi-minor axis (m):'};
                            dlgtitle         = ['Enter geometric properties of division ',num2str(i)];
                            definput         = {'0.3','5','0','0','0.04','0.02','0.04','0.02'};
                            opts.Interpreter = 'tex';
                            answer2 = inputdlg(prompt,dlgtitle,[1 70],definput,opts);
                            
                            l_p    = str2double(answer2{1});
                            N_GQ   = str2num(answer2{2});
                            ai_p   = str2double(answer2{5});
                            bi_p   = str2double(answer2{6});
                            af_p   = str2double(answer2{7});
                            bf_p   = str2double(answer2{8});
                            
                            syms X1
                            a_sym = ai_p+X1*(af_p-ai_p);
                            b_sym = bi_p+X1*(bf_p-bi_p);

                            if has(a_sym,X1)
                                a_fn = matlabFunction(a_sym);
                            else
                                a_fn = str2func(['@(X1)' num2str(ai_p)]);
                            end

                            if has(b_sym,X1)
                                b_fn = matlabFunction(b_sym);
                            else
                                b_fn = str2func(['@(X1)' num2str(bi_p)]);
                            end
                            
                            a{i} = a_fn;
                            b{i} = b_fn;
                            
                            [Xs_p,Ws_p,nGauss_p] = GaussQuadrature(N_GQ);
                            
                            a_nGauss = zeros(nGauss_p,1);
                            b_nGauss = zeros(nGauss_p,1);
                            for ii=1:nGauss_p
                                a_nGauss(ii) = a_fn(Xs_p(ii));
                                b_nGauss(ii) = b_fn(Xs_p(ii));
                            end
                            
                            Iy     = (pi/4)*a_nGauss.*(b_nGauss.^3);
                            Iz     = (pi/4)*(a_nGauss.^3).*b_nGauss;
                            Ix     = Iy+Iz;
                            A      = pi*a_nGauss.*b_nGauss;
                            
                        end
                        
                        Ms_p = zeros(6*nGauss_p,6);
                        Es_p = zeros(6*nGauss_p,6);
                        Gs_p = zeros(6*nGauss_p,6);
                        
                        for ii=1:nGauss_p
                            Ms_p((ii-1)*6+1:ii*6,:) = Rho*diag([Ix(ii),Iy(ii),Iz(ii),A(ii),A(ii),A(ii)]);
                            Es_p((ii-1)*6+1:ii*6,:) = diag([G*Ix(ii),E*Iy(ii),E*Iz(ii),E*A(ii),G*A(ii),G*A(ii)]);
                            Gs_p((ii-1)*6+1:ii*6,:) = Eta*diag([Ix(ii),3*Iy(ii),3*Iz(ii),3*A(ii),A(ii),A(ii)]);
                        end
                        
                        Ms{i}     = Ms_p;
                        Es{i}     = Es_p;
                        Gs{i}     = Gs_p;
                        lp{i}     = l_p;
                        cy        = str2double(answer2{3});
                        cz        = str2double(answer2{4});
                        gi{i}     = [eye(3),[0;cy;cz];[0 0 0 1]];
                        gf{i}     = [eye(3),[0;-cy;-cz];[0 0 0 1]];
                        Xs{i}     = Xs_p;
                        Ws{i}     = Ws_p;
                        nGauss{i} = nGauss_p;
                        
                    end
                    
                    Li.L      = sum(cell2mat(lp));
                    Li.lp     = lp;
                    Li.r      = r;
                    Li.w      = w;
                    Li.h      = h;
                    Li.a      = a;
                    Li.b      = b;
                    Li.gi     = gi;
                    Li.gf     = gf;
                    Li.Ms     = Ms;
                    Li.Es     = Es;
                    Li.Gs     = Gs;
                    Li.Xs     = Xs;
                    Li.Ws     = Ws;
                    Li.nGauss = nGauss;

            end

                Li.n_l  = 10;
                if Li.CS~='R'
                    Li.n_r = 18;
                end
                Li.color = rand(1,3);

        end
    end
    
    methods
        %% Set and property update functions
        
        %Changing Poi:
        function set.Poi(SorosimLink, value)
            SorosimLink.Poi = value;
            SorosimLink.UpdateG();
            SorosimLink.Update()
        end
        
        %Changing E:
        function set.E(SorosimLink, value)
            SorosimLink.E = value;
            SorosimLink.UpdateG();
            SorosimLink.Update()
        end
        
        %Changing Eta:
        function set.Eta(SorosimLink, value)
            SorosimLink.Eta = value;
            SorosimLink.Update()
        end
        
        %Changing G:
        function set.G(SorosimLink, value)
            SorosimLink.G = value;
            SorosimLink.Update()
        end
        
        %Changing r:
        function set.r(SorosimLink, value)
            SorosimLink.r = value;
            SorosimLink.Update()
        end
        
        %Changing w:
        function set.w(SorosimLink, value)
            SorosimLink.w = value;
            SorosimLink.Update()
        end

        %Changing h:
        function set.h(SorosimLink, value)
            SorosimLink.h = value;
            SorosimLink.Update()
        end
        
        %Changing a:
        function set.a(SorosimLink, value)
            SorosimLink.a = value;
            SorosimLink.Update()
        end

        %Changing b:
        function set.b(SorosimLink, value)
            SorosimLink.b = value;
            SorosimLink.Update()
        end
        
        %Changing L: (rigid or constant CS soft)
        function set.L(SorosimLink, value)
            SorosimLink.L = value;
            SorosimLink.Update()
        end
        
        %Changing Rho:
        function set.Rho(SorosimLink, value)
            SorosimLink.Rho = value;
            SorosimLink.Update()
        end

        %Changing lp:
        function set.lp(SorosimLink, value)
            SorosimLink.lp = value;
            SorosimLink.Update()
            SorosimLink.UpdateL()
        end
        
        %Changing nGauss:
        function set.nGauss(SorosimLink, value)
            SorosimLink.nGauss = value;
            SorosimLink.UpdateNGauss();
        end
        
        function set.jointtype(SorosimLink, value)
            SorosimLink.jointtype = value;
            SorosimLink.UpdateKj();
        end
        
        %Updating the properties
        function Update(SorosimLink)
            if isempty(SorosimLink.color)
                return
            end
            if SorosimLink.linktype=='s'
                [SorosimLink.Ms,SorosimLink.Es,SorosimLink.Gs] = LinkPropUpdate(SorosimLink);
            else
                SorosimLink.Ms = LinkPropUpdate(SorosimLink);
            end
        end
        
        %Update nGauss
        function UpdateNGauss(SorosimLink)
            if isempty(SorosimLink.color)
                return
            end
            for i=1:SorosimLink.npie-1
                [Xs_p,Ws_p] = GaussQuadrature(SorosimLink.nGauss{i}-2);
                SorosimLink.Xs{i}  = Xs_p;
                SorosimLink.Ws{i}  = Ws_p;
            end
            [SorosimLink.Ms,SorosimLink.Es,SorosimLink.Gs] = LinkPropUpdate(SorosimLink);
        end
        
        %Update E or Poi
        function UpdateG(SorosimLink)
            if isempty(SorosimLink.color)
                return
            end
            SorosimLink.G=SorosimLink.E/(2*(1+SorosimLink.Poi));
        end
        
        %Update lp
        function UpdateL(SorosimLink)
            SorosimLink.L = sum(cell2mat(SorosimLink.lp));
        end
        function UpdateKj(SorosimLink)
            SorosimLink.Kj = 0;
        end
        
        
    end
    
end