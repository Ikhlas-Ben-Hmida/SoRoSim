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
        Dj        %Joint Damping Matrix
        
        M        %Inertia matrix (only for rigid body)
        
        %Plot Properties
        
        n_l       %Number of cross sections per division. (default value: 10)
        n_r       %Number of radial points if the cross section is circular or ellipsoidal (default value: 18)
        color     %Color of link (random by default)
        
        CPF       %Custom plot function for rigid bodies (logical 1 or 0)
        PlotFn    %Handle of function to plot the geometry (for rigid link)
        Lscale    %Scaling factor for plotting symbols or axes
        
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
                        
                        prompt           = {'Joint type (Enter (R) for Revolute,(P) for Prismatic, (H) for Helical, (C) for Cylindrical, (A) for Planar, (S) for Spherical, (F) for Free motion and (N) for Fixed): ',...
                                            'Cross-section shape (C for circular, R for rectangular, E for elliptical):','Joint stiffness (dof x dof matrix) [Nm/rad or N/m]:',...
                                            'Density [kg/m^3]:','Length of rigid link [m]:'};
                        dlgtitle         = 'Rigid Link Properties';
                        definput         = {'N','C','0','1000','0.3'};
                        opts.Interpreter = 'tex';
                        answer           = inputdlg(prompt,dlgtitle,[1 60],definput,opts);
                        
                        jointtype    = answer{1};
                        CS           = answer{2};
                        Kj           = str2num(answer{3});
                        Rho          = str2double(answer{4}); 
                        L            = str2double(answer{5});
                        
                        badanswer = LinkInputCheck(jointtype,CS,Kj,badanswer);
                        
                        if L<0
                            uiwait(msgbox('In a classical universe, length cannot be negative','Error','error'));
                            badanswer = true;
                        end
                        if ~badanswer
                            
                            Li.jointtype = jointtype;
                            if ~(jointtype=='N')
                                Li.Kj = Kj;
                                Li.Dj = 0;%default
                            end
                            Li.npie      = 1;
                            Li.Rho       = Rho;
                            Li.CS        = CS;
                            Li.L         = L;
                        end
                    
                    end
                    if L>0
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

                            if V~=0
                                cx   = M_area/V;
                            else
                                if L==0
                                    cx = 0;
                                else
                                    cx = L/2;
                                end
                            end

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

                            if V~=0
                                cx   = M_area/V;
                            else
                                if L==0
                                    cx = 0;
                                else
                                    Vh      = 0;
                                    Mh_area = 0;
                                    Vw      = 0;
                                    Mw_area = 0;
                                    for ii=2:nGauss-1
                                        h_here = h(Xs(ii));
                                        w_here = w(Xs(ii));
                                        Vh      = Vh + Ws(ii)*L*h_here;
                                        Mh_area = Mh_area+ Ws(ii)*L*h_here*Xs(ii)*L;
                                        Vw      = Vw + Ws(ii)*L*w_here;
                                        Mw_area = Mw_area+ Ws(ii)*L*w_here*Xs(ii)*L;
                                    end
                                    if Vh==0
                                        cx   = Mw_area/Vw;
                                    else
                                        cx   = Mh_area/Vh;
                                    end
                                end
                            end

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

                            if V~=0
                                cx   = M_area/V;
                            else
                                if L==0
                                    cx = 0;
                                else
                                    Va      = 0;
                                    Ma_area = 0;
                                    Vb      = 0;
                                    Mb_area = 0;
                                    for ii=2:nGauss-1
                                        a_here = a(Xs(ii));
                                        b_here = b(Xs(ii));
                                        Va      = Va + Ws(ii)*L*a_here;
                                        Ma_area = Ma_area+ Ws(ii)*L*a_here*Xs(ii)*L;
                                        Vb      = Vb + Ws(ii)*L*b_here;
                                        Mb_area = Mb_area+ Ws(ii)*L*b_here*Xs(ii)*L;
                                    end
                                    if Va==0
                                        cx   = Mb_area/Vb;
                                    else
                                        cx   = Ma_area/Va;
                                    end
                                end
                            end

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
                        Ix = Ix*Rho;
                        Iy = Iy*Rho;
                        Iz = Iz*Rho;

                        M   = double([Ix 0 0 0 0 0;0 Iy 0 0 0 0;0 0 Iz 0 0 0;...
                                       0 0 0 mass 0 0;0 0 0 0 mass 0;0 0 0 0 0 mass]);
                        gi = [eye(3),[cx;0;0];[0 0 0 1]];
                        gf = [eye(3),[L-cx;0;0];[0 0 0 1]];
                        if V~=0
                            Lscale = V^(1/3);
                        else
                            Lscale = 0.1;
                        end
                    else
                        M = zeros(6,6);
                        gi = eye(4);
                        gf = eye(4);
                        Lscale = 0.1;
                    end
                    Li.gi = gi;
                    Li.gf = gf;
                    Li.M = M;
                    
                    %%
                case 'Soft'
                    
                    Li.linktype = 's';
                    badanswer = true;
                    
                    while badanswer
                        
                        prompt           = {'Joint type (Enter (R) for Revolute,(P) for Prismatic, (H) for Helical, (C) for Cylindrical, (A) for Planar, (S) for Spherical, (F) for Free motion and (N) for Fixed):',...
                                            'Cross-section shape (C for circular, R for rectangular, E for elliptical):','Joint stiffness (matrix) [Nm/rad or N/m]:',...    
                                            'Number of divisions:',...
                                            'Density (kg/m^3):','Youngs Modulus (N/m^2)','Poissons Ratio:','Material Damping (Pa.s)(for dynamic simulation):'};
                        dlgtitle         = 'Soft Link Properties';

                        definput         = {'N','C','0','1','1000','1e6','0.5','1e4'};
                        opts.Interpreter = 'tex';
                        answer           = inputdlg(prompt,dlgtitle,[1 60],definput,opts);

                        jointtype     = answer{1};
                        CS            = answer{2};
                        Kj            = str2num(answer{3});
                        ndiv          = str2double(answer{4});
                        Rho           = str2double(answer{5});
                        E             = str2double(answer{6});
                        Poi           = str2double(answer{7});
                        Eta           = str2double(answer{8});
                        
                        badanswer = LinkInputCheck(jointtype,CS,Kj,badanswer);
                        
                        if ~badanswer
                            
                            Li.jointtype = jointtype;
                            if jointtype=='N'
                                Li.Kj = [];
                            else
                                Li.Kj = Kj;
                                Li.Dj = 0;
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
                    gi   = cell(1,ndiv); gf = cell(1,ndiv);
                    
                    for i=1:ndiv
                        
                        if CS=='C'
                            
                            prompt           = {'Length (m):','Initial radius (m):','Final radius (m):'};
                            dlgtitle         = ['Enter geometric properties of division ',num2str(i)];
                            definput         = {'0.3','0.03','0.03'};
                            opts.Interpreter = 'tex';
                            answer2          = inputdlg(prompt,dlgtitle,[1 70],definput,opts);
                            
                            l_p   = str2double(answer2{1});
                            ri_p  = str2double(answer2{2});
                            rf_p  = str2double(answer2{3});
                                                     
                            syms X1
                            r_sym = ri_p+X1*(rf_p-ri_p);

                            if has(r_sym,X1)
                                r_fn = matlabFunction(r_sym);
                            else
                                r_fn = str2func(['@(X1)' num2str(ri_p)]);
                            end

                            r{i} = r_fn;
                            
                            if i==1
                                A0 = pi*ri_p^2;
                            end
                            
                        elseif CS=='R'
                            
                            prompt           = {'Length (m):','Initial height (m):','Initial width (m):','Final height (m):','Final width (m):'};
                            dlgtitle         = ['Enter geometric properties of division ',num2str(i)];
                            definput         = {'0.3','0.02','0.02','0.02','0.02'};
                            opts.Interpreter = 'tex';
                            answer2 = inputdlg(prompt,dlgtitle,[1 70],definput,opts);
                            
                            l_p    = str2double(answer2{1});
                            hi_p   = str2double(answer2{2});
                            wi_p   = str2double(answer2{3});
                            hf_p   = str2double(answer2{4});
                            wf_p   = str2double(answer2{5});
                            
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
                            if i==1
                                A0 = hi_p*wi_p;
                            end
                                                        
                         elseif CS=='E'
                            
                            prompt           = {'Length (m):','Initial semi-major axis (m):','Initial semi-minor axis (m):','Final semi-major axis (m):','Final semi-minor axis (m):'};
                            dlgtitle         = ['Enter geometric properties of division ',num2str(i)];
                            definput         = {'0.3','0.04','0.02','0.04','0.02'};
                            opts.Interpreter = 'tex';
                            answer2 = inputdlg(prompt,dlgtitle,[1 70],definput,opts);
                            
                            l_p    = str2double(answer2{1});
                            ai_p   = str2double(answer2{2});
                            bi_p   = str2double(answer2{3});
                            af_p   = str2double(answer2{4});
                            bf_p   = str2double(answer2{5});
                            
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
                            if i==1
                                A0 = pi*ai_p*bi_p;
                            end
                        end
                        
                        lp{i}     = l_p;
                        gi{i}     = eye(4);
                        gf{i}     = eye(4);
                        
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
                    Lscale = (A0*Li.L).^(1/3);

            end

                Li.n_l  = 25;
                if Li.CS~='R'
                    Li.n_r = 18;
                end
                Li.color  = rand(1,3);
                Li.CPF    = false;
                Li.PlotFn = @(g) CustomShapePlot(g);
                Li.Lscale = Lscale;
                
        end
    end
    
    methods
        %% Set and property update functions
        
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
        
        %Changing E:
        function set.E(SorosimLink, value)
            SorosimLink.E = value;
            SorosimLink.UpdateG();
        end
        
        %Changing Poi:
        function set.Poi(SorosimLink, value)
            SorosimLink.Poi = value;
            SorosimLink.UpdateG();
        end
        
        %Changing lp:
        function set.lp(SorosimLink, value)
            SorosimLink.lp = value;
            SorosimLink.UpdateL();
        end
        
        function set.jointtype(SorosimLink, value)
            SorosimLink.jointtype = value;
            SorosimLink.UpdateKD();
        end
        
        %Updating the properties
        function Update(SorosimLink)
            if isempty(SorosimLink.color)
                return
            end
            if SorosimLink.linktype=='r'
                SorosimLink.M = LinkPropUpdate(SorosimLink);
            end
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
            if isempty(SorosimLink.color)
                return
            end
            SorosimLink.L = sum(cell2mat(SorosimLink.lp));
        end
        
        function UpdateKD(SorosimLink)
            if isempty(SorosimLink.color)
                return
            end
            SorosimLink.Kj = 0;
            SorosimLink.Dj = 0;
        end
        
        
    end
    
end