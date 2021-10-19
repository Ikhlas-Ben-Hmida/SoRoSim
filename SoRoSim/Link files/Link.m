%Class to define a link
%A Link consists of a rigid joint and a body.
%Last modified by Anup Teejo Mathew - 20/05/2021
classdef Link < handle
    
    properties
        
        %General Properties
        
        jointtype       %Type of joint used to connect the link (lumped DoF). (R) for Revolute,(P) for Prismatic, (H) for Helical, (U) for Universal, (C) for Cylindrical, (A) for Planar, (S) for Spherical, (F) for Free motion and (N) for Fixed
        linktype        %'s' for soft or 'r' for rigid
        CS              %Cross section: 'R' for rectangular 'C' for circular
        npie            %Number of pieces. For a rigid link, npie = 1. For a soft link, npie=1+number of divisions 
        nGauss          %Number of Gaussian Quadrature points + 2
        Xs              %Gaussian Quadrature points inc 0 and 1
        Ws              %Gaussian Quadrature weights
        
        
        %Geometric Properties
        
        lp        %Length of each divisions of the link (soft link) [m]
        L         %Total length of the link [m]
        r         %Radius as function of X1 (X1=X/L, X1 varies from 0 to 1) [m]
        h         %Height as function of X1 (X1=X/L, X1 varies from 0 to 1) [m]
        w         %Width as function of X1 (X1=X/L, X1 varies from 0 to 1) [m]
        cx        %x coordinate of origin wrt. previous frame [m]
        cy        %y coordinate of origin wrt. previous frame [m]
        cz        %z coordinate of origin wrt. previous frame [m]
        
        %Material
        
        E         %Young's modulus [Pa]
        Poi       %Poisson's ratio [-]
        G         %Shear modulus [Pa]
        Eta       %Material Damping [Pa.s]
        Rho       %Density [kg/m^3]
        
        Ms        %Screw inertia
        Es        %Screw stiffness
        Gs        %Screw Damping
        
        %Plot Properties
        
        n_l       %Number of cross sections per division.
        n_r       %Number of radial points if the cross section is circular
        color     %color of link
        
    end
    
    
    %%
    methods
        function Link = Link
            
            quest = 'Select link type:';
            answ  = questdlg(quest,'Link type',...
                'Soft','Rigid','Soft');
            
            switch answ
                case 'Rigid'
                    
                    Link.linktype='r';
                    
                    prompt           = {'Joint type (Enter (R) for Revolute,(P) for Prismatic, (H) for Helical, (U) for Universal, (C) for Cylindrical, (A) for Planar, (S) for Spherical, (F) for Free motion and (N) for Fixed): ',...
                                        'Density (kg/m^3):','Cross-section shape (C for circular, R for rectangular):','y coordinate of origin (m) wrt. previous frame:','z coordinate of origin (m) wrt. previous frame:','Length of rigid link (m):'};
                    dlgtitle         = 'Rigid Link Properties';
                    definput         = {'N','1000','C','0','0','0.3'};
                    opts.Interpreter = 'tex';
                    answer           = inputdlg(prompt,dlgtitle,[1 60],definput,opts);
                    
                    Link.jointtype = answer{1};
                    Link.npie      = 1;
                    Rho            = str2double(answer{2});
                    Link.Rho       = Rho;
                    CS             = answer{3};
                    Link.CS        = CS;
                    Link.cy        = str2double(answer{4});
                    Link.cz        = str2double(answer{5});
                    L              = str2double(answer{6});
                    Link.L         = L;
                    
                    if CS=='C'
                        
                        prompt           = {'Initial radius (m):','Final radius (m):'};
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
                        
                        Link.r = r;
                        
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
                        
                        Iy  = 0;
                        for ii=2:nGauss-1
                            r_here = r(Xs(ii));
                            Iy     = Iy+Ws(ii)*L*((pi/4)*r_here^4+pi*r_here^2*(Xs(ii)*L-cx)^2);
                        end
                        Iz = Iy;
                        Ix = Iy+Iz;
                        
                        
                    elseif CS=='R'
                        
                        prompt           = {'Initial height (m):','Initial width (m):','Final height (m):','Final width (m):'};
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
                        
                        Link.h = h;
                        Link.w = w;
                        
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
                        
                        Iy  = 0;
                        Iz  = 0;
                        for ii=2:nGauss-1
                            h_here = h(Xs(ii));
                            w_here = w(Xs(ii));
                            Iy     = Iy+Ws(ii)*L*((Xs(ii)*L-cx)^2*h_here*w_here+h_here*w_here^3/12);
                            Iz     = Iz+Ws(ii)*L*((Xs(ii)*L-cx)^2*h_here*w_here+h_here^3*w_here/12);
                        end
                        Ix = Iy+Iz;
                    end
                    
                    
                    
                    Ms   = double([Ix 0 0 0 0 0;0 Iy 0 0 0 0;0 0 Iz 0 0 0;...
                        0 0 0 mass 0 0;0 0 0 0 mass 0;0 0 0 0 0 mass]);
                    
                    Link.cx = cx;
                    Link.Ms = Ms;
                    
                    %%
                case 'Soft'
                    
                    Link.linktype = 's';
                    
                    prompt           = {'Joint type (Enter (R) for Revolute,(P) for Prismatic, (H) for Helical, (U) for Universal, (C) for Cylindrical, (A) for Planar, (S) for Spherical, (F) for Free motion and (N) for Fixed):',...
                        'Number of divisions:','Cross-section shape (C for circular, R for rectangular):'...
                        ,'Density (kg/m^3):','Youngs Modulus (N/m^2)','Poissons Ratio:','Material Damping (Pa.s)(for dynamic simulation):'};
                    dlgtitle         = 'Soft Link Properties';
                    
                    definput         = {'N','1','C','1000','1e6','0.5','0.112e5'};
                    opts.Interpreter = 'tex';
                    answer           = inputdlg(prompt,dlgtitle,[1 60],definput,opts);
                    
                    Link.jointtype  = answer{1};
                    ndiv            = str2double(answer{2});
                    Link.npie       = ndiv+1;
                    CS              = answer{3};
                    Link.CS         = CS;
                    Rho             = str2double(answer{4});
                    Link.Rho        = Rho;
                    E               = str2double(answer{5});
                    Link.E          = E;
                    Poi             = str2double(answer{6});
                    Link.Poi        = Poi;
                    Eta             = str2double(answer{7});
                    Link.Eta        = Eta;
                    G               = E/(2*(1+Poi));
                    Link.G          = G;
                    
                    lp   = cell(1,ndiv);
                    r    = cell(1,ndiv);
                    h    = cell(1,ndiv); w  = cell(1,ndiv);
                    Ms   = cell(1,ndiv); Es = cell(1,ndiv); Gs = cell(1,ndiv);
                    cy   = cell(1,ndiv); cz = cell(1,ndiv);
                    
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
                        cy{i}     = str2double(answer2{3});
                        cz{i}     = str2double(answer2{4});
                        Xs{i}     = Xs_p;
                        Ws{i}     = Ws_p;
                        nGauss{i} = nGauss_p;
                        
                    end
                    
                    Link.L=sum(cell2mat(lp));
                    Link.lp=lp;
                    Link.r=r;
                    Link.w=w;
                    Link.h=h;
                    Link.cy=cy;
                    Link.cz=cz;
                    Link.Ms=Ms;
                    Link.Es=Es;
                    Link.Gs=Gs;
                    Link.Xs=Xs;
                    Link.Ws=Ws;
                    Link.nGauss=nGauss;

            end

                Link.n_l  = 10;
                if Link.CS=='C'
                    Link.n_r = 18;
                end
                Link.color = 'b';

        end
    end
    
    methods
        %% Set and property update functions
        
        %Changing Poi:
        function set.Poi(Link, value)
            Link.Poi = value;
            Link.UpdateG();
            Link.Update()
        end
        
        %Changing E:
        function set.E(Link, value)
            Link.E = value;
            Link.UpdateG();
            Link.Update()
        end
        
        %Changing Eta:
        function set.Eta(Link, value)
            Link.Eta = value;
            Link.Update()
        end
        
        %Changing G:
        function set.G(Link, value)
            Link.G = value;
            Link.Update()
        end
        
        %Changing r:
        function set.r(Link, value)
            Link.r = value;
            Link.Update()
        end
        
        %Changing w:
        function set.w(Link, value)
            Link.w = value;
            Link.Update()
        end

        %Changing h:
        function set.h(Link, value)
            Link.h = value;
            Link.Update()
        end
        
        %Changing L: (rigid or constant CS soft)
        function set.L(Link, value)
            Link.L = value;
            Link.Update()
        end
        
        %Changing Rho:
        function set.Rho(Link, value)
            Link.Rho = value;
            Link.Update()
        end

        %Changing lp:
        function set.lp(Link, value)
            Link.lp = value;
            Link.Update()
            Link.UpdateL()
        end
        
        %Changing nGauss:
        function set.nGauss(Link, value)
            Link.nGauss = value;
            Link.UpdateNGauss();
        end
        
        %Updating the properties
        function Update(Link)
            if isempty(Link.color)
                return
            end
            if Link.linktype=='s'
                [Link.Ms,Link.Es,Link.Gs] = PropUpdate(Link);
            else
                Link.Ms = PropUpdate(Link);
            end
        end
        
        %Update nGauss
        function UpdateNGauss(Link)
            if isempty(Link.color)
                return
            end
            for i=1:Link.npie-1
                [Xs_p,Ws_p] = GaussQuadrature(Link.nGauss{i}-2);
                Link.Xs{i}  = Xs_p;
                Link.Ws{i}  = Ws_p;
            end
            [Link.Ms,Link.Es,Link.Gs] = PropUpdate(Link);
        end
        
        %Update E or Poi
        function UpdateG(Link)
            if isempty(Link.color)
                return
            end
            Link.G=Link.E/(2*(1+Link.Poi));
        end
        
        %Update lp
        function UpdateL(Link)
            Link.L = sum(cell2mat(Link.lp));
        end
        
        
    end
    
end