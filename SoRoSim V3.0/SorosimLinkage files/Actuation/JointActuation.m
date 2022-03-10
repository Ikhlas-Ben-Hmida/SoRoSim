%Function to identify if a joint is actuated either by torque or joint
%coordintates. Outputs are the number of joint actuators, indices of links at
%which the joint is actuated, corresponding indices in Qspace, variable to
%save if the joint is wrench controlled or not, and the basis for 1 dof
%joints. 
%Last modified by Anup Teejo Mathew - 23/05/2021

function [n_jact,i_jact,i_jactq,WrenchControlled,Bqj1] = JointActuation(Tr)

    if any(Tr.jointtype=='R')%revolute joint
        [n_Ract,i_Ract,i_Ractq,WrenchControlledR,BqR] = RevoluteJointActuation(Tr);
    else
        n_Ract            = 0;
        i_Ract            = [];
        i_Ractq           = [];
        WrenchControlledR = [];
        BqR               = [];
    end

    n_jact             = n_Ract;
    i_jact             = i_Ract;
    i_jactq            = i_Ractq;
    WrenchControlled   = WrenchControlledR;
    Bqj1               = BqR;


    %Prismatic joint
    if any(Tr.jointtype=='P')
        [n_Pact,i_Pact,i_Pactq,WrenchControlledP,BqP] = PrismaticJointActuation(Tr);
    else
        n_Pact            = 0;
        i_Pact            = [];
        i_Pactq           = [];
        WrenchControlledP = [];
        BqP               = [];
    end

    n_jact              = n_jact+n_Pact;
    i_jact              = [i_jact i_Pact];
    i_jactq             = [i_jactq i_Pactq];
    WrenchControlled    = [WrenchControlled WrenchControlledP];
    Bqj1                = [Bqj1 BqP];

    %Helical joint
    if any(Tr.jointtype=='H')
        [n_Hact,i_Hact,i_Hactq,WrenchControlledH,BqH] = HelicalJointActuation(Tr);
    else
        n_Hact            = 0;
        i_Hact            = [];
        i_Hactq           = [];
        WrenchControlledH = [];
        BqH               = [];
    end

    n_jact              = n_jact+n_Hact;
    i_jact              = [i_jact i_Hact];
    i_jactq             = [i_jactq i_Hactq];
    WrenchControlled    = [WrenchControlled WrenchControlledH];
    Bqj1                = [Bqj1 BqH];

    %Universal joint
    if any(Tr.jointtype=='U')
        [n_Uact,i_Uact,i_Uactq,WrenchControlledU] = UniversalJointActuation(Tr);
    else
        n_Uact            = 0;
        i_Uact            = [];
        i_Uactq           = [];
        WrenchControlledU = [];
    end

    n_jact              = n_jact+n_Uact;
    i_jact              = [i_jact i_Uact];
    i_jactq             = [i_jactq i_Uactq];
    WrenchControlled    = [WrenchControlled WrenchControlledU];

    %Cylindrical joint
    if any(Tr.jointtype=='C')
        [n_Cact,i_Cact,i_Cactq,WrenchControlledC] = CylindricalJointActuation(Tr);
    else
        n_Cact            = 0;
        i_Cact            = [];
        i_Cactq           = [];
        WrenchControlledC = [];
    end

    n_jact              = n_jact+n_Cact;
    i_jact              = [i_jact i_Cact];
    i_jactq             = [i_jactq i_Cactq];
    WrenchControlled    = [WrenchControlled WrenchControlledC];

    %Planar joint 
    if any(Tr.jointtype=='A')
        [n_Aact,i_Aact,i_Aactq,WrenchControlledA] = PlanarJointActuation(Tr);
    else
        n_Aact            = 0;
        i_Aact            = [];
        i_Aactq           = [];
        WrenchControlledA = [];
    end

    n_jact              = n_jact+n_Aact;
    i_jact              = [i_jact i_Aact];
    i_jactq             = [i_jactq i_Aactq];
    WrenchControlled    = [WrenchControlled WrenchControlledA];

    %Spherical joint 
    if any(Tr.jointtype=='S')
        [n_Sact,i_Sact,i_Sactq,WrenchControlledS] = SphericalJointActuation(Tr);
    else
        n_Sact            = 0;
        i_Sact            = [];
        i_Sactq           = [];
        WrenchControlledS = [];
    end

    n_jact              = n_jact+n_Sact;
    i_jact              = [i_jact i_Sact];
    i_jactq             = [i_jactq i_Sactq];
    WrenchControlled    = [WrenchControlled WrenchControlledS];

    %Free joint 
    if any(Tr.jointtype=='F')
        [n_Fact,i_Fact,i_Factq,WrenchControlledF] = FreeJointActuation(Tr);
    else
        n_Fact            = 0;
        i_Fact            = [];
        i_Factq           = [];
        WrenchControlledF = [];
    end

    n_jact              = n_jact+n_Fact;
    i_jact              = [i_jact i_Fact];
    i_jactq             = [i_jactq i_Factq];
    WrenchControlled    = [WrenchControlled WrenchControlledF];
end

