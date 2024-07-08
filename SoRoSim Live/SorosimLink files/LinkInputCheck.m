function badanswer = LinkInputCheck(jointtype,CS,Kj,badanswer)

if ~any(strcmp(jointtype,{'R','P','H','C','A','S','F','N'}))
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