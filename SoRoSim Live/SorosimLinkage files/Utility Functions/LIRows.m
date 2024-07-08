%Code to get linearly independent rows and their indices from a matrix
%using QR decomposition technique
%Modified the code 'licols'
%Last modified by Anup Teejo Mathew 02.03.2022

function [Xsub,idx]=LIRows(X,tol)
%Extract a linearly independent set of rows of a given matrix X
%
%    [Xsub,idx]=licols(X)
%
%in:
%
%  X: The given input matrix
%  tol: A rank estimation tolerance. Default=1e-10
%
%out:
%
% Xsub: The extracted columns of X
% idx:  The indices (into X) of the extracted columns
X=X';
     if ~nnz(X) %X has no non-zeros and hence no independent columns
         Xsub=[]; idx=[];
         return
     end
     if nargin<2, tol=1e-6; end
       [Q, R, E] = qr(X,0); 
       if ~isvector(R)
        diagr = abs(diag(R));
       else
        diagr = R(1);   
       end
       %Rank estimation
       r    = find(diagr >= tol*diagr(1), 1, 'last'); %rank estimation
       idx  = sort(E(1:r));
       Xsub = X(:,idx);
       Xsub = Xsub';
end