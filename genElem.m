function [E,numE,numN]=genElem(E0,m,n,posE,posN, mat,prop)
% E0 - element connectivity matrix
% m - number of elements in x-direction
% n - number of elements in y-direction
% posE - element position/number
% posN - node position/number
% mat - material
% prop - property
% Returns complete element connectivity matrix and
% the number of the element and nodes

for row=1:n
    for col=1:m
        ielem=(row-1)*m+col; % index of elements
        % 4-node elements (n1,2,3,4)
        n1=(row-1)*m+col+row-1+posN;
        n2=(row-1)*m+col+1+row-1+posN;
        n3=row*m+2+col+row-1+posN;
        n4=row*m+1+col+row-1+posN;
        % Regenerates the nodes index of ELEM matrix
        E0(ielem+posE,2:7)=[n1 n2 n3 n4 mat prop];
    end
end

numE=posE+m*n;
numN=posN+(m+1)*(n+1);
E=E0;