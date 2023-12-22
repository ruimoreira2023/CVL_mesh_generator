function [E,numE,numN]=genElemExtrude(E0,m,n,posE,posN,nnode)
% E0 - element connectivity matrix
% m - number of elements in x-direction
% n - number of elements in y-direction
% posE - element position/number
% posN - node position/number
% nnode - 
% mat - material
% prop - property
% Returns complete element connectivity matrix and
% the number of the element and nodes

% genElemExtrude
% nelem of surface only


for row=1:n
    for col=1:m
        ielem=(row-1)*m+col; % index of elements
        % 4-node elements (n1,2,3,4)
        n1=(row-1)*m+col+row-1+posN; n5=n1+nnode;
        n2=(row-1)*m+col+1+row-1+posN; n6=n2+nnode;
        n3=row*m+2+col+row-1+posN; n7=n3+nnode;
        n4=row*m+1+col+row-1+posN; n8=n4+nnode;
        % Regenerates the nodes index of ELEM matrix
        % ELEM=[ielem node1 node2 node3 node4 node5 node6 node7 node8 Prop Mat];
        E0(ielem+posE,8:11)=[n5 n6 n7 n8];
    end
end

numE=posE+m*n;
numN=posN+(m+1)*(n+1);
E=E0;