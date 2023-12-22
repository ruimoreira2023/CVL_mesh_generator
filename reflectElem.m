function [E]=reflectElem(E0,nelem,nnode)
% E0 - element connectivity matrix
% nelem - number of surface module elements
% nnode - number of surface module nodes

for ielem = 1:nelem/2
    %n1, n2, n3, n4
    E0(ielem+nelem/2,2:5) = E0(ielem,2:5)+nnode/2;
    E0(ielem+nelem/2,6:7) = E0(ielem,6:7);
end

E=E0;