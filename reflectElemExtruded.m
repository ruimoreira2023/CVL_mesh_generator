function [E]=reflectElemExtruded(E0,nelem,nnode)
% nelem is the nof surface module elements


for ielem = 1:nelem/2
    %n1, n2, n3, n4
    E0(ielem+nelem/2,8:11) = E0(ielem,8:11)+nnode/2;
end

E=E0;