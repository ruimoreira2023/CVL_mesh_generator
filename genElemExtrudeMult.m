function [E]=genElemExtrudeMult(E0,nz,nelem,nnode)
% E0 - element connectivity matrix
% m - number of elements in x-direction
% n - number of elements in y-direction
% posE - element position/number
% posN - node position/number
% mat - material
% prop - property
% Returns complete element connectivity matrix and
% the number of the element and nodes

% genElemExtrude multiple times
% nelem of surface only

% index = 0:(nz-1); % [0:9]
index = 1:(nz-1); % [1:9]

for layer = 1:nz-1
    for ielem = 1:nelem
        % First 1:9
        % Second 0:8
%         E0(ielem+(index(layer)+1)*nelem,2:5) = E0(ielem+index(layer)*nelem,8:11);
        E0(ielem+index(layer)*nelem,2:5) = E0(ielem+(index(layer)-1)*nelem,8:11);
        E0(ielem+index(layer)*nelem,8:11) = E0(ielem+(index(layer)-1)*nelem,8:11)+nnode;
        
        E0(ielem+(index(layer))*nelem,6:7) = E0(ielem,6:7);
        
        %         E0(ielem+1*nelem,2:5) = E0(ielem+0*nelem,8:11);
        %         E0(ielem+2*nelem,2:5) = E0(ielem+1*nelem,8:11);
        %         E0(ielem+3*nelem,2:5) = E0(ielem+2*nelem,8:11);
    end
end


E = E0;