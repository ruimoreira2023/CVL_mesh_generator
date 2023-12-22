function [ELEM]=copyElem(E0,nelemExtrude,nnodeExtrude,numCopies)

elemID = nelemExtrude+1:numCopies*nelemExtrude;
E0(nelemExtrude+1:numCopies*nelemExtrude,1) = elemID;

for idx = 1:numCopies-1
    for ielem = 1:nelemExtrude
        
        E0(idx*nelemExtrude+1:(idx+1)*nelemExtrude,2:5) = E0(1:nelemExtrude,2:5) + idx*nnodeExtrude;
        E0(idx*nelemExtrude+1:(idx+1)*nelemExtrude,8:11) = E0(1:nelemExtrude,8:11) + idx*nnodeExtrude;

    end
end

% mat & prop
M = numCopies-1; % Ns
N = 1;
E0(nelemExtrude+1:numCopies*nelemExtrude,6:7) = repmat(E0(1:nelemExtrude,6:7), [M,N]);

ELEM = E0;