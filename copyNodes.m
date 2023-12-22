function [NODES]=copyNodes(N0,nnodeExtrude,numCopies,mod)

nodeID = nnodeExtrude+1:numCopies*nnodeExtrude;
N0(nnodeExtrude+1:numCopies*nnodeExtrude,1) = nodeID;

for idx = 1:numCopies-1
    N0(idx*nnodeExtrude+1:(idx+1)*nnodeExtrude,2) = N0(1:nnodeExtrude,2) + idx*2*mod;
end

M = numCopies-1; % Ns
N = 1;
N0(nnodeExtrude+1:numCopies*nnodeExtrude,3) = repmat(N0(1:nnodeExtrude,3), [M,N]);
N0(nnodeExtrude+1:numCopies*nnodeExtrude,4) = repmat(N0(1:nnodeExtrude,4), [M,N]);


NODES = N0;