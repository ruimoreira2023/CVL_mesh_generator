function [NODES2,ELEM] = CheckCoincidentNodes(NODES, ELEM, error)

% verifies if nodes in NODES are coincident (if distance between a node
% pair is <= "error" value



nnode=size(NODES,1); % returns number of rows of NODES
nelem=size(ELEM,1); % returns number of rows of ELEM

for inode=1:nnode
    for jnode=inode+1:nnode
        ori=NODES(inode,2:4);
        dup=NODES(jnode,2:4);

        if abs(ori-dup)<=error %these nodes are coincident
            
            for ielem=1:nelem
                for celem=2:5
                    if ELEM(ielem,celem)==NODES(jnode,1)
                        ELEM(ielem,celem)=NODES(inode,1);
                    end
                end
%                 for celem=8:11
%                     if ELEM(ielem,celem)==NODES(jnode,1)
%                         ELEM(ielem,celem)=NODES(inode,1);
%                     end
%                 end
            end
            NODES(jnode,1)= 0 ;%NODES(inode,1);
        else         
        end
        
    end
end

% remove zero lines from NODES
idx_nonzerolines = NODES(:,1)>0 ;
NODES2 = NODES(idx_nonzerolines,:) ;


end