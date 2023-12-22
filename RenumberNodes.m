function [NODES,ELEM] = RenumberNodes(NODES, ELEM, criteria)

% criteria:
% 'X+' - sort according to x-coordinate (lowest to highest)
% 'X-' -                                  (highest to lowest)
% 'Y+' - sort according to y-coordinate (lowest to highest)    
% 'Y-' -                                  (highest to lowest)
% 'Z+' - sort according to z-coordinate (lowest to highest)
% 'Z-' -                                  (highest to lowest)
% '+'  - sort based on node id number 1-->++
% '-'  - sort based on node id number ++->1


shift_index=max(NODES(:,1))+10; 

switch upper(criteria)
    case 'X+'
        NODES=sortrows(NODES,[2 3 4]);
    case 'X-'
        NODES=sortrows(NODES,[-2 3 4]);
    case 'Y+'
        NODES=sortrows(NODES,[ 3 2 4]);
    case 'Y-'
        NODES=sortrows(NODES,[-3 2 4]);
    case 'Z+'
        NODES=sortrows(NODES,[4 2 3]);
    case 'Z-'
        NODES=sortrows(NODES,[-4 2 3]); 
    case '+'
        NODES=sortrows(NODES,1);
    case '-'
        NODES=sortrows(NODES,-1);
end

NodesInELEM=ELEM(:,2:5);

for inode=1:size(NODES,1)
    [index]=find(NodesInELEM==NODES(inode));
    ELEM(index+size(ELEM,1))=inode;
end


NODES(:,1)=1:size(NODES,1);
