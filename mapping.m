function [X,Y] = mapping(Xi,Yi,m,n)
% Xi, Yi - coordinates of nodes
% m - number of elements of the region/x dim (A,B,C)
% n - number of elements of the layer/y dim (1,2,3)
% Returns coordinates of the nodes in each region/layer line

% Region map
% Generates a vector with m+1 number of nodes in the x-direction
% equally spaced
t = linspace(0,1,m+1);

% Generate grid 
% Alocates nodes (x,y) in the x-direction in line 2-4 and line 1-3
L1 = [Xi(2)+(Xi(4)-Xi(2))*t ; Yi(2)+(Yi(4)-Yi(2))*t] ;
L2 = [Xi(1)+(Xi(3)-Xi(1))*t ; Yi(1)+(Yi(3)-Yi(1))*t] ;

% Mesh
X = zeros(m+1,n+1) ;
Y = zeros(m+1,n+1)  ;
% Generates a vector with n+1 number of nodes in the y-direction
% equally spaced
t = linspace(0,1,n+1) ;

% Now, it alocates nodes in the y-direction for each mi node
for i = 1:m+1
    X(i,:) = L1(1,i)+(L2(1,i)-L1(1,i))*t ;
    Y(i,:) = L1(2,i)+(L2(2,i)-L1(2,i))*t ;
end

end