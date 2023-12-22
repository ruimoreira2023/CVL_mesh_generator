function [Zm]=coordToColumn(Zm,nz)

nlines = size(Zm,1)*size(Zm,2);

Zm = reshape(Zm,nlines,nz+1);