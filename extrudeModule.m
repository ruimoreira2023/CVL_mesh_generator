function [Zm]=extrudeModule(Xm,extrudeLength,nz)

Zm = zeros(size(Xm,1),size(Xm,2),nz);

for izs = 2:nz+1
    Zm(:,:,izs) = -(izs-1)*(extrudeLength/nz);
end
