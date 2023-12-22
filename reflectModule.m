function [XmReflected,YmReflected]=reflectModule(Xm,Ym,module)
% Given the coordinates matrix (Xm, Ym)
% This function returns de reflected coordinates across the module plane
% Module is the length of half geometry

XmReflected = zeros(size(Xm,1), size(Xm,2));
YmReflected = zeros(size(Ym,1), size(Ym,2));
distToMod = zeros(size(Xm,1), size(Xm,2));

for i = 1:size(Xm,1)
    for j = 1:size(Xm,2)
        distToMod(i,j) = module - Xm(i,j);
        XmReflected(i,j) = Xm(i,j) + 2*distToMod(i,j);
        YmReflected(i,j) = Ym(i,j);
    end
end
