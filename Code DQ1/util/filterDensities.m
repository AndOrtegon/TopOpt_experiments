function xFil = filterDensities(x)
global ft He Hs rmin nelx nely problem

if ft == 1
    xFil = x;
elseif ft == 2
%     xFil(:) = (H*x(:))./Hs;
    xFil = He*x(:) ;
    xFil = reshape(xFil,size(x)) ;
end
end