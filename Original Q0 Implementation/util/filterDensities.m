function xFil = filterDensities(x)
    global ft He
    if ft == 1
        xFil = x;
    elseif ft == 2
        xFil = He*x(:) ;
        xFil = reshape(xFil,size(x)) ;
    end
end