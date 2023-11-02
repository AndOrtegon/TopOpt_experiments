function xp = ord2pos(x)
    xp = zeros(2*size(x,2),2*size(x,3)) ;
    xp(2:2:end,1:2:end) = x(1,:,:) ;
    xp(2:2:end,2:2:end) = x(2,:,:) ;
    xp(1:2:end,2:2:end) = x(3,:,:) ;
    xp(1:2:end,1:2:end) = x(4,:,:) ;
end