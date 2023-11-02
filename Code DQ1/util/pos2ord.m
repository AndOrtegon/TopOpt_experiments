function x = pos2ord(xp) 
    x = zeros(4,size(xp,1)/2,size(xp,2)/2) ;
    x(1,:,:) = xp(2:2:end,1:2:end) ;
    x(2,:,:) = xp(2:2:end,2:2:end) ;
    x(3,:,:) = xp(1:2:end,2:2:end) ;
    x(4,:,:) = xp(1:2:end,1:2:end) ;
end