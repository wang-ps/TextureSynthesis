function [ zi, zj ] = Idx2Coordinate( idx, w, n )

 	zj = mod(idx, n-2*w);
	if zj == 0
        zj = w;
    end
    zi = fix((idx-zj)/(n-2*w) + 1);
    zi = zi + w;
    zj = zj + w;
end

