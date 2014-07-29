function [ i, j ] = Idx2Coordinate( idx, w, n )

	j = mod(idx, n-2*w);
    if j == 0
        j = w;
    end
    i = fix((idx-j)/(n-2*w) + 1);
    i = i+w;
    j = j+w;

end

