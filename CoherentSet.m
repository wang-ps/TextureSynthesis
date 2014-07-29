function [ Sq ] = CoherentSet( X, w, num)

	[m, n, c] = size(X);
    X = double(X);
	k = c*(2*w+1)^2;
    
	% build kdt
 	% nearest neighbour data
    XN = zeros((m-2*w)*(n-2*w), k);
    for i = 1 : m-2*w
        for j = 1 : n-2*w
            idx = (i-1)*(n-2*w) + j;
            XN(idx, :) = reshape(X(i:i+2*w, j:j+2*w, :), 1, k);
        end
    end
    %create kdtree
    kdt = createns(XN,'nsmethod','kdtree');

    Sq = cell(m, n);
    for i = w+1 : m-w
    	for j = w+1 : n-w
            xv = reshape(X(i-w:i+w, j-w:j+w, :), 1, k);
            idx = knnsearch(kdt, xv, 'k', num);
            zp = ones(num, 2);
            for ii = 1 : num
            	[zi, zj] = Idx2Coordinate(idx(ii), w, n);
                zp(ii, :) = [zi, zj];
            end
            Sq{i, j} = zp;    
        end       
    end

end

