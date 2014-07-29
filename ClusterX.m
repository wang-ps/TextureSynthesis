function [ Xc, cp ] = ClusterX( X, w, num )


	[m, n, c] = size(X);
	% image data
	XN = zeros((m-2*w)*(n-2*w), k);
    for i = 1 : m-2*w
        for j = 1 : n-2*w
            idx = (i-1)*(n-2*w) + j;
            XN(idx, :) = reshape(Z(i:i+2*w, j:j+2*w, :), 1, k);
        end
    end

    [IDX, ~, ~, D] = kmeans(X, num);
    nidx = size(IDX, 1);
    Xc = zeros(num, 2);
    distance = 1.0e30 * ones(num, 1);
    for i = 1 : nidx
    	d = D(i, IDX(i));
    	if d < distance(IDX(i))
    		distance(IDX(i)) = d;
    		Xc(IDX(i), :) = Idx2Coordinate(i, w, n);
    	end
    end

    cp = zeros(m, n, 2);
    for i = 1 : m-2*w
        for j = 1 : n-2*w
        	ii =  (i-1)*(n-2*w) + j;
    		cp(i, j, :) = Xc(IDX(ii));
    	end
    end
end
