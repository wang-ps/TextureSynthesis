function [ Sq ] = CoherentSet( X, w, num)
%COHERENTSET Summary of this function goes here
%   Detailed explanation goes here
	[m, n, c] = size(X);
	k = c*(2*w+1)^2;

	%h = fspecial('Gaussian');
	% % build Gaussian Pyramid
	% % assert level >= 1
	% Xg = cell(1, level);
	% Xg{1} = X;
	% for i = 2 : level
	% 	X = Xg{i-1};
	% 	X = imfilter(X, h, 'replicate');
	% 	Xg{i} = X(1:2:end, 1:2:end, :);
	% end
	% for i = level : -1 : 1
	% 	% build kdt		
	% end

	% build kdt
 	% nearest neighbour data
    XN = zeros((m-2*w)*(n-2*w), k);
    for i = 1 : m-2*w
        for j = 1 : n-2*w
            idx = (i-1)*(n-2*w) + j;
            XN(idx, :) = reshape(Z(i:i+2*w, j:j+2*w, :), 1, k);
        end
    end
    %create kdtree
    kdt = createns(XN,'nsmethod','kdtree');

    Sq = cell(m, n);

    for i = w+1 : m-w
    	for j = w+1 : n-w
            xv = reshape(X(i-w:i+w, j-w:j+w, :), 1, k);
            idx = knnsearch(kdt, xp, 'k', num);
            nidx = size(idx, 1);
            zp = [];
            for ii = 1 : nidx
            	zj = mod(idx(ii), nz-2*w);
                if zj == 0
                    zj = w;
                end
                zi = fix((idx(ii)-zj)/(nz-2*w) + 1);
                zp = [zp; [zi+w zj+w]];
            end
            Sq{i, j} = zp;    
        end       
    end

end

