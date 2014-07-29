function [ Z ] = InverseTextureSynthesis(X, w, sq, Xc, cp, Z)
% Inverse Texture Synthesis

    
	[mx, nx, c] = size(X);
	[mz, nz, c] = size(Z);
    k = c*(2*w+1)^2;
    X = double(X);
    Z = double(Z);

	alpha = 0.01;
	sample_rate_inv = ceil(w/2);
	sample_rate_for = ceil(w/4);
	weight_inv = sample_rate_inv^2 /(mx*nx)  ;
	weight_for = alpha *  sample_rate_for^2 /(mz*nz);

	nPixel = c*(2*w+1)^2; % #pixel in a window
	%Zp = zeros(mx, nx, 2); % Xp -> Zp, inverse item
	%Xq = zeros(mz, nz, 2); % Zq -> Xq, forward item
	Zp(:, :, 1) = randi([w+1, mz-w], mx, nx);
	Zp(:, :, 2) = randi([w+1, nz-w], mx, nx);
	Xq(:, :, 1) = randi([w+1, mx-w], mz, nz);
	Xq(:, :, 2) = randi([w+1, nx-w], mz, nz);

% 	% calculate c(p)
% 	cp = zeros(mx, nx);
% 	nXc = size(Xc, 1);
% 	Xcv = size(nXc, nPixel);
% 	for i = 1 : nXc
% 		[ixc, jxc] = Xc(i);
% 		Xcv(i, :) = reshape(X(ixc-w:ixc+w, jxc-w:jxc+w, :), 1, nPixel);
% 	end
% 	for i = 1 : mx-2*w
% 		for j = 1 : nx-2*w
% 			Xv = reshape(X(i-w:i+w, j-w:jxc+w, :), 1, nPixel);
% 		end
% 	end

    iter_num = 10;
	for it = 1:iter_num	
		%% z E-step
		% calc k-coherence set
% 		KCHS = cell(mz, nz);
% 		dis(:, :, 1) = diag(-w:w)*ones(2*w+1);
% 		dis(:, :, 2) = ones(2*w+1)*diag(-w:w);	% displacement
% 
% 		for i = w+1 : mz-w
% 			for j = w+1 : nz-w						
% 				chp = Xq(i-w:i+w, j-w:j+w, :) - dis; 	% coherent pixel
% 				chs = [];								% coherent set
% 				for ic = 1 : 2*w+1
% 					for jc = 1 : 2*w+1   
% 						chs = [chs; sq{chp(ic, jc, 1), chp(ic, jc, 2)}];
% 					end
% 				end
% 				KCHS{i, j} = unique(chs, 'rows');	
% 			end 
% 		end		

        KCHS  = KCoherentSet( mz, nz, w, Xq, sq );
        
		ZX1 = cell(mz, nz);
		% for the forward item
		for i = w+1 : sample_rate_for : mz-w
			for j = w+1 : sample_rate_for : nz-w
				for iw = -w : w
					for jw = -w : w
						ZX1{i + iw, j+jw} = [ZX1{i + iw, j+jw}; Xq(i, j) + [iw, jw]];
					end
				end
			end
		end
		ZX2 = cell(mz, nz);
		% for the inverse item
		for i = w+1 : sample_rate_inv : mx-w
			for j = w+1 : sample_rate_inv : nx-w
				for iw = -w : w
					for jw = -w : w
						izp = Zp(i, j, 1);
						jzp = Zp(i, j, 2);
						ZX2{izp + iw, jzp + jw} = [ZX2{izp + iw, jzp + jw}; [i, j] + [iw, jw]];
					end
				end
			end
		end

		% discrete optimization
		for i = w+1 : mz-w
			for j = w+1 : nz-w
				energy = 1.0e30;
				inum = size(KCHS{i, j}, 1);
                ZXPos1 = ZX1{i, j};% corresponding pixel
                ZXPos2 = ZX2{i, j};

                nZ1 = size(ZXPos1, 1);
                nZ2 = size(ZXPos2, 1);

                ZXP1 = zeros(nZ1, 3);
                for iz1 = 1 : nZ1
                	ZXP1(iz1, :) = X(ZXPos1(iz1, 1), ZXPos1(iz1, 2), :);
                end

                ZXP2 = zeros(nZ2, 3);
                for iz2 = 1 : nZ2
                	ZXP2(iz2, :) = X(ZXPos2(iz2, 1), ZXPos2(iz2, 2), :);
                end
                
                ZXn1 = size(ZXP1, 1);
				ZXn2 = size(ZXP2, 1);              
					
				for ik = 1 : inum
					ic = KCHS{i, j}(ik, 1);
					jc = KCHS{i, j}(ik, 2);

					e1 = 0;
					for izx = 1 : ZXn1
						dzxp = ZXP1(izx, :) - reshape(X(ic, jc, :), 1, 3); % difference of pixel color
						e1 = e1 + dzxp * dzxp';
					end
					e1 = e1 * weight_for;
                    
                    e2 = 0 ;
					for izx = 1 : ZXn2
						dzxp = ZXP2(izx, :) - reshape(X(ic, jc, :), 1, 3);
						e2 = e2 + dzxp * dzxp';
					end
					e2 = e2 * weight_inv;

                    e = e1 + e2;
					if e < energy
						energy = e;
						Z(i, j, :) = X(ic, jc, :);
						Xq(i, j, :) = [ic, jc];
					end
				end				
			end 
		end
		%% Inverse M-step
		% nearest neighbour data
	    ZN = zeros((mz-2*w)*(nz-2*w), nPixel);
	    for i = 1 : mz-2*w
	        for j = 1 : nz-2*w
	            idx = (i-1)*(nz-2*w) + j;
	            ZN(idx, :) = reshape(Z(i:i+2*w, j:j+2*w, :), 1, nPixel);
	        end
	    end
	    %create kdtree
	    kdt = createns(ZN,'nsmethod','kdtree');

		nXc = size(Xc, 1);
		Zp1 = zeros(nXc, 2);
		for i = 1 : nXc
			ixc = Xc(i, 1);
            jxc = Xc(i, 2);
			Xv = reshape(X(ixc-w:ixc+w, jxc-w:jxc+w, :), 1, nPixel);
			idx = knnsearch(kdt, Xv);
            [zi, zj] = Idx2Coordinate(idx, w, nz);
            Zp1(i, :) = [zi zj];
		end

		for i = w+1 : mx-w
			for j = w+1 : nx-w                
				Zp(i, j, :) = Zp1(cp(i, j), :);
			end
        end
        
		%% Forward M-step
        
        % calculate K-coherent set
        KCHS  = KCoherentSet( mz, nz, w, Xq, sq );
		for i = w+1 : mz-w
			for j = w+1 : nz-w
                energy = 1.0e30;
				% exhaustive search
                chs = KCHS{i, j};
				nc = size(chs, 1);
				for ic = 1:nc
					ci = chs(ic, 1);
                    cj = chs(ic, 2);
                    e = 0;
                    for ie = 1: c
                        e = e + norm(X(ci-w:ci+w, cj-w:cj+w, ie) - Z(i-w:i+w, j-w:j+w, ie), 'fro');
                    end
					if e < energy
						energy = e;
						Xq(i, j, :) = [ci, cj];
					end
                end
			end
        end
	end
end