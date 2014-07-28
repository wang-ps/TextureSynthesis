function [ Z ] = TextureSynthesis(X, w, iter_num, sq, Xc, cp, Xp)
% Inverse Texture Synthesis

    
	[mx, nx, c] = size(X);
	[mz, nz, c] = size(Xp);

	alpha = 0.01;
	sample_rate_inv = floor(w/2);
	sample_rate_for = floor(w/4);
	weight_inv = (mx*nx) / sample_rate_inv^2 ;
	weight_for = alpha * (mz*nz) / sample_rate_for^2;

	nPixel = c*(2*w+1)^2; % #pixel in a window
	Zp = zeros(mx, nx, 2); % Xp -> Zp, inverse item
	%Xq = zeros(mz, nz, 2); % Zq -> Xq, forward item
	X = zeros(mz, nz, c);


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

    
	for it = 1:iter_num
	
		% ---------
		%% z E-step
		% calc k-coherence set
		KCHS = cell(mz, nz);
		dis(:, :, 1) = diag(-w:w)*ones(2*w+1);
		dis(:, :, 2) = ones(2*w+1)*diag(-w:w);	% displacement

		for i = w+1 : mz-w
			for j = w+1 : nz-w			
				
				chp = Xq(i-w:i+w, j-w:j+w, :) - dis; 	% coherent pixel
				chs = [];								% coherent set
				for ic = 1 : 2*w+1
					for jc = 1 : 2*w+1
						chs = [chs; sq(chp(ic, jc))];
					end
				end
				KCHS{i, j} = unique(chs, 'rows');			

			end 
		end		

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
			for j = w+1 : nx-w
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
				energy = 1.0*e30;
				inum = size(KCHS{i, j}, 1);
				ZXP1 = X(ZX1{i, j}, :); % corresponding pixel color
				ZXP2 = X(ZX2{i, j}, :);
				ZXn1 = size(ZX1, 1);
				ZXn2 = size(ZX2, 1);
					
				for ik = 1 : inum
					ic = KCHS{i, j}(ik, 1);
					jc = KCHS{i, j}(ik, 2);

					e = 0;
					for izx = 1 : ZXn1
						dzxp = ZXP1(izx, :) - X(ic, jc, :); % difference of pixel color
						e = e + dzxp' * dzxp;
					end
					e = e * weight_for;

					for izx = 1 : ZXn2
						dzxp = ZXP2(izx, :) - X(ic, jc, :);
						e = e + dzxp' * dzxp;
					end
					e = e * weight_inv;

					if e < energy
						energy = e;
						Z(i, j, :) = X(ic, jc, :);
						Xq(i, j, :) = [ic, jc];
					end
				end				
			end 
		end

		% ----------
		%% Inverse M-step
		% nearest neighbour data
	    ZN = zeros((mz-2*w)*(nz-2*w), k);
	    for i = 1 : mz-2*w
	        for j = 1 : nz-2*w
	            idx = (i-1)*(nz-2*w) + j;
	            ZN(idx, :) = reshape(Z(i:i+2*w, j:j+2*w, :), 1, nPixel);
	        end
	    end
	    %create kdtree
	    kdt = createns(ZN,'nsmethod','kdtree');

		nXc = size(Xc, 1);
		Zp1 = size(nXc, 2);
		for i = 1 : nXc
			ixc = Xc(i, 1);
            jxc = Xc(i, 2);
			Xv = reshape(X(ixc-w:ixc+w, jxc-w:jxc+w, :), 1, nPixel);
			idx = knnsearch(kdt, Xv);
			zj = mod(idx, nz-2*w);
			if zj == 0
                zj = w;
            end
            zi = fix((idx-zj)/(nz-2*w) + 1);
            Zp1(i, :) = [zi+w zj+w];
		end

		for i = 1 : mx-2*w
			for j = 1 : nx-2*w
				Zp(i, j, :) = Zp1(cp(i, j), :);
			end
		end

		% ----------
		%% Forward M-step
		for i = w+1 : zm-w
			for j = w+1 : zn-w
				% calc k-coherence set
				dis(:, :, 1) = diag(-w:w)*ones(2*w+1);
				dis(:, :, 2) = ones(2*w+1)*diag(-w:w);	% displacement
				chp = Xp(i-w:i+w, j-w:j+w, :) - dis; 	% coherent pixel
				chs = [];					% coherent set
				for ic = 1 : 2*w+1
					for jc = 1 : 2*w+1
						chs = [chs; sq(chp(i, j))];
					end
				end
				chs = unique(chs, 'rows');

				% exhaustive search
				nc = size(chs, 1);
				for ic = 1:nc
					ci = chs(ic, 1);
                    cj = chs(ic, 2);
					e = norm(X(ci-w:ci+w, cj-w:cj+w, :) - Z(i-w:i+w, j-w:j+w, :), 'fro');
					if e < energy
						energy = e;
						Xp(i, j, :) = [ci, cj];
					end
				end

			end
		end
	end
end