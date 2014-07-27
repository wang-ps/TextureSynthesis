function [ Z ] = TextureSynthesis(X, w, iter_num, KC, Xc, Z0)
% Inverse Texture Synthesis

	[mx, nx, c] = size(X);
	[mz, nz, c] = size(Z0);

	Zq = zeros(mx, nx, 2);
	Xp = zeros(mz, nz, 2);

	for it = 1:iter_num
		% z E-step
		for i = w+1 : zm-w
			for j = w+1 : zn-w
				
				% calc k-coherence set
				dis(:, :, 1) = diag(-w:w)*ones(2*w+1);
				dis(:, :, 2) = ones(2*w+1)*diag(-w:w);
				chp = Xp(i-w:i+w, j-w:j+w, :) - dis; % coherent pixel
				chs = KC(chp(1, 1));				 % coherent set
				for ic = 1 : 2*w+1
					for jc = 1 : 2*w+1
						chs = union(chs, KC(chp(i, j)));
					end
				end

				% exhaustive search
				energy = 1.0e30;
				[nc] = size(chs, 1);
				for ic = 1:nc
					[ci cj] = chs(ic);
					e = norm(X(ci-w:ci+w, cj-w:cj+w, :) - Z(i-w:i+w, j-w:j+w, :), 'fro');
					if e < energy
						energy = e;
						Xp(i, j, :) = [ci, cj];
						Z(i, j, :) = X(ci, cj, :);
					end
				end

			end 
		end

		  
	end



end