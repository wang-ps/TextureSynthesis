function [ KCHS ] = KCoherentSet( mz, nz, w, Xq, sq )

    KCHS = cell(mz, nz);
    dis(:, :, 1) = diag(-w:w)*ones(2*w+1);
    dis(:, :, 2) = ones(2*w+1)*diag(-w:w);	% displacement

    for i = w+1 : mz-w
        for j = w+1 : nz-w						
            chp = Xq(i-w:i+w, j-w:j+w, :) - dis; 	% coherent pixel
            chs = [];								% coherent set
            for ic = 1 : 2*w+1
                for jc = 1 : 2*w+1   
                    chs = [chs; sq{chp(ic, jc, 1), chp(ic, jc, 2)}];
                end
            end
            KCHS{i, j} = unique(chs, 'rows');	
        end 
    end		
end

