function [ X ] = TextureSynthesisR(Z, m, n, w, iter_num,  X )
% Texture Synthesis

    figure;
    subplot(1, 2, 1); imshow (Z);
    Z = double(Z);
    X = double(X);
    
    % algorithm parameters
    [mz, nz, c] = size(Z);
    sample_rate = ceil(w/2);
    k = c*(2*w+1)^2; % #pixel in a window
   
    % initialize zp
    zp = zeros(m, n, 2);

    % nearest neighbour data
    ZN = zeros((mz-2*w)*(nz-2*w), k);
    for i = 1 : mz-2*w
        for j = 1 : nz-2*w
            idx = (i-1)*(nz-2*w) + j;
            ZN(idx, :) = reshape(Z(i:i+2*w, j:j+2*w, :), 1, k);
        end
    end
    %create kdtree
    kdt = createns(ZN,'nsmethod','kdtree');

    % start iteration
    for it = 1 : iter_num
        % synthesis sample matrix
        Smp(:, :, 1) = randi([w+1, m-w], m, n);
        Smp(:, :, 2) = randi([w+1, n-w], m, n);
        
        % step 1
        % calc nearest neighbour
        t1 = tic;
        zp_old = zp;
        for i = w+1 : sample_rate : m-w
            for j = w+1 : sample_rate : n-w
                is = Smp(i, j, 1);
                js = Smp(i, j, 2);
                xp = reshape(X(is-w:is+w, js-w:js+w, :), 1, k);
                idx = knnsearch(kdt, xp);
                [zi, zj] = Idx2Coordinate(idx, w, nz);
                zp(i, j, :) = [zi , zj];
            end
        end
        
        % teiminate?
        zp_old = zp_old-zp;
        if (sum(abs(zp_old(:))) < 200)
            break;
        end
        t2 = toc(t1); 
        
        % step 2
        % minimize the energy
        X_old = X;
        X = zeros(m, n, 3);
        Xf = zeros(m, n);
        Xa = zeros(m, n);
        fs = fspecial('gaussian', 2*w+1, 2*w+1);
        for i = w+1 : sample_rate : m-w
            for j = w+1 : sample_rate : n-w
                zi = zp(i, j, 1);
                zj = zp(i, j, 2);
                
                is = Smp(i, j, 1);
                js = Smp(i, j, 2);
                
                % calc weight
                dx = X_old(is-w:is+w, js-w:js+w, :) - Z(zi-w:zi+w, zj-w:zj+w, :);
                dx = dx.^2;
                weight = (1.0e-5 + sum(dx(:)))^(-0.6);
                
                % gaussian blur
                dz = Z(zi-w:zi+w, zj-w:zj+w, :);
                for h = 1:3
                    dz(:, :, h) = dz(:, :, h) .* fs;
                end
                X(is-w:is+w, js-w:js+w, :) = X(is-w:is+w, js-w:js+w, :) + dz.*weight;
                Xa(is-w:is+w, js-w:js+w) = Xa(is-w:is+w, js-w:js+w) + weight*fs;
                Xf(is-w:is+w, js-w:js+w) = ones(2*w+1);
            end
        end		
        % normalize
        for i = 1:m
            for j = 1:n
                if Xf(i, j) ~= 0
                    X(i, j, :) =  X(i, j, :)./Xa(i, j);
                else
                    X(i, j, :) = X_old(i, j, :);
                end                    
            end
        end
        
        % result
        str = num2str(t2);
        disp(['The nearest neighbors:' str]);
        str = [num2str(m) '_' num2str(n) '_' num2str(w) '_' num2str(it)];
        imwrite(uint8(X), ['.\rst\' str '.jpg']);
        disp(['iteration:' str]);
    end

    X = uint8(X);
    subplot(1, 2, 2); imshow (X);
end