  %% input sample texture 
  %   Z = imread('.\rst\Texture-01.png');
  Z = imread('rst.jpg');
  
  % build pyramid
  Zp{1} = Z;
  Z2 = imfilter(Z, fspecial('gaussian'), 'replicate');
  Zp{2} = Z2(1:2:end, 1:2:end, :);
  figure;
  subplot(1, 2, 1); imshow(Zp{1}); 
  subplot(1, 2, 2); imshow(Zp{2});
  
  %% level 2
  mx = 64;
  nx = 64;
  w = 3;
  
  disp('-----------');
  disp('level 2, low resolution');
  
  % initialize X for level 2
  %   X = imresize(Z(1:50, 1:50, :), [100, 100]);  
  [nz, mz, ~] = size(Zp{2});
  X = zeros(mx/2, nx/2, 3);
  for i = w+1 : w : mx/2-w
      for j = w+1 : w: nx/2-w
          ii = randi([w+1, mz-w], 1);
          jj = randi([w+1, nz-w], 1);
          X(i-w:i+w, j-w:j+w, :) = Zp{2}(ii-w:ii+w, jj-w:jj+w, :);
      end
  end
  
  X = TextureSynthesis(Zp{2}, mx/2, nx/2, 4, 3, X);
  X = TextureSynthesis(Zp{2}, mx/2, nx/2, 2, 3, X);

  %% level 1
  disp('-----------');
  disp('level 1, high resolution');

  % initialize X for level 2
  X1 = zeros(mx, nx, 3);
  X1(1:2:end, 1:2:end, :) = X;
  X = 4*imfilter(X1, fspecial('gaussian'), 'replicate');

%   X = TextureSynthesis(Z, mx, nx, 8, 3, X);
  X = TextureSynthesis(Z, mx, nx, 4, 3, X);
  X = TextureSynthesis(Z, mx, nx, 2, 5, X);
  X = TextureSynthesis(Z, mx, nx, 2, 5, X);
  