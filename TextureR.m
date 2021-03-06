  %% input sample texture 
    Z = imread('.\Image\Texture-01.png');
   % Z = imread('.\rst\InvTexture-02_rst.jpg');
  
  % build pyramid
  Zp{1} = Z;
  Z2 = imfilter(Z, fspecial('gaussian'), 'replicate');
  Zp{2} = Z2(1:2:end, 1:2:end, :);
  figure;
  subplot(1, 2, 1); imshow(Zp{1}); 
  subplot(1, 2, 2); imshow(Zp{2});
  
  %% level 2 ( low resolution )
  mx = 200;
  nx = 200;
  w = 5;
  
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
  
  % start level 2
  X = TextureSynthesisR(Zp{2}, mx/2, nx/2, 10, 5, X);
  X = TextureSynthesisR(Zp{2}, mx/2, nx/2, 5, 5, X);

  %% level 1
  disp('-----------');
  disp('level 1, high resolution');

  % initialize X for level 2
  X1 = zeros(mx, nx, 3);
  X1(1:2:end, 1:2:end, :) = X;
  X = 4*imfilter(X1, fspecial('gaussian'), 'replicate');

%   X = TextureSynthesisR(Z, mx, nx, 8, 3, X);
  X = TextureSynthesisR(Z, mx, nx, 12, 5, X);
  X = TextureSynthesisR(Z, mx, nx, 8, 5, X);
  X = TextureSynthesisR(Z, mx, nx, 4, 10, X);
  X = TextureSynthesisR(Z, mx, nx, 2, 10, X);