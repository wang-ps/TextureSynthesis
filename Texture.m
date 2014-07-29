%% input sample texture 
Z = imread('.\rst\Texture-01.png');
  
  % build pyramid
  Zp{1} = Z;
  Z2 = imfilter(Z, fspecial('gaussian'), 'replicate');
  Zp{2} = Z2(1:2:end, 1:2:end, :);
  figure;
  subplot(1, 2, 1); imshow(Zp{1}); 
  subplot(1, 2, 2); imshow(Zp{2});
  %% level 2
  disp('-----------');
  disp('level 2, low resolution');

  % initialize X for level 2
%   X = imresize(Z(1:50, 1:50, :), [100, 100]);
  w = 12;
  [nz, mz, ~] = size(Zp{2});
  X = zeros(100, 100, 3);
  for i = w+1 : w : 100-w
      for j = w+1 : w: 100-w
          ii = randi([w+1, mz-w], 1);
          jj = randi([w+1, nz-w], 1);
          X(i-w:i+w, j-w:j+w, :) = Zp{2}(ii-w:ii+w, jj-w:jj+w, :);
      end
  end
  
  X = TextureSynthesis(Zp{2}, 100, 100, 8, 3, X);
  X = TextureSynthesis(Zp{2}, 100, 100, 4, 3, X);

  %% level 1
  disp('-----------');
  disp('level 1, high resolution');

  % initialize X for level 2
  X1 = zeros(200, 200, 3);
  X1(1:2:end, 1:2:end, :) = X;
  X = 4*imfilter(X1, fspecial('gaussian'), 'replicate');

  X = TextureSynthesis(Z, 200, 200, 8, 3, X);
  X = TextureSynthesis(Z, 200, 200, 4, 3, X);
  X = TextureSynthesis(Z, 200, 200, 2, 5, X);
  