  Z = imread('texture.png');
  
  % level 1
  disp('-----------');
  disp('level 1');
  X = imresize(Z(1:50, 1:50, :), [100, 100]);
  X = TextureSynthesis(Z, 100, 100, 12, X);
  X = TextureSynthesis(Z, 100, 100, 8, X);
  X = TextureSynthesis(Z, 100, 100, 4, X);

  % level 2
  disp('-----------');
  disp('level 1');
  X = imresize(X, [200, 200]);
  X = TextureSynthesis(Z, 200, 200, 12, X);
  X = TextureSynthesis(Z, 200, 200, 8, X);
  X = TextureSynthesis(Z, 200, 200, 4, X);
  