filename = 'InvTexture-02';
X = imread([ '.\Image\', filename, 'png']);
X = X(50:250, 50:250, :);

% initial Z
mz = 50; 
nz = 50;
Z = X(1:mz, 1:nz, :);
w = 3;
[Xc, cp] = ClusterX(X, w, mz*nz);
Sq = CoherentSet(X, w, 5);

Z = InverseTextureSynthesis(X, w, Sq, Xc, cp, Z);
imwrite(Z(w+1:mz-w, w+1:nz-w, :), ['.\rst\' filename '_rst.jpg']);
