X = imread('Texture-01.png');

% initial Z
mz = 32; 
nz = 32;
Z = X(1:mz, 1:nz, :);
w = 3;
[Xc, cp] = ClusterX(X, w, mz*nz);
Sq = CoherentSet(X, w, 5);

Z = InverseTextureSynthesis(X, w, Sq, Xc, cp, Z);
