X = imread('Texture-01.png');

% initial Z
mz = 16; 
nz = 16;
Z = X(1:16, 1:16, :);

[Xc, cp] = ClusterX(X, 1, mz*nz);
Sq = CoherentSet(X, 1, 5);

Z = InverseTextureSynthesis(X, 1, Sq, Xc, cp, Z);
