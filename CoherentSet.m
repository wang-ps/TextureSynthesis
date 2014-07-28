function [ Sp ] = CoherentSet( X, w, level)
%COHERENTSET Summary of this function goes here
%   Detailed explanation goes here
	[m, n] = size(z);

	h = fspecial('Gaussian');

	% build Gaussian Pyramid
	% assert level >= 1
	Xg = cell(1, level);
	Xg{1} = X;
	for i = 2 : level
		X = Xg{i-1};
		X = imfilter(X, h, 'replicate');
		Xg{i} = X(1:2:end, 1:2:end, :);
	end

	for i = level : -1 : 1
		% build kdt
		
	end




end

