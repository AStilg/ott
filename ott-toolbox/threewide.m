function wide_vector = threewide( a )
% threewide.m - converts an input vector (either row or column
%               vector) into a column vector repeated in
%               three columns.
% Usage:
% wide_vector = threewide(original_vector);
%
% You might find this useful for multiplying a vector of scalars
% with a column vector of 3-vectors.
%
% PACKAGE INFO

a = a(:);
wide_vector = [ a a a ];

return

