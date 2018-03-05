function [B,C,P] = vsh(n,m,theta,phi)% vsh.m : vector spherical harmonics%% Usage:% [B,C,P] = vsh(n,m,theta,phi)%% Scalar n,m for the moment.% B,C,P are arrays of size length(theta,phi) x 3% theta and phi can be vectors (of equal length) or scalar.%% The three components of each vector are [r,theta,phi]%% "Out of range" n and m result in return of [0 0 0]%% This file is part of the package Optical tweezers toolbox 1.0
% Copyright 2006 The University of Queensland.
% See README.txt or README.m for license and details.
%
% http://www.physics.uq.edu.au/people/nieminen/software.html% Convert a scalar theta or phi to a vector to match a vector% partner[theta,phi] = matchsize(theta,phi);input_length = length(theta);B = zeros(input_length,3);C = zeros(input_length,3);P = zeros(input_length,3);[Y,Ytheta,Yphi] = spharm(n,m,theta,phi);B(:,2) = Ytheta;B(:,3) = Yphi;C(:,2) = Yphi;C(:,3) = - Ytheta;P(:,1) = Y;return
