% Example file for finding spring constants
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

% Make warnings less obtrusive
ott_warning('once');
change_warnings('off');

% Specify refractive indices
n_medium = 1.34;
n_particle = 1.59;
n_relative = n_particle/n_medium;

% If you want to give all measurements in wavelengths in the surrounding
% medium, then:
wavelength = 1;
% wavelength = wavelength0 / n_medium;
% else you can give it in any units you want. Only k times lengths matters
k = 2*pi/wavelength;

radius = 1.5;
Nmax = ka2nmax(k*radius);

if Nmax < 12
    Nmax = 12;
end

% a Gaussian beam: w0 = 2/(k*tan(theta))
beam_angle = 50; % Convergence half-angle of 50 degrees

% Polarisation. [ 1 0 ] is plane-polarised along the x-axis, [ 0 1 ] is
% y-polarised, and [ 1 -i ] and [ 1 i ] are circularly polarised.
polarisation = [ 1 0 ];

[a,b] = bsc_pointmatch_farfield(Nmax,1,[ 0 0 beam_angle 1 polarisation 90 ]);

% If you're going to do a range of particles, then the T-matrix has to
% calculated inside the loop.

% To search for a refractive index, I recommend a bisection search. For an
% example of bisection search, see find_axial_equilibrium.m

T = tmatrix_mie(Nmax,k,k*n_relative,radius);

[z,k] = axial_equilibrium(T,a,b)

