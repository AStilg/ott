function [force,torque] = force_z(ibeam, sbeam)
% force_z.m: Finds z component of optical force and torque
% TODO: Clean up documentation
%
% Usage:
% [Fz,Tz] = force_z(n,m,a,b,p,q)
% OR
% [Fz,Tz] = force_z(n,m,ab,pq)
%
% What units are you using for a,b,p,q?
% If you have simple units like (using incoming/outgoing):
%     power = sum( abs(a).^2 ... )
% then divide by c to get newtons, divide by omega to get N.m
% If you have
%    sum( abs(a).^2 + abs(b).^2 ) = 1
% then the force and torque are in units of the momentum per photon
% and hbar per photon.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

warning('ott:force_z:depreciated', ...
    'This file will be removed in the next release');

% Ensure beams are the same size
if ibeam.Nmax > sbeam.Nmax
  sbeam.Nmax = ibeam.Nmax;
elseif ibeam.Nmax < sbeam.Nmax
  ibeam.Nmax = sbeam.Nmax;
end

% Ensure the beam is incomming-outgoing
sbeam = sbeam.toOutgoing(ibeam);

% Get the relevent beam coefficients
[a, b] = ibeam.getCoefficients();
[p, q] = sbeam.getCoefficients();
[n, m] = ibeam.getModeIndices();

force = forcez(n,m,a,b) - forcez(n,m,p,q);
torque = sum( m.*( abs(a).^2 + abs(b).^2 - abs(p).^2 - abs(q).^2 ) );

return


% Find z-component of force
% Magic formula from Crichton
function fz = forcez(n,m,a,b)

import ott.*
import ott.utils.*

ci = combined_index(n,m);

aa = zeros(max(ci),1);
bb = zeros(max(ci),1);
aap = zeros(max(ci),1);
bbp = zeros(max(ci),1);

[nn,mm] = combined_index((1:max(ci))');

aa(ci) = a;
bb(ci) = 1i*b;

n1 = find( n>1 & n>abs(m));

ci1 = combined_index(n(n1)-1,m(n1));

aap(ci1) = a(n1);
bbp(ci1) = 1i*b(n1);

fz = 2 * mm ./ nn ./ (nn+1) .* imag( conj(aa) .* bb ) ...
    - 2 ./ (nn+1) .* sqrt( nn .* (nn+2) .*  (nn-mm+1) .* (nn+mm+1) ./ (2*nn+1) ./ (2*nn+3) ) ...
    .* imag( aa.*conj(aap) + bb.*conj(bbp) );

fz = sum(fz);

return

