function [force,torque,spin] = forcetorque(ibeam, sbeam)
% forcetorque.m
%
% Finds optical force and torque
% TODO: Clean up documentation
%
% Usage:
% force = forcetorque(n,m,a,b,p,q)
% or
% [force,torque] = forcetorque(n,m,a,b,p,q)
% or
% [force,torque,spin] = forcetorque(n,m,a,b,p,q)
% or
% [force,torque,spin] = forcetorque(n,m,ab,pq)
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

warning('ott:forcetorque:depreciated', ...
    ['This file will be replaced with ' ...
    'force_torque_farsund.m in the next release']);

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

force = [ 0 0 0 ];
torque = [ 0 0 0 ];
spin = [ 0 0 0 ];

% z-component is easiest
force(3) = forcez(n,m,a,b) - forcez(n,m,p,q);
torque(3) = sum( m.*( abs(a).^2 + abs(b).^2 - abs(p).^2 - abs(q).^2 ) );
spin(3) = spinz(n,m,a,b) - spinz(n,m,p,q);

% Now find x,y components by rotating by 90 degrees and re-using the
% z-component formulae

% First, rotate x axis onto z axis

[~, D] = ibeam.rotateY(pi/2);

a2 = D*a;
b2 = D*b;
p2 = D*p;
q2 = D*q;

force(1) = forcez(n,m,a2,b2) - forcez(n,m,p2,q2);
torque(1) = sum( m.*( abs(a2).^2 + abs(b2).^2 - abs(p2).^2 - abs(q2).^2 ) );
spin(1) = spinz(n,m,a2,b2) - spinz(n,m,p2,q2);

% Finally, rotate (original) y axis onto z axis

[~, D] = ibeam.rotateX(pi/2);

a2 = D*a;
b2 = D*b;
p2 = D*p;
q2 = D*q;

force(2) = forcez(n,m,a2,b2) - forcez(n,m,p2,q2);
torque(2) = sum( m.*( abs(a2).^2 + abs(b2).^2 - abs(p2).^2 - abs(q2).^2 ) );
spin(2) = spinz(n,m,a2,b2) - spinz(n,m,p2,q2);

end


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

% Also magic formula for z-component of spin
function sz = spinz(n,m,a,b)

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

sz = mm ./ nn ./ (nn+1) .* ( abs(aa).^2 + abs(bb).^2 ) ...
    - 2 ./ (nn+1) .* sqrt( nn .* (nn+2) .*  (nn-mm+1) .* (nn+mm+1) ./ (2*nn+1) ./ (2*nn+3) ) ...
    .* real( aa.*conj(bbp) - bb.*conj(aap) );

sz = sum(sz);

return

