classdef Bsc
%Bsc abstract class representing beam shape coefficients
%
% Bsc properties:
%   a               Beam shape coefficients a vector
%   b               Beam shape coefficients b vector
%
% Bsc methods:
%   translateZ      Translates the beam along the z axis
%   translateXyz    Translation to xyz using rotations and z translations
%   translateRtp    Translation to rtp using rotations and z translations
%   farfield        Calculate fields in farfield
%   emfieldXyz      Calculate fields at specified locations
%   set.Nmax        Resize the beam shape coefficient vectors
%   get.Nmax        Get the current size of the beam shape coefficient vectors
%
% Static methods:
%   make_beam_vector    Convert output of bsc_* functions to beam coefficients
%
% Abstract methods:
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    a           % Beam shape coefficients a vector
    b           % Beam shape coefficients b vector
    type        % Coefficient type (incomming, outgoing or regular)

    dz          % Distance the beam has been translated
  end

  properties (Dependent)
    Nmax        % Current size of beam vectors
  end

  methods (Abstract)
  end

  methods (Static)
    function [a, b, n, m] = make_beam_vector(a, b, n, m, Nmax)
      %MAKE_BEAM_VECTOR converts output of bsc_* functions to sparse vectors

      if nargin < 5
        Nmax = max(n);
      end

      total_orders = ott.utils.combined_index(Nmax, Nmax);
      ci = ott.utils.combined_index(n, m);

      a = sparse(ci, 1, a, total_orders, 1);
      b = sparse(ci, 1, b, total_orders, 1);

      [n, m] = ott.utils.combined_index(1:Nmax^2+2*Nmax);
      n = n.';
      m = m.';
    end
  end

  methods (Access=protected)
    function beam = Bsc()
      %BSC construct a new beam object
      beam.dz = 0.0;
    end
  end

  methods

    function [E, H] = farfield(beam, theta, phi)
      %FARFIELD finds far field at locations thieta, phi.

      [theta,phi] = ott.utils.matchsize(theta,phi);

      [theta_new,~,indY]=unique(theta);
      [phi_new,~,indP]=unique(phi);

      E = zeros(length(theta),3);
      Etheta=zeros(length(theta),1);
      Ephi=zeros(length(theta),1);

      H = zeros(length(theta),3);
      Htheta=zeros(length(theta),1);
      Hphi=zeros(length(theta),1);

      if beam.type == 'incomming'

        a = beam.a;
        b = beam.b;
        p = zeros(size(beam.a));
        q = zeros(size(beam.b));

      elseif beam.type == 'outgoing'

        a = zeros(size(beam.a));
        b = zeros(size(beam.a));
        p = beam.a;
        q = beam.b;

      else

        % TODO: Can we convert from regular to incomming + outgoing?
        error('Unsupported beam type');

      end

      a = ott.utils.threewide(a);
      b = ott.utils.threewide(b);
      p = ott.utils.threewide(p);
      q = ott.utils.threewide(q);

      [n,m]=ott.utils.combined_index(find(abs(beam.a)|abs(beam.b)));

      for nn = 1:max(n)

        vv=find(n==nn);
        [Y,Ytheta,Yphi] = ott.utils.spharm(nn,m(vv), ...
            theta_new,zeros(size(theta_new)));

        [M,PHI]=meshgrid(m(vv),phi_new);

        expimphi=repmat(exp(1i*M.*PHI),[1,3]);

        %this makes the vectors go down in m for n.
        % has no effect if old version code.
        Nn = 1/sqrt(nn*(nn+1));

        for ii=1:length(vv)
          index=nn*(nn+1)+m(vv(ii));

          TEMP=Nn*(((1i)^(nn+1) * a(index)+(-1i)^(nn+1) ...
              * p(index)).*Yphi(:,ii)+((1i)^nn * b(index) ...
              +(-1i)^nn * q(index)).*Ytheta(:,ii));
          Etheta=Etheta+TEMP(indY).*expimphi(indP,ii);

          TEMP=Nn*(-((1i)^(nn+1) * a(index)+(-1i)^(nn+1) ...
              * p(index)).*Ytheta(:,ii)+((1i)^nn * b(index) ...
              +(-1i)^nn * q(index)).*Yphi(:,ii));
          Ephi=Ephi+TEMP(indY).*expimphi(indP,ii);

          if nargout>1
            TEMP=Nn*(((1i)^(nn+1) * b(index)+(-1i)^(nn+1) ...
                * q(index)).*Yphi(:,ii)+((1i)^nn * a(index) ...
                +(-1i)^nn * p(index)).*Ytheta(:,ii));
            Htheta=Htheta+TEMP(indY).*expimphi(indP,ii);

            TEMP=Nn*(-((1i)^(nn+1) * b(index)+(-1i)^(nn+1) ...
                * q(index)).*Ytheta(:,ii)+((1i)^nn * a(index) ...
                +(-1i)^nn * p(index)).*Yphi(:,ii));
            Hphi=Hphi+TEMP(indY).*expimphi(indP,ii);
          end
        end
      end

      E=[zeros(size(Etheta)),Etheta,Ephi];
      H=[zeros(size(Htheta)),Htheta,Hphi];

      % SI-ify units of H
      %H = H * -1i;
    end

    % TODO: emfieldXyz

    function nmax = get.Nmax(beam)
      %get.Nmax calculates Nmax from the current size of the beam coefficients
      nmax = ott.utils.combined_index(length(beam.a));
    end

    function beam = set.Nmax(beam, nmax)
      %set.Nmax resizes the beam vectors
      total_orders = ott.utils.combined_index(nmax, nmax);
      if length(beam.a) > total_orders

        amagA = full(sum(sum(abs(beam.a).^2)));
        bmagA = full(sum(sum(abs(beam.b).^2)));

        beam.a = beam.a(1:total_orders);
        beam.b = beam.b(1:total_orders);

        amagB = full(sum(sum(abs(beam.a).^2)));
        bmagB = full(sum(sum(abs(beam.b).^2)));

        aapparent_error = abs( amagA - amagB )/amagA;
        bapparent_error = abs( bmagA - bmagB )/bmagA;

        warning_error_level = 1e-6;
        if aapparent_error > warning_error_level || ...
            bapparent_error > warning_error_level
          warning('ott:Bsc:set.Nmax:truncation', ...
              ['Apparent errors of ' num2str(aapparent_error) ...
                  ', ' num2str(bapparent_error) ]);
        end
      else
        [arow_index,acol_index,aa] = find(beam.a);
        [brow_index,bcol_index,ba] = find(beam.b);
        beam.a = sparse(arow_index,acol_index,aa,total_orders,1);
        beam.b = sparse(brow_index,bcol_index,ba,total_orders,1);
      end
    end

    function beam = translateZ(beam, z)
      %TRANSLATEZ translate a beam along the z-axis

      % Add a warning when the beam is translated outside nmax2ka(Nmax) 
      % The first time may be OK, the second time does not have enough
      % information.
      % TODO: Planewave beams are not valid the first time either
      if dz > ott.utils.nmax2ka(beam.Nmax)/beam.k_medium
        warning('ott:Bsc:translateZ:outside_nmax', ...
            'Repeated translation of beam outside Nmax region');
      end
      dz = dz + abs(z);

      [A, B] = ott.utils.translate_z(beam.Nmax, z);
      beam = blkdiag(A, B) * beam;
    end

    function beam = translateXyz(beam, x, y, z)
      %TRANSLATEXYZ translate the beam given Cartesian coordinates
      if nargin == 2
        rtp = xyz2rtp(x);
      else
        rtp = xyz2rtp(x, y, z);
      end
      beam = beam.translateRtp(rtp);
    end

    function beam = translateRtp(beam, r, theta, phi)
      %TRANSLATERTP translate the beam given spherical coordinates
      if nargin == 2
        [r, theta, phi] = r(1:3);
      end
      [~, D] = beam.rotateYz(theta, phi);
      beam = D * beam;
      beam.translateZ(r);
      beam = D' * beam;
    end

    function [beam, D] = rotate(beam, R)
      %ROTATE apply the rotation matrix R to the beam coefficients
      D = ott.utils.wigner_rotation_matrix(beam.Nmax, R);
      beam = D * beam;
    end

    function [beam, D] = rotateX(beam, angle)
      %ROTATEX rotates the beam about the x-axis an angle in radians
      [beam, D] = beam.rotate(rotx(angle));
    end

    function [beam, D] = rotateY(beam, angle)
      %ROTATEX rotates the beam about the y-axis an angle in radians
      [beam, D] = beam.rotate(roty(angle));
    end

    function [beam, D] = rotateZ(beam, angle)
      %ROTATEX rotates the beam about the z-axis an angle in radians
      [beam, D] = beam.rotate(rotz(angle));
    end

    function [beam, D] = rotateXy(beam, anglex, angley)
      %ROTATEX rotates the beam about the x then y axes
      [beam, D] = beam.rotate(roty(angley)*rotx(anglex));
    end

    function [beam, D] = rotateXz(beam, anglex, anglez)
      %ROTATEX rotates the beam about the x then z axes
      [beam, D] = beam.rotate(rotz(anglez)*rotx(anglex));
    end

    function [beam, D] = rotateYz(beam, angley, anglez)
      %ROTATEX rotates the beam about the y then z axes
      [beam, D] = beam.rotate(rotz(anglez)*roty(angley));

    function [beam, D] = rotateXyz(beam, anglex, angley, anglez)
      %ROTATEX rotates the beam about the x, y then z axes
      [beam, D] = beam.rotate(rotz(anglez)*roty(angley)*rotx(anglex));
    end
  end
end
