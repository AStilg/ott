classdef TmatrixMie < ott.Tmatrix
%TmatrixMie construct T-matrix from Mie scattering coefficients
%
% TmatrixMie properties:
%   radius            The radius of the sphere the T-matrix represents
%   k_medium          Wavenumber in the trapping medium
%   k_particle        Wavenumber of the particle
%
% This class is based on tmatrix_mie.m and tmatrix_mie_layered.m from ottv1.
%
% This file is part of the optical tweezers toolbox.
% See LICENSE.md for information about using/distributing this file.

  properties (SetAccess=protected)
    radius            % Radius of particle
    k_medium          % Wavenumber of medium
    k_particle        % Wavenumber of particle
  end

  methods (Access=protected)
    function T = tmatrix_mie(tmatrix, Nmax, internal)
      %TMATRIX_MIE code from tmatrix_mie.m

      n=[1:Nmax];

      m = tmatrix.k_particle/tmatrix.k_medium;

      r0 = tmatrix.k_medium * tmatrix.radius;
      r1 = tmatrix.k_particle * tmatrix.radius;

      indexing=ott.utils.combined_index(1:Nmax^2+2*Nmax)';

      import ott.utils.sbesselj
      import ott.utils.sbesselh1

      j0 = (sbesselj(n,r0)).';
      j1 = (sbesselj(n,r1)).';
      h0 = (sbesselh1(n,r0)).';
      j0d = (sbesselj(n-1,r0) - n.*sbesselj(n,r0)/r0).';
      j1d = (sbesselj(n-1,r1) - n.*sbesselj(n,r1)/r1).';
      h0d = (sbesselh1(n-1,r0) - n.*sbesselh1(n,r0)/r0).';

      if internal == false
        % Calculate external T-matrix
        b = -( j1d.*j0 - m*j0d.*j1 ) ./ ( j1d.*h0 - m*h0d.*j1 );
        a = -( j0d.*j1 - m*j1d.*j0 ) ./ ( h0d.*j1 - m*j1d.*h0 );
        T=sparse([1:2*(Nmax^2+2*Nmax)],[1:2*(Nmax^2+2*Nmax)], ...
            [a(indexing);b(indexing)]);
      else
        % Calculate internal T-matrix
        d = -( h0d.*j0 - j0d.*h0 ) ./ ( m*j1.*h0d - j1d.*h0 );
        c = -( j0d.*h0 - h0d.*j0 ) ./ ( m*j1d.*h0 - h0d.*j1 );
        T=sparse([1:2*(Nmax^2+2*Nmax)],[1:2*(Nmax^2+2*Nmax)], ...
            [c(indexing);d(indexing)]);
      end
    end

    function T = tmatrix_mie_layered(tmatrix, Nmax, internal)
      %TMATRIX_MIE code from tmatrix_mie_layered.m

      n_layer=[tmatrix.k_particle,tmatrix.k_medium]/2/pi;

      n=[1:Nmax]; %n for nmax

      [ha_0,hb_0]=tmatrix.layerrecursion(n,n_layer, ...
          tmatrix.radius,1,1/n_layer(1),1);

      m_0 = n_layer(end-1);
      m_1 = n_layer(end);

      r1 = 2*pi*m_1 * tmatrix.radius(end);

      import ott.utils.sbesselj
      import ott.utils.sbesselh1

      j1 = (sbesselj(n,r1)).';
      h1 = (sbesselh1(n,r1)).';
      j1d = (sbesselj(n-1,r1) - n.*sbesselj(n,r1)/r1).';
      h1d = (sbesselh1(n-1,r1) - n.*sbesselh1(n,r1)/r1).';

      d1_1=j1d./j1;
      d3_1=h1d./h1;
      r_1=j1./h1;

      % %TM/TE coeffs...
      al_1=r_1.*(m_1*ha_0-m_0*d1_1)./(m_1*ha_0-m_0*d3_1);
      bl_1=r_1.*(m_0*hb_0-m_1*d1_1)./(m_0*hb_0-m_1*d3_1);
      % aL_1=((ha_0/m_0+n.'/r1).*j1-sbesselj(n-1,r1).') ...
      %   ./((ha_0/m_0+n.'/r1).*h1-sbesselh1(n-1,r1).');
      % bL_1=((m_0*hb_0+n.'/r1).*j1-sbesselj(n-1,r1).') ...
      %   ./((m_0*hb_0+n.'/r1).*h1-sbesselh1(n-1,r1).');

      %swap the modes to TE/TM... and negatize.
      b = -al_1;
      a = -bl_1;

      %t-matrix indices...
      indexing=ott.utils.combined_index(1:Nmax^2+2*Nmax)';

      if internal == false
        T=sparse([1:2*(Nmax^2+2*Nmax)],[1:2*(Nmax^2+2*Nmax)], ...
            [a(indexing);b(indexing)]);
      else
        error('This release does not calculate internal coefficients.')
        %because I have to use different ratio functions to find c/d

        % TODO: Calculate internal T-matrix for layered spheres
      end
    end
  end

  methods (Static, Access=private)
    function [ha_0, hb_0]=layerrecursion(n,n_layer, ...
        radius,recursionnumber,ha_n,hb_n)
      %LAYERRECURSION function used by tmatrix_mie_layered.
      %
      % Internal function which should never be bare. This function has a lot
      % of overhead but takes very little time to complete so I'm happy with
      % it as is. I made the recursion start from 1 because it should...

      import ott.utils.sbesselj
      import ott.utils.sbesselh1

      m_0 = n_layer(recursionnumber);
      r0 = 2*pi*m_0 * radius(recursionnumber);
      j0 = (sbesselj(n,r0)).';
      h0 = (sbesselh1(n,r0)).';
      j0d = (sbesselj(n-1,r0) - n.*sbesselj(n,r0)/r0).';
      h0d = (sbesselh1(n-1,r0) - n.*sbesselh1(n,r0)/r0).';

      d1_0=j0d./j0;
      d3_0=h0d./h0;
      r_0=j0./h0;

      if recursionnumber>1
        m_n = n_layer(recursionnumber-1);
        rn = 2*pi*m_0 * radius(recursionnumber-1);
        jn = (sbesselj(n,rn)).';
        hn = (sbesselh1(n,rn)).';
        jnd = (sbesselj(n-1,rn) - n.*sbesselj(n,rn)/rn).';
        hnd = (sbesselh1(n-1,rn) - n.*sbesselh1(n,rn)/rn).';

        d1_n=jnd./jn;
        d3_n=hnd./hn;
        r_n=jn./hn;
      else
        m_n=1;
        rn=0;
        jn=0;
        hn=0;
        jnd=0;
        hnd=0;

        d1_n=0;
        d3_n=0;
        r_n=0;
      end

      g0=m_0*ha_n-m_n*d1_n;
      g1=m_0*ha_n-m_n*d3_n;
      g0b=m_n*hb_n-m_0*d1_n;
      g1b=m_n*hb_n-m_0*d3_n;
      q_0=r_n./r_0;

      ha_0=(g1.*d1_0-q_0.*g0.*d3_0)./(g1-q_0.*g0);
      hb_0=(g1b.*d1_0-q_0.*g0b.*d3_0)./(g1b-q_0.*g0b);

      %doesn't calculate last layer with this recursion...
      if recursionnumber<length(radius)
        recursionnumber=recursionnumber+1;
        [ha_0,hb_0]=layerrecursion(n,n_layer,radius,recursionnumber,ha_0,hb_0);
      end
    end
  end

  methods
    function tmatrix = TmatrixMie(radius, varargin)
      %TMATRIXMIE construct a new Mie T-matrix for a sphere with size radius.
      %
      % If radius and k_particle (see bellow) are vectors, calculates the
      % coefficients for a leyered sphere with radius and k_particle
      % specifying the radius and wavenumber of each layer starting
      % from the core and going out.
      %
      % Note: Will not compute for particles of 500 layers or more. Not all of
      % the recursions in the reference are implemented because the original
      % author was lazy. For layered sphere, Nmax of 100 is numerically stable.
      %
      %  TMATRIXMIE(..., 'Nmax', Nmax) specifies the size of the
      %  T-matrix to use.  If not specified, the size is calculated
      %  from ott.utils.ka2nmax(radius*k_medium) or set to 100 for
      %  layered spheres.
      %
      %  TMATRIXMIE(..., 'k_medium', k)
      %  or TMATRIXMIE(..., 'wavelength_medium', wavelength)
      %  or TMATRIXMIE(..., 'index_medium', index)
      %  specify the wavenumber, wavelength or index in the medium.
      %
      %  TMATRIXMIE(..., 'k_particle', k)
      %  or TMATRIXMIE(..., 'wavelength_particle', wavelength)
      %  or TMATRIXMIE(..., 'index_particle', index)
      %  specify the wavenumber, wavelength or index in the particle.
      %
      %  TMATRIXMIE(..., 'wavelength0', wavelength) specifies the
      %  wavelength in the vecuum, required when index_particle or
      %  index_medium are specified.
      %
      %  TMATRIXMIE(..., 'internal', internal) if true, calculates the
      %  T-matrix for the internal coefficients.
      %
      % Layered sphere reference:
      % "Improved recursive algorithm for light scattering by a multilayered
      % sphere", Wen Yang, Applied Optics 42(9), 2003

      tmatrix = tmatrix@ott.Tmatrix();

      % Parse inputs
      p = inputParser;
      p.addParameter('Nmax', []);
      p.addParameter('k_medium', []);
      p.addParameter('wavelength_medium', []);
      p.addParameter('index_medium', []);
      p.addParameter('k_particle', []);
      p.addParameter('wavelength_particle', []);
      p.addParameter('index_particle', []);
      p.addParameter('wavelength0', []);
      p.addParameter('internal', false);
      p.parse(varargin{:});

      % Store inputs: radius, k_medium, k_particle
      tmatrix.radius = radius;
      tmatrix.k_medium = tmatrix.parser_k_medium(p);
      tmatrix.k_particle = tmatrix.parser_k_particle(p);

      % Check radius and k_particle are similar lengths
      if numel(tmatrix.radius) ~= numel(tmatrix.k_particle)
        error('radius and k_particle must be the same length');
      end

      % Check number of layers
      if numel(tmatrix.radius) >= 100
        warning('May not work well for particles with >= 100 layers');
      end

      % If Nmax not specified, choose a good Nmax
      if isempty(p.Results.Nmax)
        if numel(tmatrix.radius) == 1
          Nmax = ott.utils.ka2nmax(tmatrix.k_medium*radius);
        else
          Nmax = 100;
        end
      else
        Nmax = p.Results.Nmax;
      end

      % TODO: We should store if it is internal or external
      % TODO: Is this a total or scattered T-matrix?  Store this too.
      % TODO: For layered spheres, can we resize Nmax after calculation?

      % Calculate the T-matrix and store it
      if numel(tmatrix.radius) == 1
        tmatrix.data = tmatrix.tmatrix_mie(Nmax, p.Results.internal);
      else
        tmatrix.data = tmatrix.tmatrix_mie_layered(Nmax, p.Results.internal);
      end
    end
  end
end
