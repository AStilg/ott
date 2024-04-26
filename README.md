ott - Optical Tweezers Toolbox
==============================

[![DOI](https://zenodo.org/badge/123386773.svg)](https://zenodo.org/badge/latestdoi/123386773)
[![Documentation Status](https://readthedocs.org/projects/ott/badge/?version=latest)](https://ott.readthedocs.io/en/latest/?badge=latest)
[![View Optical Tweezers Toolbox on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://au.mathworks.com/matlabcentral/fileexchange/73541-optical-tweezers-toolbox)

The optical tweezers toolbox can be used to calculate optical forces
and torques of particles using the T-matrix formalism in a vector
spherical wave basis.
The toolbox includes codes for calculating T-matrices, beams described
by vector spherical wave functions, functions for calculating forces
and torques, simple codes for simulating dynamics and examples.

We are currently working on [documentation](https://ott.readthedocs.io/)
and we welcome feedback/suggestions/comments.  Additional documentation
can be found via the Matlab `help` command or in the source code.

Installation and usage
----------------------

There are several methods for installing the toolbox.
If using Matlab, the easiest method is to launch Matlab and
navigate to Home -> Addons -> Get-Addons and search for
"Optical Tweezers Toolbox".  Then, simply click the
"Add from GitHub" button to automatically download the package and
add it to the path.
Alternatively, you can download the toolbox directly from the
[GitHub repository](https://github.com/ilent2/ott) or select a
specific release; if using this method you will need to add the
toolbox directory to Matlab's path using

```matlab
addpath('<download-path>/ott');
help ott   % Test that ott was found, should display ott Contents page
```

Regardless of the method you acquired OTT with, once installed you
should be able to access the toolbox functions and classes contained
in the `ott.*` package, for example

```matlab
beam = ott.BscPmGauss();
```

or for the graphical user interface

```matlab
ott.ui.Launcher
```

More detailed instructions can be found in the Getting Started
section of the [documentation](https://ott.readthedocs.io/).

Dependencies
------------

The toolbox runs on recent versions of Matlab, most functionality
should work on at least R2016b but the graphical user interface might
need R2018a or newer.
Most functionality should work with
[GNU Octave](https://www.gnu.org/software/octave/), however this
has not been tested recently and performance is optimised for Matlab.

Some functionality may require additional dependences including
additional Matlab products.
We are currently working on a full list; however, if you encounter any
difficulties with missing dependencies, please let us know and we may
be able to find a workaround.

Quick-start guide
-----------------

The toolbox has changed a lot since previous releases, the most notable
change is addition of a graphical user interface (still a work-in-progress)
and moving from a folder structure to a Matlab package structure.
To get started, take a look at
the [documentation over on read the docs](https://ott.readthedocs.io).
The documentation source can be found in the [docs](docs) directory or
you can download a PDF copy of the documentation with the latest release.
Alternatively, take a look at the [examples directory](examples).

To calculate the force on a spherical particle, you need to setup a beam
object, setup a particle and run the `ott.forcetorque` function.
For example:

```matlab
beam = ott.BscPmGauss();
tmatrix = ott.TmatrixMie(0.5, 'k_medium', 2*pi, 'k_particle', 2*pi*1.3);
fz = ott.forcetorque(beam, tmatrix, 'position', [0;0;1].*linspace(-8, 8));
figure(), plot(fz.')
```

This example assumes everything is in units of the wavelength,
and creates a Gaussian beam with the default parameters.

Recent changes
--------------

There have been many changes since the previous release, mainly the switch
to object orientated programming.  Beams and T-matrices are now represented
by objects, we have added shape objects and moved everything into packages.

T-matrices are represented by `Tmatrix` objects.  For simple shapes,
the `Tmatrix.simple` method can be used to construct T-matrices for
a variety of common objects.
More complex T-matrices can be generated by inheriting the T-matrix
class, for an example, take a look at TmatrixMie and TmatrixPm.

Beams are represented by a `Bsc` objects.  A beam can be multiplied
by T-matrices or other matrix/scalar values to generate new beams.
For Gaussian type beams, including Hermite-Gauss, Ince-Gauss, and
Laguarre-Gaussian beams, the `BscPmGauss` class provides the
equivalent of `bsc_pointmatch_farfield` in the previous release.

The new implementation hides `Nmax`, most routines have a default
choice of `Nmax` based on the beam/particle size.  `Nmax` can still
be accessed and changed manually, but in most cases the automatic
choice of `Nmax` should be fine.
Beams can T-matrices can be multiplied without needing to
worry about the having equal `Nmax`, the beam/T-matrix will be
expanded to match the maximum `Nmax`.
If repeated calculations are being done, it may be faster to first
ensure the `Nmax` of the beam and T-matrix match, this is done in
`forcetorque` when the position or rotation arguments are used.

Upcoming release
----------------

* Version 2 will introduces a focus on simulating particles in
  optical traps rather than just focussing on calculating optical
  forces and torques.  The plan is also to introduce geometric
  optics and other methods not requiring a T-matrix.  The toolbox
  will be more automated and include a graphical user interface.
* Version 1.6 we may move Beams and Tmatrices to a beams and tmatrix
  sub-package.  We might also add drag calculation codes.

Licence
-------

Except where otherwise noted, this toolbox is made available under the
Creative Commons Attribution-NonCommercial 4.0 License.
For full details see LICENSE.md.
For use outside the conditions of the license, please contact us.
The toolbox includes some third-party components, information about
these components can be found in the documentation and corresponding
file in the thirdparty directory.

This version of the toolbox can be referenced by citing the following paper

> T. A. Nieminen, V. L. Y. Loke, A. B. Stilgoe, G. Knöner, A. M. Branczyk, N. R. Heckenberg, and H. Rubinsztein-Dunlop,
> "Optical tweezers computational toolbox",
> [Journal of Optics A 9, S196-S203 (2007)](http://iopscience.iop.org/1464-4258/9/8/S12/)

or by directly citing the toolbox

> I. C. D. Lenton, T. A. Nieminen, V. L. Y. Loke, A. B. Stilgoe,
> Y. Hu, G. Knöner, A. M. Branczyk, N. R. Heckenberg, and H. Rubinsztein-Dunlop,
> "Optical tweezers toolbox", https://github.com/ilent2/ott

and the respective Bibtex entry

```latex
@misc{Lenton2020,
  author = {Lenton, Isaac C. D. and Nieminen, Timo A. and Loke, Vincent L. Y. and Stilgoe, Alexander B. and Y. Hu and Kn{\ifmmode\ddot{o}\else\"{o}\fi}ner, Gregor and Bra{\ifmmode\acute{n}\else\'{n}\fi}czyk, Agata M. and Heckenberg, Norman R. and Rubinsztein-Dunlop, Halina},
  title = {Optical Tweezers Toolbox},
  year = {2020},
  publisher = {GitHub},
  howpublished = {\url{https://github.com/ilent2/ott}},
  commit = {A specific commit or version (optional)}
}
```

Contact us
----------

The best person to contact for inquiries about the toolbox or lincensing
is [Isaac Lenton](mailto:isaac.lenton@ist.ac.at)

File listing
------------

```
README.md     - Overview of the toolbox (this file)
LICENSE.md    - License information for OTSLM original code
AUTHORS.md    - List of contributors (pre-GitHub)
CHANGES.md    - Overview of changes to the toolbox
TODO.md       - Changes that may be made to the toolbox
thirdparty/   - Third party licenses (multiple files)
examples/     - Example files showing different toolbox features
tests/        - Unit tests to verify toolbox features function correctly
docs/         - Sphinx documentation for the project
+ott/         - The toolbox
```

The +ott package, as well as tests/ and examples/ directories
and sub-directories contain Contents.m files which list the files
and packages in each directory.
These files can be viewed in Matlab by typing `help ott`
or `help ott.subpackage`.

