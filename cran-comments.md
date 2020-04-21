## Test environments

* OS X High Sierra 10.13.6 local install, R 3.5.3
* OS X High Sierra 10.13.6, R 3.6.2 (on Travis CI)
* win-builder (devel and release)
* rhub (Ubuntu Linux 16.04 LTS, R-release, GCC, Fedora Linux, R-devel, clang, gfortran, Windows Server 2008 R2 SP1, R-devel, 32/64 bit)

## R CMD check results

There were no ERRORs, WARNINGs, or NOTES on most environments. 

There was one WARNING on rhub Windows Server 2008 R2 SP1, R-devel, 32/64 bit ("'qpdf' is needed for checks on size reduction of PDFs"), but the build was still successful.


## Downstream dependencies

There are currently no downstream dependencies for this package.


