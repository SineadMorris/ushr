## Resubmission

This is a resubmission. In response to the initial submission comments I have:

* removed "in R" from the title (DESCRIPTION file).

* edited the description to include method references (DESCRIPTION file).

* changed "Year"" to "YEAR" (LICENSE file).

  
## Test environments

* OS X High Sierra 10.13.6 local install, R 3.5.3
* OS X High Sierra 10.13.3, R 3.6.1 (on Travis CI)
* ubuntu 16.04.6, R 3.6.1  (on Travis CI)
* win-builder (devel and release)


## R CMD check results

There were no ERRORs or WARNINGs on any platform. 

There was one NOTE on win-builder:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Sinead E. Morris <sinead.morris@columbia.edu>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  Perelson (12:183)
  al (12:195)
  antiretroviral (12:71)
  biphasic (12:114)
  et (12:192)

Response to this NOTE: 
- the words antiretroviral and biphasic are the correct technical terms.
- the remaining words, Perelson et al., are correct and refer to the author list of one of the added references. 


## Downstream dependencies

There are currently no downstream dependencies for this package.

