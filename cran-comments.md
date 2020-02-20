## Test environments
* Ubuntu Linux 16.04 LTS, R-release, GCC (R-hub)
* Fedora Linux, R-devel, clang, gfortran (R-hub)
* Debian Linux, R-devel, GCC ASAN/UBSAN (R-hub)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (R-hub)
* Windows Server 2008 R2 SP1, R-release, 32/64 bit ***
* local OS X install (release and devel)

## R CMD check results

There were no ERRORs or WARNINGs.

NOTE for OS X (release and devel); Windows (devel); :

  * Non-standard files/directories found at top level:
    ‘README.Rmd’ ‘cran-comments.md’

These are commonly-used files. We include an example in our README and so we generate this with R Markdown. Submission comments are in cran-comments.md to indicate the results of check() on different platforms.

NOTES for Ubuntu Linux 16.04 LTS, R-release, GCC; Fedora Linux, R-devel, clang, gfortran; and Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  
  * New submission
  
This is a new submission
  
  * Possibly mis-spelled words in DESCRIPTION:
  NODF (6:61, 11:26)
  Nestedness (2:36)
  nestedness (7:43)

These are technical terms that are spelled correctly.

There was 1 additional NOTE for Windows (devel):

  * checking for non-standard things in the check directory ... NOTE
Found the following files/directories:
  'maxnodf-Ex_x64.Rout' 'tests_i386' 'tests_x64'
  'examples_i386' 'examples_x64' 'maxnodf-Ex_i386.Rout'

This seems to be a result of the CMD check rather than an actual problem.

## Downstream dependencies
There are currently no known downstream dependencies for this package.