---
title: "CRAN-comments"
output: html_document
---

## Resubmission
This is a resubmission.  In this version I have:

* Added a reference for the method to the 'Description' field of 
the DESCRIPTION file
* Unwrapped examples previously wrapped in `\dontrun{}` and modified these so
that they can be executed in <5 seconds

## Test enviornments 
* local Windows 10, R 3.4.3
* ubuntu 14.04 (on travis-ci), R 3.4.3
* win-builder (devel)

## R CMD check results
There were no ERRORs, WARNINGs.

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Rheanna Mainzer <rheanna.mainzer@unimelb.edu.au>'
  New submission
  
  This is my first submission.

## Downstream dependencies
There are currently no downstream dependencies for this package.
