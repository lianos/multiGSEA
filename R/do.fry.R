## ------------------------------------------------------------------------
##   r102075 | smyth@wehi.edu.au | 2015-04-07 22:56:02 -0700 (Tue, 07 Apr 2015)
##
## 8 April 2015: limma 3.23.13
##
## - new function fry(), which provides a fast version of mroast()
## when there is little heteroscedasticity between genes.
##
## - minor change to mroast() output so that FDR values are never
## smaller than p-values when midp=TRUE.
