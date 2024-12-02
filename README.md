moonboot
========

This is an R package for m-out-of-n bootstrap that provides functions for

- computation of confidence intervals
- estimation of the scaling factor tau
- different methods for estimating m

Usage
-----

A typical usage is

    # application of estimator to the subset indices
    boot.stat <- function(indices, dat) {
      my.stat(dat[indices])
    }

    # apply m-out-of-n bootstrap
    boot.out <- mboot(x, boot.stat, m=2*sqrt(length(x)))

    # compute 95% confidence interval
    ci <- mboot.ci(boot.out, type="basic")
    print(ci)

Beware that the last command estimates the scaling factor, which can be quite
unreliable. It is thus better to provide the scaling factor if it is known,
e.g., for a root-n consistent estimator:

    ci <- mboot.ci(boot.out, tau=function(n) { n^0.5 }, type="basic")


Installation
------------

If you have cloned the repository to the directory *moonboot*, you can install
it with

    R CMD INSTALL moonboot

If you have downloaded a file release *moonboot-X.Y.Z.tar.gz*, you can install
it with

    R CMD INSTALL moonboot-X.Y.Z.tar.gz


Authors & Copyright
-------------------

 - Christoph Dalitz, 2024 <https://lionel.kr.hsnr.de/~dalitz/>
 - Felix Lögler, 2024

For licensing information, see the file LICENSE for details.
