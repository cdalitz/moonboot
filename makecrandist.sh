#!/bin/sh

#
# Script for creating CRAN package
# necessary due to different handling of LICENSE file on CRAN
#

BUILDDIR=build

mkdir $BUILDDIR
rsync -av ./ $BUILDDIR/ --exclude $BUILDDIR --exclude makecrandist.sh
mv $BUILDDIR/lic_cran $BUILDDIR/LICENSE
R CMD build $BUILDDIR
rm -rf $BUILDDIR
