# AuGMod

This notebook tests the functionality and features of the augmod.py file. This package includes the methods:
- cutOut : make a postage stamp size cutout centered at given coordinates
- bulkCutOut : cutOut but run over a list of coordinates given on a text file
- getPixelMask : create a pixel mask to exclude pixels whose values exceed the given level
- galfit : Create a multiextension fits file containing:
  - Original image
  - Model
  - Residual image

For more information on using galfit see [here](https://users.obs.carnegiescience.edu/peng/work/galfit/README.pdf "Galfit User's Manual")
