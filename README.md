
# UVA Forest Model Enhanced

[![DOI](https://zenodo.org/badge/145996324.svg)](https://zenodo.org/badge/latestdoi/145996324)

---------------------------------

This repository holds the source code for UVAFME. UVAFME is an individual-based forest gap model which simulates the establishment, growth, and mortality of individual trees on independent gaps or "plots" of a forested landscape.

For more information about UVAFME see the [How To Guide](https://github.com/UVAFME/UVAFME_model/blob/main/Running_UVAFME.pdf) or visit the [UVAFME website](https://uvafme.github.io/).

Please contact *Adrianna Foster* at adrianna.foster@nau.edu with any questions.

## Updates :
------------------------------

### January 2021

_Major Updates_

The UVAFME TreeData type now no longer extends from the SpeciesData type. This was done in order to
reduce memory load. The TreeData type now contains a pointer to the correct SpeciesData object in each
site's SpeciesData array.

Updates to the allometric equations, including the inclusion of a species-specific input *beta* parameter.
See the Running_UVAFME document for more information.

Updates to the Soil module subroutines to improve performance in boreal North America.

Updates to number of input parameters for the UVAFME2018_specieslist.csv and the
UVAFME2018_site.csv. See the Running_UVAFME document for more information.


## How To Cite:
------------------------------

To cite this model please use:

Adrianna Foster, Jacquelyn Shuman, & Hank Shugart. (2020, June 5). UVAFME/UVAFME_model: BorealNA (Version v1.0). Zenodo. http://doi.org/10.5281/zenodo.3879074 <br/>
and <br/>
Foster et al. (2019) Ecological Modelling 409 https://doi.org/10.1016/j.ecolmodel.2019.108765
