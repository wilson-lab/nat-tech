
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nat-tech

<!-- badges: start -->
<!-- badges: end -->

The goal of `nat-tech` is to provide R client utilities for streamlining
registration of light-level microscopy data to existing [template
brains](https://www.janelia.org/open-science/jrc-2018-brain-templates)
for for *D. melanogaster* neurons. This tool utilizes the [Computational
Morphometry Toolkit (CMTK)](https://www.nitrc.org/projects/cmtk/) to
write CMTK registration commands.

# Data sets

The major EM dataset at the time when this package was built was the
[hemibrain
connectome](https://www.janelia.org/project-team/flyem/hemibrain).
Connectome data can be seen using the [neuPrint
website](https://neuprint.janelia.org/help/videos?dataset=hemibrain) and
accessed programmatically in R using
[neuprintr](https://github.com/natverse/neuprintr). This tool also will
match the registered data to the corresponding neuron in the hemibrain
connectome.

# How have these registrations been performed?

In the past, light-level microscopy data could be registrated using the
[fiji-cmtk-gui](https://github.com/jefferis/fiji-cmtk-gui). This tool
automates this process more without haveing to open the GUI and also
adds functionality by co-plotting these neurons with connectome date.

In the following, we detail some of `nat-tech`’s functionality

# Example Functions

``` r
#takes a .nrrd file and plots it with its counterpart in hemibrain connectome
nrrd_to_hemibrain()

#converts a .swc neuron trace and plots it with its counterpart in hemibrain connectome
neuron_to_hemibrain()

#converts a hemibrain neuron into a .nrrd file (to compare to a light-level image image)
hemibrain_to_nrrd()

#converts a flywire neuron into a .nrrd file (to compare to a light-level image image)
flywireid_to_nrrd()

#writes the cmtk registration command based on the Registration file location
write_cmtkreg()
```

To run the pipeline

# Acknowledgements

This package was created by Emily Kellogg and Alexander Bates while in
the group of [Dr. Rachel
Wilson](https://en.wikipedia.org/wiki/Rachel_Wilson_(neurobiologist)).

# References

-   **The hemibrain connectome (hemibrain:v1.2.1)**: Scheffer, L.K., Xu,
    C.S., Januszewski, M., Lu, Z., Takemura, S.-Y., Hayworth, K.J.,
    Huang, G.B., Shinomiya, K., Maitlin-Shepard, J., Berg, S., et
    al. (2020). A connectome and analysis of the adult Drosophila
    central brain. Elife 9. [doi:
    https://doi.org/10.7554/eLife.57443](https://doi.org/10.7554/eLife.57443)

-   **JRC2018F brain and VNC templates**: Bogovic, J.A., Otsuna, H.,
    Heinrich, L., Ito, M., Jeter, J., Meissner, G.W., Nern, A.,
    Colonell, J., Malkesman, O., Ito, K., et al. (2018). An unbiased
    template of the Drosophila brain and ventral nerve cord. bioRxiv.
    [doi:
    https://doi.org/10.1101/376384](https://doi.org/10.1101/376384)
