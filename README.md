# Building the PPREY matrix for Atlantis GOA

This is a series of markdown documents and R scripts that generate diet composition information for GOA species from a number of data sources. 

There are a number of exploratory scripts and utility files, but the scripts numbered 1-7 arrange diet data from different sources into a format that then is used by `write_pprey.R` to build the pre-calibrated PPREY matrix that goes into the `biol.prm` file for Atlantis GOA.

Refer to individual scripts labeled 1-7 for data sources. Most of the data comes from the [REEM Program](https://apps-afsc.fisheries.noaa.gov/refm/reem/webdietdata/dietdataintro.php), which is diet data from the bottom trawl surveys.

The resulting PPREY matrix is the __pre-calibration__ version, so the actual diet preferences in Atlantis GOA will be different. 