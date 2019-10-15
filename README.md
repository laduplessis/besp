# The Bayesian Epoch Sampling Skyline Plot

This [BEAST 2](https://beast2.org) package contains implementations of the the Bayesian Skyline Plot (BSP) and the Bayesian Epoch Sampling Skyline Plot (BESP) in BEAST2. 


## Citation

KV Parag, L du Plessis, OG Pybus (2019) _Jointly inferring the dynamics of population size and sampling intensity from molecular sequences_ [https://doi.org/10.1101/686378](https://doi.org/10.1101/686378)


## Installation

To install BESP:

1. Download and install [BEAST2](www.beast2.org) (at least version 2.5.0).
2. Launch the BEAUti application distributed with BEAST2.
3. From the File menu select "Manage Packages".
4. Click the "Package repositories" button at the bottom of the Package Manager dialog box.
5. Select "Add URL" and enter the following URL: `https://laduplessis.github.io/besp/package.xml`.
7. Click the "Done" button, then select "besp" from the list of packages.
8. Click the "Install/Upgrade" button. Once installation is complete, XML files with BESP classes can be used in BEAST2.


## Building package from source

To build this package from source, ensure that  you have the following installed:

- Java JDK v1.8
- Apache Ant v1.9 or later
- An internet connection

The internet connection is required since the build script downloads the most recent version of the BEAST 2 source to build the package against. Assuming both Java and Ant are on your execution path and your CWD is the root of this repository, simply type "ant build" from the command line to build the package. This may take up to a minute due to the script fetching the BEAST2 source, and the resulting binary will be left in the `/dist` directory. To run the unit tests, use "ant test".


## License

This software is free (as in freedom). You are welcome to use it, modify it, and distribute your modified versions provided you extend the same courtesy to users of your modified version. Specifically, it is made available under the terms of the GNU General Public License version 3.


## Acknowledgements

Work on this project was supported by 

- The European Commission Seventh Framework Programme (FP7/2007--2013)/European Research Council grant agreement 614725-PATHPHYLODYN
- The UK MRC and Department for International Development under grant reference MR/R015600/1
- The Oxford Martin School

---

Louis du Plessis, 2019