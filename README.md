# HADES-CrAP
**HADES** **C**o**r**elation **A**nalysis of **P**rotons

## Prequisites

- HYDRA version 7.0 or greater
- vae24

## Directory Structure

The important part of this repo looks as follows:

```
.
├── FemtoMixer/
│   ├──  EventCandidate.hxx
│   ├──  EventUtils.hxx
│   ├──  JJUtils.hxx
│   ├──  MdcWires.hxx
│   ├──  PairCandidate.hxx
│   ├──  PairUtils.hxx
│   ├──  Target.hxx
|   └──  TrackCandidate.hxx
├── Includes.h
├── LICENSE
├── loopDST/
│   ├──  analysis.cc
│   ├──  listsApr12/
|   └──  sendScript_SL.sh
├── macros/
|   └── *.cc
├── newFemtoAnalysis.cc
├── newQaAnalysis.cc
└── README.md
```

- FemtoMixer/ - Contains my classes which are used for structuring and slecting events, tracks, and pairs, as well as helper functions used in my macros
- Includes.h - HYDRA include list, this file is included in my HYDRA-related analysis macros (newFemtoAnalysis.cc and newQaAnalysis.cc).
- loopDST/ - Contains macros and file lists for submitting an analysis which uses standard DSTs.
- macros/ - Contains all of my macros, some of them are obsolete. I will remove them eventually (yeah sure I will...).
- newFemtoAnalysis.cc - My currently used macro fro runnig femtoscopic analysis.
- newQaAnalysis.cc - My currently used macro fro runnig QA analysis (a lot of duplicate code with newFemtoAnalysis.cc).
- README.md - What you're reading right now.

### 1. FemtoMixer/

- EventCandidate.hxx, TrackCandidate.hxx, PairCandidate.hxx - classes which define structures of event, track and pair of tracks respectively
- EventUtils.hxx, PairUtils.hxx - classes wchic contain functions used for the process of femtoscopic mixing, which is needed for my [mixer](https://github.com/JohnTheBlindMilkman/JJFemtoMixer)
- MdcWires.hxx - collection of functions, objects and constants used for wrapping information from HParticleMetaMatcher regarding the MDC wires
- Target.hxx - compile-time information of target selection in X, Y, and Z for Apr12 and Feb24
- JJUtils.hxx - collection of helper functions used in my macros

> [!NOTE]
> All of the classes and functions are documented using Doxygen, which is also configured here. Just run `doxygen` in this directory to generate an HTML documentation

### 2. loopDST/

Standard implementation as in any other HADES analysis.
- analysis.cc - main file for running my macros with analysis (newFemtoAnalysis.cc and newQaAnalysis.cc)
- sendScript_SL.sh - shell script used for sending the compiled analysis.cc binaries to the job farm
- listsApr12/ - directory containing paths to every DST and simulation file (day-by-day and all of them at once)

To run type `make` and once compilation is successfull run the shell script with `. sendScript_SL.sh`.

### 3. macros/

Not all of the macros listed in this directory are in current use. The ones which are consist of:

- femtoMerge1D.cc and femtoMerge3D.cc - macros for summing signal and background histograms (appends `_processed` to the end of the file)
- drawProton1DInteg.cc and drawProton3DIntegrated.cc - macros for drawing correlation function integrated over whole phase space (appends `_Integ` to the end of a file)
- drawProton1DMultiDIff.cc and drawProton3DMultiDiff.cc - macros for drawing phase space differential correlation functions (e.g. kT-differential)
- makeMomResCorrection1D.cc - macro for correcting the correlation functions for imperfect momentum resolution (appends `_momRes` to the end of the fule)
- fromHADEStoHAL.cc - macro for convering the X axis of the correlation function histograms from LCMS to PRF (appends `_forHAL` to the end of the file)
- plotHBTscaling.cc - macro for drawing the radii obtained from the fit

As an example, the process of running the macros to generate results for 1D analysis in 0-10% centrality is as follows:
1. `root -q -b femtoMerge1D.cc++g`
2. `root -q -b drawProton1DInteg.cc++g`
3. `root -q -b drawProton1MultiDiff.cc++g`
4. `hadd -f 1Dcorr_0_10_cent_Purity.root 1Dcorr_0_10_cent_Integ.root 1Dcorr_0_10_cent.root` (this step is a substitution for purity correction which I don't make)
5. `root -q -b makeMomResCorrection1D.cc++g`
6. `root -q -b fromHADEStoHAL.cc++g`
7. Run your fitter
8. `root -q -b plotHBTscaling.cc++g`