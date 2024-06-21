# HADES-CrAP
**HADES** **C**o**r**elation **A**nalysis of **P**rotons

## General Information

Blah blah blah, something about PhD.

## Prequisites

- HYDRA trunk
- vae23
- ???

## Directory Structure

The important part of this repo looks as follows:

```
.
├── FemtoMixer/
│   ├──  EventCandidate.hxx
│   ├──  JJFemtoMixer.hxx
│   ├──  PairCandidate.hxx
|   └── TrackCandidate.hxx
├── Includes.h
├── LICENSE
├── loopDST/
│   ├──  analysis.cc
│   ├──  build/
│   ├──  jobScript_SL.sh
│   ├──  listsApr12/
│   ├──  listsApr12Custom/
│   ├──  listsFeb24/
│   ├──  Makefile
|   └── sendScript_SL.sh
├── loopRAW/
│   ├──  allParam_APR12_gen8_10072015.root
│   ├──  analysisDST
│   ├──  analysisDST_2step.cc
│   ├──  build/
│   ├──  jobScript_SL.sh
│   ├──  lists/
│   ├──  Makefile
│   ├──  out/
│   ├──  sendScript.sh
│   ├──  treeFilter.h
|   └── wrap.sh
├── macros/
|   └── *.cc
├── myCrap/
├── newFemtoAnalysis.cc
├── newQaAnalysis.cc
└── README.md
```

- FemtoMixer/ - Contains my mixer class for femtoscopic analysis, as well as classes for structuring and slectinc events, tracks, and pairs.
- Includes.h - HYDRA include list, this file is included in my HYDRA-related analysis macros (newFemtoAnalysis.cc and newQaAnalysis.cc).
- loopDST/ - Contains macros and file lists for submitting an analysis which uses standard DSTs.
- loopDST/ - Contains macros and file lists for submitting an analysis which creates custom DSTs.
- macros/ - Contains all of my macros, some of them are obsolete. I will remove them eventually (yeah sure I will...).
- myCrap/ - **Obsolete** library which was supposed to be included into HYDRA to use my femtoscopic mixing class (I think, but I'm not sure).
- newFemtoAnalysis.cc - My currently used macro fro runnig femtoscopic analysis.
- newQaAnalysis.cc - My currently used macro fro runnig QA analysis (a lot of duplicate code with newFemtoAnalysis.cc).
- README.md - What you're reading right now.

### 1. FemtoMixer/

Lorem Ipsum

### 2. loopDST/

Lorem Ipsum

### 3. loopRAW/

Lorem Ipsum

### 4. macros/

Lorem Ipsum