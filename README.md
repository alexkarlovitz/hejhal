# hejhal
Algorithms for extending Hejhal's algorithm to infinite volume.

Hejhal's algorithm (see description in section 2.1 of [1]) attempts to compute the eigenvalue and Fourier coefficients of a Maass form. The original algorithm considered Maass forms on Fuchsian groups with *finite* covolume. In this project, we consider extensions of this algorithm to groups with infinite covolume.

## Algorithms
This folder contains PARI code (see https://pari.math.u-bordeaux.fr/) for running our extensions of Hejhal's algorithm. Files with "cusp" in the name make use of a cuspidal Fourier expansion, and files with "disk" or "flare" in the name similarly use disk or flare expansions, resp. To run this code, one needs to download PARI (https://pari.math.u-bordeaux.fr/download.html); then, you can open PARI in the terminal using the command "gp". This opens an interactive PARI session. To run algorithms in flare_Schottky.pari, you simply read the file with "\r flare_Schottky.pari" to gain access to those functions.

## Notes
This folder contains LaTeX files describing the algorithms, mathematical background, experimental results, etc.

**Warning:** Some of the notation changes across documents. Sorry in advance for any confusion!

## Visuals
This folder contains python code for creating visuals of fundamental domains in the upper half plane and disk models. These files produced the visuals in the LaTeX documents in Notes.

## References

[1] Booker, Andrew R., Andreas Strömbergsson, and Akshay Venkatesh. "Effective computation of Maass cusp forms." International mathematics research notices 2006 (2006).
