# Digest2
Version 0.17

## Contents

1. [Program overview](#markdown-header-)
2. [Building from source](#markdown-header-_1)
3. [Obtaining a solar system model](#markdown-header-_2)
4. [Command line usage](#markdown-header-_3)
5. [Configuring file locations](#markdown-header-_4)
6. [File formats](#markdown-header-_5)
7. [Reproducing MPC output](#markdown-header-_6)
8. [Algorithm outline](#markdown-header-_7)

##
## 1.  Program overview

Digest2 uses statistical ranging techniques to compute chances that an
object is of various orbit classes.  Input is a file of 80 column MPC-format
observations, with at least two observations per object.  Output is orbit
class scores for each object.

The MPC observation format is an ASCII encoded text format documented at
http://www.minorplanetcenter.net/iau/info/OpticalObs.html.

The program is provided as C source code and an example make file.
The make file may require minor modifications for your C compiler.

Sample run:

Here are a few observations of known NEOs (with made up designations.)

````
     NE00030  C2004 09 16.15206 16 13 11.57 +20 52 23.7          21.1 Vd     291
     NE00030  C2004 09 16.15621 16 13 11.34 +20 52 16.8          20.8 Vd     291
     NE00030  C2004 09 16.16017 16 13 11.13 +20 52 09.6          20.7 Vd     291
     NE00199  C2007 02 09.24234 06 08 06.06 +43 13 26.2          20.1  c     704
     NE00199  C2007 02 09.25415 06 08 05.51 +43 13 01.7          20.1  c     704
     NE00199  C2007 02 09.26683 06 08 04.80 +43 12 37.5          19.9  c     704
     NE00269  C2003 01 06.51893 12 40 50.09 +18 27 46.9          21.4 Vd     291
     NE00269  C2003 01 06.52850 12 40 50.71 +18 27 46.1          21.8 Vd     291
     NE00269  C2003 01 06.54359 12 40 51.68 +18 27 42.5          21.9 Vd     291
````

You put them in a file, say fmo.obs, then type "digest2 fmo.obs" and 
get the following output:

```
Digest2 version 0.17 -- Released May 15 2015 -- Compiled May 15 2015
Desig.    RMS Int NEO N22 N18 Other Possibilities
NE00030  0.15 100 100  36   0
NE00199  0.56  98  97  17   0 (MC 2) (JFC 1)
NE00269  0.42  18  18   3   0 (MC 9) (Hun 4) (Pho 27) (MB1 <1) (Han <1) (MB2 30) (MB3 12) (JFC 1)
```

Only considering these short arcs, digest2 predicts that the first two
objects are almost certain to be NEOs.  The last one, with a NEO score of
only 18 shows little chance of being a NEO.  Pho 27, MB2 30, and MB3 12
show that it is much more likely to have a Phocaea-like or main belt orbit.
The scores may not add to 100 because they are computed independently and
because all possible orbits are not categorized.  Orbit classes are based
on orbital elements and are either known dynamic populations like "Hungarias"
or popular classifications like "NEO."

The RMS figure is a root-mean-square of residuals in arc seconds of the
observations against a great circle fit.  A high RMS could indicate either bad
observations or significant great circle departure over the observation arc.
It can thus be used as a quick check that scores are meaningful.  If the RMS
is low, scores are meaningful.  Scores are still meaningful if the RMS is high
due to great circle departure, but digest2 offers no way of distinguishing this
case from one where observations are bad.

##
## 2.  Building from source

Source files are available at http://bitbucket.org/mpcdev/digest2/downloads.
The downloads page has three tabs, "Downloads," "Tags," and "Branches."
Click Tags.  Read across the row the latest version; on the right are options
"zip," "gz," "bz2."  Click one of them and your browser should download the
source archive.  Alternatively, with a right click your browser may give you
the option to name the archive file.  By default the archive may be named
with a commit hash.

In any case, you should be able to unpack the archive with a command like
`tar xjf <filename>` on most Unix-like systems.

Cd into the directory created by unpacking the archive and you should find
the following files:

    README.md     this file
    LICENSE       declaration that software is public domain
    Makefile      script for building executable
    digest2.c     C source and header files
    d2cli.c
    d2model.c
    d2modelio.c
    d2math.c
    d2mpc.c
    digest2.h
    d2model.h

Type make, and the digest2 executable should be built.

There is no additional installation command.  The executable is built in the
current directory and is usable as it is.  If you would like to manually
relocate it, read below under configuration.

##
## 3.  Obtaining a solar system model

Digest2 requires a solar system model which is obtained separately from the program.  The file digest2.model.csv can be downloaded from the [download page](https://bitbucket.org/mpcdev/digest2/downloads) of this repository.  Alternatively, if you wish to duplicate exactly the results from the MPC [Neo Rating](http://mpc.cfa.harvard.edu/iau/NEO/PossNEO.html) page, follow instructions there to download an archive of the model and associated configuration file currently in use by the MPC.

##
## 4.  Command line usage

If you successfully build a `digest2` executable and downloaded a `digest2.model.csv` to the same directory, you should be able to run it from this directory.  Invoking the program without command line
arguments (or with invalid arguments) shows this usage prompt.

```
Usage: digest2 [options] <obs file>    score observations in file
       digest2 [options] -             score observations from stdin
       digest2 -m <binary model file>  generate binary model from CSV
       digest2 -h or --help            display help and quick reference
       digest2 -v or --version         display version and copyright

Options:
       -c or --config <config file>
       -m or --model <binary model file>
       -o or --obscodes <obscode file>
       -p or --config-path <path>
       -u or --cpu <n-cores>
```

The help information lists a quick reference to keywords and orbit classes
allowed in the configuration file.  The configuration file is explained
below under File Formats.

The options config, model, obscodes and config-path are described in the
next section.

The cpu option allows you to specify the number of threads to use.
The default is the number returned by the C function sysconf, typically
the total number of cores in the computer.

If you create `fmo.obs` as described above in section 1, program overview,
the trial run `digest2 fmo.obs` should give results similar to those shown above.

If you had a successful trial run, digest2 will have created two additional files.
It will have accessed the MPC web site and downloaded observatory code data
and it will have written a binary form of the model file.  Digest2 downloads
observatory code data only if finds it missing.  Similarly it writes a new
binary model file only if the csv file is updated.


##
## 5.  Configuring file locations

The digest2 executable can be copied to another location, a bin directory
for example, and you can access it as you would any other binary executable.
In this case, you can use the -p command line option to specify a path
where digest2 can find its associated files, `digest2.config`, `digest2.obscodes`,
`digest2.model`, and `digest2.model.csv`.  (A shell alias can be useful here.)

For greater control, these files can be specified individually:

	File               Command line option
	digest2.obscodes   -o
	digest2.model      -m
	digest2.config     -c

A configuration file is required to be present if -c is used.

If you specify -p in combination with -c, -o, or -m, the path specified
with the -c, -o, or -m option takes precedence.  That is, the path specified
with -p is not joined with with a file name specified with -c, -o, or -m.

The -o and -m options can also be used in a form of the digest2 command without
input observations.  With -o, the action is to get a fresh copy of obscode data
from the MPC web site and store it to the specified file.  With -m, the action
is to read `digest2.model.csv` and write the binary equivalent to the specified
file.

##
## 6.  File formats

Observations, whether supplied in a file or through stdin, should contain
observations in the MPC 80 column observation format and nothing else.
The observations should be sorted first by designation and then by time
of observation, and there should be at least two observations of each object.

digest2.obscodes is a text file containing observatory codes in the standard
MPC format.

digest2.model.csv is a comma separated text file containing the solar system model
used by digest2.

digest2.model is a binary encoding of digest2.model.csv

digest2.config, the optional configuration file, is a text file with a simple
format.  Empty lines and lines beginning with # are ignored.  Other lines must
contain either a keyword or an orbit class.

Allowable keywords:

    headings
    noheadings
    rms
    norms
    raw
    noid
    repeatable
    random
    obserr
    poss

Headings and the rms column can be turned off if desired.

Keywords raw and noid determine the score produced as described below
under Algorithm outline.  The default is noid.  If both keywords are present,
both scores are output.

The keywords repeatable and random determine if program output is
strictly repeatable or can vary slightly from one run to the next.
The program uses a Monte Carlo method.  By default, the pseudo random
number generator is seeded randomly.  When the keyword repeatable is
used, it is reseeded with a constant value for each tracklet, yielding
repeatable scores.

Keyword obserr specifies the amount of observational error that the algorithm
should allow for.  It is specified in arc seconds as in,

```
obserr=0.7
```

The default, if no obserr is specified, is 1.0 arc seconds.
Obserr may be specified for individual observatory codes as in,

```
obserrF51=.3
obserr 704 = 1
```

As shown, white space is optional.

The keyword poss specifies to output the "Other Possibilities" column.
By default, other possibilities are suppressed if orbit classes are
explicitly specified.

Orbit classes:

    Abbr.  Long Form
    ---    -------------
    Int    MPC interest.
    NEO    NEO(q < 1.3)
    N22    NEO(H <= 22)
    N18    NEO(H <= 18)
    MC     Mars Crosser
    Hun    Hungaria gr.
    Pho    Phocaea group
    MB1    Inner MB
    Pal    Pallas group
    Han    Hansa group
    MB2    Middle MB
    MB3    Outer MB
    Hil    Hilda group
    JTr    Jupiter tr.
    JFC    Jupiter Comet

Listing an orbit class limits scoring to only the listed classes.
Other possibilities are not computed or listed.  Either abbreviations or
long forms may be used.  In any case they must be spelled exactly as
shown.

Example 1:

```
Int
Neo
N22
N18
poss
```

This is equivalent to default program behavior without a config file.

Example 2:

```
# just three
NEO
Hun
JTr
```

program output is

```
Digest2 version 0.17 -- Released May 15 2015 -- Compiled May 15 2015
Desig.    RMS NEO Hun JTr
NE00030  0.15 100   0   0
NE00199  0.56  97   0   0
NE00269  0.42  18   4   0
```

The program runs considerably faster in this case, as it computes scores for
only these three classes and not all possible classes.

Example 3:

```
noheadings
norms
N22
```

output:

```
NE00030  37
NE00199  18
NE00269   3
```

This might be useful for generating results to be analyzed by another program.

Known quirk:  The order of designations may not match the input.

##
## 7.  Reproducing MPC output

Digest2 output scores can vary somewhat due to three causes: 1) the random sampling performed by the algorithm, 2) the configured observational error allowances, and 3) the asteroid catalog on which the solar system model is based.  These sources of variability can be eliminated by configuring repeatable mode, configuring observational error allowances, and using a reference model file.  The MPC is now running digest2 configured in repeatable mode and providing a download of the reference model and configuration file.  The [Neo Rating](http://mpc.cfa.harvard.edu/iau/NEO/PossNEO-test.html) page has a link to download an archive containing the model and configuration file in use.

To reproduce MPC output, you should check that you are running the same version of digest2, that you have downloaded the model and configuration from the MPC, and that you are using these files.

You should follow these steps only if you feel a compelling reason to reproduce MPC scores.  In general, the recommended practice is to run digest2 in random mode (the default) with observational error allowances you feel are appropriate, and with any reasonably current model.

##
## 8.  Algorithm outline

1.  For each object, the program computes a motion vector from the
first and last observation, and computes a V magnitude from whatever
photometry is provided.

2.  It then generates a number of orbits that are consistent with
the computed motion, complete with absolute magnitudes consistent with
the computed V magnitude.

3.  Each orbit is located within a binned, or histogram, model of the
solar system.  The model is binned in the four dimensions of q, e, i, and H.
As the bin is determined for each orbit, a tag is set for that bin.
Additionally, each orbit is tested for each configured class and a separate
tag is set, by bin, for each class.

4.  A search algorithm is used to generate orbits covering the entire range
of possible orbits, and tag corresponding possible bins.  As the algorithm
generates variant orbits, it checks if the orbits are yielding new bin
taggings, either in general or of specific orbit classes.  The algorithm
terminates after reaching a point of diminishing returns in finding bins.

5.  The histogram contains object population counts by bin.  For each bin
there is a population of each orbit class, and also a complete population
count.  After orbit search, the sum of complete population of tagged bins
represents the number of possible candidate objects in the solar system.
The population sum of tagged bins of a specified class represents the number
of possible candidates of that class.  The class sum divided by the complete
sum represents the fraction of candidate objects that are of that class.
This fraction is multiplied by 100 and output as the "raw" percentage.

6.  No-ID scores are computed similarly with a parallel histogram.
In it, population counts are not of the expected complete population of the
solar system, but of the expected yet-undiscovered population.  This
population is computed by reducing the modeled complete population by known
discoveries.  As the intended context of the no-ID score is after attempted
object identification, the selected known population is that which is readily
identifiable.  The current criteria used is sky uncertainty < 1' arc.
The uncertainty parameter selected for this comparison is field 24 of
astorb.dat, which is a peak ephemeris uncertainty over a 10 year period.