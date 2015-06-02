# Matching MPC results exactly

Some users have expressed interest in getting a local install of digest2
to reproduce MPC results exactly.  This is possible by using the same
configuration as the MPC.  The critical parts of the configuration are

* The program version
* The Solar System model
* Configuring a "repeatable" mode
* Configuring observational error allowances

The MPC [NEO Rating](http://mpc.cfa.harvard.edu/cgi-bin/test/neorate.cgi) page
lists the version of digest2 in use at the MPC and the date of the Solar System
model in use.  The digest2 program, the Solar System model, and copies of
the MPC configuration file can all be downloaded from this Bitbucket
repository.

After compiling according to INSTALLING.md, `digest2 -v` shows the program
version and also the date and filesize of the Solar System model.  These
should match the values on the NEO Rating page.  Included in `d2model.tar.bz2`,
downloaded from the
[Downloads page](https://bitbucket.org/mpcdev/digest2/downloads), is the file
`MPC.config`.  This file contains the configuration settings for repeatable
mode and observational error allowances as used at the MPC.  Use this file as
the digest2 config file, for example by renaming it to digest2.config, and you
should be able to reproduce exactly the results from the NEO Rating page with
the "Int class" button selected.

For example, with no configuration file, digest2 produces output as shown
in OPERATION.md, which should also be much like output from the MPC NEO Rating
page with the "All orbit classes" button selected:

```
Desig.    RMS Int NEO N22 N18 Other Possibilities
NE00030  0.15 100 100  36   0
NE00199  0.56  98  97  17   0 (MC 2) (JFC 1)
NE00269  0.42  18  18   3   0 (MC 9) (Hun 4) (Pho 27) (MB1 <1) (Han <1) (MB2 30) (MB3 12) (JFC 1)
```

With `MPC.config` however, digest2 produces output matching the NEO Rating
page with the "Int class" button selected:

```
Desig.    RMS Int
NE00030  0.15 100
NE00199  0.56  98
NE00269  0.42  24
```

With matching program version and Solar System model as shown by `digest2 -v`,
results should match exactly.
