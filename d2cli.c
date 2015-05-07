// d2cli.c
//
// Public domain.

// posix source added for fileno()
#define _POSIX_SOURCE

#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "digest2.h"

char msgVersion[] =
    "Digest2 version 0.16 -- Released February 20 2015 -- Compiled %s\n";
char msgCopyright[] =
    "Public domain.";

// stuff used for parsing command line
//-----------------------------------------------------------------------------
char sOpt[] = "hvc:m:o:p:u:";
struct option lOpt[] = {
    {"help", no_argument, 0, 'h'},
    {"version", no_argument, 0, 'v'},
    {"config", required_argument, 0, 'c'},
    {"model", required_argument, 0, 'm'},
    {"obscodes", required_argument, 0, 'o'},
    {"config-path", required_argument, 0, 'p'},
    {"cpu", required_argument, 0, 'u'},
    {0, 0, 0, 0}
};
_Bool configSpec = 0;
_Bool modelSpec = 0;
_Bool ocdSpec = 0;
_Bool pathSpec = 0;
_Bool cpuSpec = 0;
_Bool classPossible;
_Bool raw, noid;
_Bool headings, rms, repeatable;
int nClassCompute;
int nClassColumns;
int classCompute[D2CLASSES];
int classColumn[D2CLASSES];


// "constant" global data
//-----------------------------------------------------------------------------
// global data that is constant, at least after being initialized, and so
// can be shared by threads.

char *fnConfig = "digest2.config";
char *fnModel = "digest2.model";
char *fnCSV = "model.csv";
char *fnOcd = "digest2.obscodes";
char *fpConfig = "";

// some stuff just used for CSV.  could be dynamically allocated...
char line[400];
char field[40];

char msgConfig[] = "Unrecognized line in config file:  %s\n";
char msgMemory[] = "Memory allocation failed.\n";
char msgOpen[] = "Open %s failed.\n";
char msgRead[] = "Read %s failed.\n";
char msgWrite[] = "Write %s failed.\n";
char msgReadInvalid[] = "Read %s failed, contents invalid.\n";
char msgStatus[] = "Internal error:  Unexpected tracklet status.\n";
char msgThread[] = "Thread creation failed.\n";
char msgOption[] = "Unknown option: %s\n";
char msgObsErr[] = "%s\nConfig file line: %s\n";
char msgUsage[] = "\
Usage: digest2 [options] <obsfile>    score observations in file\n\
       digest2 [options] -            score observations from stdin\n\
       digest2 -m <model file spec>   generate binary model from CSV\n\n\
       digest2 -h or --help           display help and quick reference\n\
       digest2 -v or --version        display version and copyright\n\
\n\
Options:\n\
       -c or --config <config-file>\n\
       -m or --model <binary model-file>\n\
       -o or --obscodes <obscode-file>\n\
       -p or --config-path <path>\n\
       -u or --cpu <n-cores>\n";

// functions
//-----------------------------------------------------------------------------

/* fatal

print message to stdout and terminate program.
*/
void fatal(char *msg)
{
    printf(msg);
    exit(-1);
}

void fatal1(char *msg, char *arg)
{
    printf(msg, arg);
    exit(-1);
}

void fatal2(char *msg, char *arg1, char *arg2)
{
    printf(msg, arg1, arg2);
    exit(-1);
}

// TODO consistency, standard exit
void csvError(char *mod, char *class, int iq, int ie, int ii, char *heading) {
    printf("CSV error %s %s\n", mod, class);
    printf("Q e i = %g %g %g\n", qpart[iq], epart[ie], ipart[ii]);
    printf("%s field: \"%s\"\n", heading, field);
    exit(-1);
}

// scans to ',' or end of string, copying to field as long as it fits.
// returns pointer to character of start that terminated the scan, ',' or
// the terminating null.  in the case where data would overflow field,
// nothing is copied and start is returned.
char *scanField(char *start) {
    char *end = strchr(start, ',');
    if (!end) {
        end = strchr(start, 0);
    }
    int n = end-start;
    if (n > sizeof(field)-1) {
        *field = 0;
        return start;
    }
    memcpy(field, start, n);
    field[n] = 0;
    return end;
}

void readCSVClass(FILE *fcsv, double pop[QX][EX][IX][HX], char *mod, char *class) {
	char *t;
    double d;
    for (int iq = 0; iq < QX; iq++)
        for (int ie = 0; ie < EX; ie++)
            for (int ii = 0; ii < IX; ii++) {
                fgets(line, sizeof(line), fcsv);
                char *p = scanField(line);
                if (*p != ',' || strcmp(field, mod) != 0) {
                    csvError(mod, class, iq, ie, ii, "Model");
                }
                p = scanField(p+1);
                if (*p != ',' || strcmp(field, class) != 0) {
                    csvError(mod, class, iq, ie, ii, "Class");
                }
                p = scanField(p+1);
                if (*p != ',' || fabs(strtod(field, &t) - qpart[iq]) > 1e-9) {
                    csvError(mod, class, iq, ie, ii, "Q");
                }
                p = scanField(p+1);
                if (*p != ',' || fabs(strtod(field, &t) - epart[ie]) > 1e-9) {
                    csvError(mod, class, iq, ie, ii, "e");
                }
                p = scanField(p+1);
                if (*p != ',' || fabs(strtod(field, &t) - ipart[ii]) > 1e-9) {
                    csvError(mod, class, iq, ie, ii, "i");
                }
                for (int ih = 0; ih < HX; ih++) {
                    p = scanField(p+1);
                    pop[iq][ie][ii][ih] = strtod(field, &t);
                    if (*p == (ih < HX-1 ? ',' : 0)) { // need delimiter
						if (t > field) continue; // and either something parsed
						if (!*t || *t == '\n') continue; // or empty field
                    }
                    sprintf(line, "H%g", hpart[ih]);
                    csvError(mod, class, iq, ie, ii, line);
                }
            }
}

void mheader() {
    strcpy(line, "Model,Class,Q,e,i");
    int j = strlen(line);
    for (int i = 0; i < HX; i++) {
        j += sprintf(line+j, ",H%g", hpart[i]);
    }
    strcpy(line+j, "\n");
}

char *CPspec(char *fn, _Bool spec) {
    if (spec || !pathSpec) {
		return fn;
	}
    char *p = malloc(strlen(fpConfig) + strlen(fn) + 2);
    sprintf(p, "%s/%s", fpConfig, fn);
    return p;
}

FILE *openCP(char *fn, _Bool spec, char *mode) {
    return fopen(CPspec(fn, spec), mode);
}

void readCSV(struct stat *buf) {
	// TODO temp code: read csv from hardcoded path
	FILE *fcsv = openCP(fnCSV, 0, "r"); // 0 because no switch for this.
	if (!fcsv)
		fatal1(msgOpen, fnCSV);

	mheader();
    char test_line[400];
	strcpy(test_line, line);
	fgets(line, sizeof(line), fcsv);
	if (strcmp(line, test_line)) {
		// TODO.  reconsider error message.  maybe just a warning?
		puts(test_line);
		puts(line);
		fatal("CSV header\n");
	}
	readCSVClass(fcsv, modelAllSS, "All", "SS");
	readCSVClass(fcsv, modelUnkSS, "Unk", "SS");
    for (int c = 0; c < D2CLASSES; c++) {
		readCSVClass(fcsv, modelAllClass[c], "All", classAbbr[c]);
		readCSVClass(fcsv, modelUnkClass[c], "Unk", classAbbr[c]);
	}
	fstat(fileno(fcsv), buf);
	fclose(fcsv);
}

void writeModel(struct stat *csv) {
	FILE *fmod = openCP(fnModel, modelSpec, "w");
	if (!fwrite(&csv->st_size, sizeof(csv->st_size), 1, fmod) ||
		!fwrite(&csv->st_mtime, sizeof(csv->st_mtime), 1, fmod) ||
    	!fwrite(modelAllSS, sizeof modelAllSS, 1, fmod) ||
        !fwrite(modelUnkSS, sizeof modelUnkSS, 1, fmod) ||
        !fwrite(modelAllClass, sizeof modelAllClass, 1, fmod) ||
        !fwrite(modelUnkClass, sizeof modelUnkClass, 1, fmod)) {
        fatal1(msgWrite, fnModel);
        fclose(fmod);
	}
}

char *parseObsErr(char *s)
{
    regmatch_t ss[3];
    int r = regexec(&rxObsErr, s, 3, ss, 0);
    if (r == REG_ESPACE)
        fatal(msgMemory);
    if (r == REG_NOMATCH)
        return "Invalid format for obserr.";
    char *endp;
    errno = 0;
    double oe = strtod(s + ss[2].rm_so, &endp);
    if (errno != 0 || endp != s + ss[2].rm_eo)
        return "Invalid obserr value.";
    if (oe > 10)
        return "Observational error > 10 arc seconds not allowed.";
    if (ss[1].rm_eo == ss[1].rm_so) {
        obsErr = oe * arcsecrad;
        return 0;
    }
    int oc = parseCod3(s + ss[1].rm_so);
    if (oc < 0)
        return "Obscode not recognized.";
    siteTable[oc].obsErr = oe * arcsecrad;
    return 0;
}

void convertCSV(FILE* fmod) {
	fatal("not implemented");
}

// read population model.
// prefer binary, with check that it matches the CSV.
void readModelStatCSV()
{
	// start with binary, if it's not not there, options are limited.
    FILE *fmod = openCP(fnModel, modelSpec, "r");
    if (!fmod) {       // no binary,
        convertCSV(0); // try csv, with no binary fallback
		return; // success reading csv
	}

	// binary opened okay, read header
	struct stat mod;
    if (!fread(&mod.st_size, sizeof mod.st_size, 1, fmod) ||
        !fread(&mod.st_mtime, sizeof mod.st_mtime, 1, fmod)) {
		fclose(fmod);                    // binary seems corrupt
		printf(msgReadInvalid, fnModel); // give warning message
		convertCSV(0); // try replacing it, but with no binary fallback
		return; // success reading csv
	}

	// binary readable so far, stat csv and compare
	struct stat csv;
	int csv_stat_err = stat(CPspec(fnCSV, 0), &csv);
	if (!csv_stat_err ||                // csv stat okay,
		csv.st_size != mod.st_size ||   // but size or time don't match
		csv.st_mtime != mod.st_mtime)
	{
		convertCSV(fmod); // try replacing it, with fmod as fallback.
		return; // success reading csv
	}

	// at least no csv inconsistency, continue with binary
    _Bool bin_ok = fread(modelAllSS, sizeof modelAllSS, 1, fmod) &&
        fread(modelUnkSS, sizeof modelUnkSS, 1, fmod) &&
        fread(modelAllClass, sizeof modelAllClass, 1, fmod) &&
        fread(modelUnkClass, sizeof modelUnkClass, 1, fmod);
	fclose(fmod);
	if (bin_ok) {
		return; // success reading binary
	}

	// last chance
	if (!csv_stat_err) {
		convertCSV(0);
		return;        // whew, success reading csv
	}

    printf(msgReadInvalid, fnModel); // message about binary
	fatal1(msgRead, fnCSV);          // and ultimate failure
}

void readModel() {
    FILE *fmod = openCP(fnModel, modelSpec, "r");
    if (!fmod) {
		fatal1(msgRead, fnModel);
	}
	struct stat mod;
    if (!fread(&mod.st_size, sizeof mod.st_size, 1, fmod) ||
        !fread(&mod.st_mtime, sizeof mod.st_mtime, 1, fmod) ||
        !fread(modelAllSS, sizeof modelAllSS, 1, fmod) ||
        !fread(modelUnkSS, sizeof modelUnkSS, 1, fmod) ||
        !fread(modelAllClass, sizeof modelAllClass, 1, fmod) ||
        !fread(modelUnkClass, sizeof modelUnkClass, 1, fmod)) {
		fatal1(msgRead, fnModel); 
	}
	fclose(fmod);
}

void readConfig()
{
    FILE *fcfg = openCP(fnConfig, configSpec, "r");

    if (!fcfg) {
        if (!configSpec)
            // config file not required to exist unless a config file name
            // is specified on the command line
            return;

        fatal1(msgOpen, fnConfig);
    }

    char line[20];
    if (!fgets(line, sizeof(line), fcfg)) {
        printf(msgRead, fnConfig);
        if (!configSpec)
            return;
        exit(-1);
    }

    _Bool rawSpec = 0;
    _Bool classSpec = 0;
    do {
        if (line[strlen(line) - 1] == '\n')
            line[strlen(line) - 1] = 0;

        if (*line == 0 || *line == '#')
            continue;           // skip empty lines and comments

        if (!strcmp(line, "headings")) {
            headings = 1;
            continue;
        }
        if (!strcmp(line, "noheadings")) {
            headings = 0;
            continue;
        }
        if (!strcmp(line, "rms")) {
            rms = 1;
            continue;
        }
        if (!strcmp(line, "norms")) {
            rms = 0;
            continue;
        }
        if (!strcmp(line, "raw")) {
            if (!rawSpec) {
                rawSpec = 1;
                noid = 0;
            }
            raw = 1;
            continue;
        }
        if (!strcmp(line, "noid")) {
            if (!rawSpec) {
                rawSpec = 1;
                raw = 0;
            }
            noid = 1;
            continue;
        }
        if (!strcmp(line, "poss")) {
            if (!classSpec) {
                classSpec = 1;
                nClassColumns = 0;
            }
            classPossible = 1;
            continue;
        }
        if (!strcmp(line, "repeatable")) {
            repeatable = 1;
            continue;
        }
        if (!strcmp(line, "random")) {
            repeatable = 0;
            continue;
        }
        if (!strncmp(line, "obserr", 6)) {
            char *errStr = parseObsErr(line + 6);
            if (errStr)
                fatal2(msgObsErr, errStr, line);
            continue;
        }
        for (int c = 0;; c++) {
            if (c == D2CLASSES)
                fatal1(msgConfig, line);
            if (!strcmp(line, classHeading[c])
                || !strcmp(line, classAbbr[c])) {
                if (!classSpec) {
                    classSpec = 1;
                    nClassColumns = 0;
                    classPossible = 0;
                }
                classColumn[nClassColumns++] = c;
                break;
            }
        }
    }
    while (fgets(line, sizeof(line), fcfg));

    // by default, we were going to compute all classes, but if classes
    // have been specified here and poss was not specified, then only
    // the specified classes need to be computed.
    if (classSpec && !classPossible)
        for (nClassCompute = 0; nClassCompute < nClassColumns;
             nClassCompute++)
            classCompute[nClassCompute] = classColumn[nClassCompute];
}

// returns obs file specified on command line.  this can be null, "-",
// or a filespec.
char *parseCl(int argc, char **argv)
{
    while (1) {
        int ox;
        int oc = getopt_long(argc, argv, sOpt, lOpt, &ox);

        switch (oc) {
        case 'h':
            puts("\
Digest2 help\n\
\n\
Digest2 uses statistical ranging techniques on short arc astrometry to\n\
compute probabilities that observed objects are of various orbit classes.\n\
Input is a file of 80 column MPC-format observations, with at least two\n\
observations per object.  Output is orbit class scores for each object.\n\
\n\
Config file keywords:\n\
   headings\n\
   noheadings\n\
   rms\n\
   norms\n\
   raw\n\
   noid\n\
   repeatable\n\
   random\n\
   poss\n\
   obserr\n\
\n\
Orbit classes:");
            for (int c = 0; c < D2CLASSES; c++)
                printf("   %-5s %-15s\n", classAbbr[c], classHeading[c]);
            puts("\nSee README for additional information.");
            exit(0);
        case 'v':
            printf(msgVersion, __DATE__);
			puts(msgCopyright);
            exit(0);
        case 'c':
            fnConfig = optarg;
            configSpec = 1;
            break;
        case 'm':
            fnModel = optarg;
            modelSpec = 1;
            break;
        case 'o':
            fnOcd = optarg;
            ocdSpec = 1;
            break;
        case 'p':
            fpConfig = optarg;
            pathSpec = 1;
            int last = strlen(fpConfig) - 1;
            if (fpConfig[last] == '/')
                fpConfig[last] = 0;
            break;
        case 'u':
            cores = atoi(optarg);
            cpuSpec = 1;
            break;
        case -1:
            // typcally one arg should be left, the input observation file.
			// it can be "-", meaning read from stdin.
			// if no args remain, -m is required.
			if (optind == argc && modelSpec) {
                return 0;
			}
			char *fnObs = argv[optind];
            if (optind == argc - 1 && (fnObs[0] != '-' || fnObs[1] == 0)) {
				return fnObs;
            }
			// else fall through
        default:
            fatal(msgUsage);
        }
    }
}