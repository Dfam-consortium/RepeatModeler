---
name: Bug report
about: Report a problem encountered with RepeatModeler.
title: ''
labels: bug
assignees: ''

---

**Describe the issue**

A concise description of the bug, including any error messages.

**Reproduction steps**

1. Steps to reproduce the behavior, including the command lines given to the program

* and links to publicly available genome assemblies and other data files (if available).

**Log output** 

Please paste or attach any and all log output, which includes useful information including data file statistics and version numbers.  An easy way to capture this is to redirect the log output to a file e.g `RepeatModeler -database mydb >& output.log`.  The log output should include the "random seed" value at the start of the run.  This number will be necessary in order to reproduce the run exactly.

**Environment (please include as much of the following information as you can find out):**

* How did you install RepeatModeler? e.g. manual installation from repeatmasker.org, bioconda, the Dfam TE Tools container, or as part of another bioinformatics tool? 

* Which version of RepeatModeler do you have?  The output of `RepeatModeler` without any options will be a help page with the version of the program displayed at the top.

* Which version of RepeatMasker is this RepeatModeler installation using? Have you installed RepBase RepeatMasker Edition for RepeatMasker, or the full Dfam database? 

* Operating system and version. The output of `uname -a` and `lsb_release -a` can be used to find this.

**Additional context**

* Add any other context you have about the problem here. Some possible examples:
  * If an older version of RepeatModeler worked before
  * If the problem only happens with specific data files
