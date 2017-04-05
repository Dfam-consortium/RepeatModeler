##
## Utility functions for the maintainer of the RepeatModeler
## project.  This is *NOT* intended to be distributed
## with a RepeatModeler release.
##
## Robert Hubley - Jan 2010
##
##---------------------------------------------------------------##
## Our project name
PROJ = RepeatModeler
## All project files
SRC = SequenceSimilarityMatrix.pm MultAln.pm \
      NeedlemanWunschGotohAlgorithm.pm  \
      TRF.pm TRFResult.pm TRFMask  BuildDatabase RepeatUtil.pm \
      Refiner RepeatModeler EdgeExtension SegregationAnalysis \
      RepeatClassifier RepModelConfig.pm RepModelConfig.pm.tmpl \
      SeedAlignment.pm SeedAlignmentCollection.pm \
      configure util/viewMSA.pl util/dfamConsensusTool.pl \
      util/renameIds.pl 
## Directory where the HTML perldoc should go ( for development )
DOCDIR = /usr/local/rmserver/html/apidocs/RepeatModeler
## Directory where RepeatModeler will be installed for local users
CLUSTER_DIR = /usr/local/RepeatModeler-Dev

DST = $(addprefix $(DOCDIR)/,$(SRC:%=%.html))

##---------------------------------------------------------------##

## Default operation -- just rebuild the perldocs
all: $(DST) $(DOCDIR)/index.html

$(DOCDIR)/%.html: %
	-perldoc -u -d tmp.pod $<
	-pod2html tmp.pod > $@
	-touch $@
	-rm tmp.pod

$(DOCDIR)/index.html: $(SRC)
	-rm $(DOCDIR)/index.html
	echo '<h2>'$(PROJ)'</h2>' >> $(DOCDIR)/index.html
	for X in $(SRC) ; do \
	    echo '<a href="'$$X'.html" target="classFrame">'$$X'</a><br>' >> $(DOCDIR)/index.html ; \
	    done

## Cleanup the PERL files using perltidy
tidy:
	for X in $(SRC) ; do \
	  /home/rhubley/Perl-Tidy-20031021/perltidy -b -i 2 -ci 4 -lp -pt 0 -sbt 0 $$X; \
	  perl -i -pe 's/^\s+$$/\n/g' $$X ; \
	  done


##
## Install a copy for the cluster ie. /usr/local/bin/RepeatMasker
##
install-cluster:
	-mkdir $(CLUSTER_DIR)
	for X in $(SRC) ; do \
	  cp $$X $(CLUSTER_DIR)/$$X ; \
	  perl -i -pe 's/^#\!.*perl.*/#\!\/usr\/bin\/perl/g' $(CLUSTER_DIR)/$$X ; \
	done
	-mkdir $(CLUSTER_DIR)/Matrices
	cp -R Matrices/* $(CLUSTER_DIR)/Matrices
	-mkdir $(CLUSTER_DIR)/Libraries
	cp -R Matrices/* $(CLUSTER_DIR)/Libraries
	-rm -rf $(CLUSTER_DIR)/CVS
	-rm -rf $(CLUSTER_DIR)/Matrices/CVS
	-rm -rf $(CLUSTER_DIR)/Libraries/CVS


##
## Make a tar'd copy of the program for a quick ( non-official ) distribution
##
quickdistro:
	-mkdir dist        
	-mkdir dist/Matrices
	-mkdir dist/util        
	for X in $(SRC) ; do \
	  cp $$X dist/$$X ; \
	  perl -i -pe 's/^#\!.*perl.*/#\!\/usr\/local\/bin\/perl/g' dist/$$X ; \
	  done
	cp RepModelConfig.pm dist/
	cp -R Matrices/* dist/Matrices/
	rm -rf dist/Matrices/CVS
	tar zcvf RMDIST.tar.gz dist


