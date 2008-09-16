# installs test suite files

PERL_SCRIPTS = MadBuild.pl MadTest.pl MadDiff.pl MadBuildAndTest.pl MadTrigTest.pl MadWorkReport.pl
XML_DOCS = TestScenario.xml ProcessScenario.xsl
HTML_DOCS = MadTestWebStyle.css index.htm

INST_SCRIPT_FILES = $(PERL_SCRIPTS) $(XML_DOCS)
INST_WEB_FILES = $(HTML_DOCS)

# temporary locations
DELIVERY_TARGET = /afs/cern.ch/user/n/nougaret/scratch0/mad-automation 
# directory in which cron starts build
WEB_FOLDER_TARGET = /afs/cern.ch/user/n/nougaret/www/mad

install: $(INSTFILES)
	chmod a+x *.pl	# just in case
	mad_install $(INST_SCRIPT_FILES) $(DELIVERY_TARGET)
	mad_install $(INST_WEB_FILES) $(WEB_FOLDER_TARGET)
