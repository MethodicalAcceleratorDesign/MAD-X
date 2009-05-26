# installs test suite files

PERL_SCRIPTS = MadBuild.pl MadTest.pl MadDiff.pl MadBuildAndTest.pl \
	MadTrigTest.pl MadWorkReport.pl MadWorkBetweenReleases.pl \
	MadWindowsCompileServer.pl MadTrigRelease.pl
 
PERL_ADMIN_TOOLS = MadTestOverwrite.pl

PYTHON_SCRIPTS = MadBuildAndTest.py Notify.py MadTrigTest.py \
	MadTestExceptions.py \
	MadBuildPy.pl MadTestPy.pl \
	MadDiff.pl 


XML_DOCS = TestScenario.xml ProcessScenario.xsl
HTML_DOCS = MadTestWebStyle.css index.htm

INST_SCRIPT_FILES = $(PERL_SCRIPTS) $(PERL_ADMIN_TOOLS) $(XML_DOCS)
INST_SCRIPT_FILES_NEW = $(PYTHON_SCRIPTS) $(XML_DOCS)

INST_WEB_FILES = $(HTML_DOCS)

# temporary locations
DELIVERY_TARGET = /afs/cern.ch/user/n/nougaret/scratch0/mad-automation 
DELIVERY_TARGET_NEW = /afs/cern.ch/user/n/nougaret/scratch1/mad-automation

# directory in which cron starts build
WEB_FOLDER_TARGET = /afs/cern.ch/user/n/nougaret/www/mad

install: $(INSTFILES)
	chmod a+x *.pl	# just in case
	chmod a+x *.py
	mad_install $(INST_SCRIPT_FILES) $(DELIVERY_TARGET)
	mad_install $(INST_WEB_FILES) $(WEB_FOLDER_TARGET)
	mad_install $(INST_SCRIPT_FILES_NEW) $(DELIVERY_TARGET_NEW)
