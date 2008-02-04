#!/usr/bin/perl

# acrontab command: 25 * * * * lxplus (/afs/cern.ch/user/n/nougaret/myMAD/madX/AUTOMATION/MadBuild.pl) > /afs/cern.ch/user/n/nougaret/MADbuilt.out 2>&1


use MIME::Lite; # to send e-mail

@extractedPackages = ('madX');


$pwd = `pwd`;
chop $pwd;
$extractDir = join("", $pwd, "/MadCvsExtract") ;
mkdir($extractDir, 0777);
chdir($extractDir);


$cvsDir = ":kserver:isscvs.cern.ch:/local/reps/madx" ;


foreach(@extractedPackages) {
    $pack = $_;
    print "Extract package $pack from CVS\n";
    `cvs -d $cvsDir checkout $pack`;
}

# build
chdir('./madX');

@targets = ("madx","madxp");
foreach $target (@targets){
    `make clean`;
    `rm $target`;
    `make $target`;
    $nbOfTargets = `ls $target | wc -w`;
    if ($nbOfTargets == 1) { $compilationOutcome{$target} = 'success'}
    else { $compilationOutcome{$target} = 'failure'}

}


# then send an e-mail to the FESA support team
  $msg = MIME::Lite->new(
			 From     => 'Jean-Luc.Nougaret@cern.ch',
			 To       => 'Jean-Luc.Nougaret@cern.ch',
			 Subject  => "Automated MAD Build $compilationOutcome{'madx'} for madx, $compilationOutcome{'madxp'} for madxp",
			 Data     => "This is an automated e-mail. Check report on\nhttp://isscvs.cern.ch/cgi-bin/viewcvs-all.cgi/?root=madx&sortby=date"
			);
  $msg->send;

