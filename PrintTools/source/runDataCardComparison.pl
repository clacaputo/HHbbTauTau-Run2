#!/usr/bin/perl 
# -w
$scriptsPath=".";
$workingPath="PrintTools/source";


sub compare {

    $channel=$_[0];
    $year=$_[1];

    $group1=$_[2];
    $path1=$_[3];
    $file1=$_[4];

    $group2=$_[5];
    $path2=$_[6];
    $file2=$_[7];
    
    system("rm -f ${channel}${year}_${group1}_${group2}_DataCardDiff.pdf");
    system("root -b  ${workingPath}/compareDataCard.C\\\(\\\"${channel}\\\",\\\"${year}\\\",\\\"${group1}\\\",\\\"${path1}\\\",\\\"${file1}\\\",\\\"${group2}\\\",\\\"${path2}\\\",\\\"${file2}\\\"\\\,\\${maxbin}\\\)");

}


#### IT 
$groupCERN="IT";
$pathTauTau2012CERN="/Users/Tita/Desktop/analysis_HH_bbTauTau/Limits/auxiliaries/shapes/Italians";
$fileTauTau2012CERN="htt_et.inputs-Hhh-8TeV_m_sv";



#### ImperialCollege - IC 
$groupMIT="IC";
$pathTauTau2012MIT="/Users/Tita/Desktop/analysis_HH_bbTauTau/Limits/auxiliaries/shapes/Imperial";
$fileTauTau2012MIT="htt_et.inputs-Hhh-8TeV_m_sv";


#### SET BIN MAXIMUM
$maxbin = 0 ;  ## 0 means full range, 26 means up to 350

compare("eleTau","2012",$groupCERN,$pathTauTau2012CERN,$fileTauTau2012CERN,$groupMIT,$pathTauTau2012MIT,$fileTauTau2012MIT,$maxbin);







