#!/usr/bin/perl
use strict;
#use diagnostics;
use Getopt::Long; 
use File::Path qw(make_path remove_tree); 
use File::Copy; 
use File::Basename; 
use Cwd 'abs_path'; 
use FindBin qw($RealBin);
use File::Spec;
use Encode;
use Getopt::Long;
use Spreadsheet::ParseExcel;
use Spreadsheet::ParseExcel::FmtUnicode;
my $USAGE = qq{
=================================================================================================
Name:
	$0
Function:
	1.create shell for the QTLseq and MutMap PipleLine
Usage:
	$0 -type QTL or MutMap -cfg rawData.info -faDir fa -faDownLoadUrl fa download Url -gff gff3 file -projID projID name -outDir Project Dir -appPathFile applicationPathConfig file
Options:
	-type	<string>	QTL or MutMap
	-cfg	<string>	config file
	-fa	<string>	fa file
	-faUrl	<string>	fa Download Url(ftp://ftp.ensemblgenomes.org/pub/release-19/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.19.dna.toplevel.fa.gz)
	-gff	<string>	gff3 file
	-anno	<string>	gene annotation file; with header; the first colum is geneID; others are its annotations
	-projID	<string>	projID name
	-outDir	<string>	project Dir
	-p1	<string>	parent 1; wild type; extreme trait A; if is reference, please input REF
	-p2	<string>	parent 2; mutant type; extreme trait B; if is reference, please input REF
	-s1	<string>	offspring 1; wild type; extreme trait A，
	-s2	<string>	offspring 2; mutant type; extreme trait B
	-appPathFile	<string>	appPath File(optional,default:/PUBLIC/software/RESEQ/Pipeline/current/BSA/QTLseq_TopV1Test/00.bin/appPath.info)
	-num	<int>	individuals numbers for pooling
	-popty <string>	population type,such as F2 or RIL 
    	-split <string> split the chr when call varaint YES or NO
	-indel <string> y or n


Version:
        v1.0;2014-04-25
=================================================================================================
};
my ($appPathFile,$projID,$cfg);
our ($type,$outDir,$fa,$gffFile,$gene_anno,$faUrl,$p1,$p2,$s1,$s2,$num,$popty,$indel,$split);
GetOptions (
	"type=s" =>\$type,
	"outDir=s" => \$outDir,
	"cfg=s" => \$cfg,
	"gff=s" => \$gffFile,
	"anno=s"=>\$gene_anno,
	"projID=s" => \$projID,
	"fa=s"=> \$fa,
	"faUrl=s"=> \$faUrl,
	"p1=s"=>\$p1,
	"p2=s"=>\$p2,
	"s1=s"=>\$s1,
	"s2=s"=>\$s2,
	"num=i"=>\$num,
	"popty=s"=>\$popty,
	"indel=s"=>\$indel,
	"split=s"=>\$split,
	"appPathFile=s" =>\$appPathFile,
);
$indel||="y";
die "$USAGE" unless ( $outDir and $fa and $cfg and $projID);
#一个文库下多个lane，而excel文件里面可能会有多条记录，这样会导致重复链接操作。因此一定要处理这种情况
sub main(){
	if(checkCfg($cfg) eq "no"){
		die "\ncfg file is not avaible\n";
	}
	if(! -d $outDir){
		make_path($outDir,{verbose=>1,mode=>0755});
	}
	$outDir = abs_path($outDir);
	$fa     = abs_path($fa);
	checkIndexForFA($fa);
	$gffFile= abs_path($gffFile);
	my $appPathHashPtr                         = getApp($appPathFile);
	my ($cfgFileInfoHashPtr,
	           $sampleNameArrPtr,$rawDataList) = parserRawDataFile($cfg);
	checkSample($sampleNameArrPtr);
	createQCPipleShell($rawDataList,$appPathHashPtr);
	createQTLseqMutMapPipleShell($projID,$appPathHashPtr);
	createReleaseAndReportsShell($projID,$appPathHashPtr);
	print "";
}
sub checkIndexForFA($){
	my $fa = $_[0];
	my $bwaIndexFlag = "false";
	my $samtoolsIndexFlag = "false";
	my $warnning = "";
	if(-f "$fa.amb" and -f "$fa.ann" and -f "$fa.bwt" and -f "$fa.pac" and -f "$fa.sa"){#check bwa index
		$bwaIndexFlag = "true";
	}else{
		$warnning .= "no available bwa index!\n";
	}
	if(-f $fa.".fai"){
		$samtoolsIndexFlag = "true";
	}else{
		$warnning .= "no samtools index!\n";
	}
	if($warnning ne ""){
		warn "In this FA(softLink will be parserd by abs_path() function) path:$fa\n";
		die $warnning;
	}
}

sub createReleaseAndReportsShell($){
	our $qtl_snpindel_pipleReportsAppName;
	our $qtl_snp_pipleReportsAppName;
	our $PythonAppName;
	my $projID = $_[0];
	my $appPathHashPtr = $_[1];
	my $newReportFile = "$outDir/02.varDetect/02.CreateReports.sh";
	my $ReportsApp=$appPathHashPtr->{$qtl_snpindel_pipleReportsAppName};
	my $rawDataDir = "$outDir/01.QC/00.rawMerged/";
	my $indVarDir = "$outDir/02.varDetect/";
	my $qcDir = "$outDir/01.QC";
	my $staDir = "$outDir/02.varDetect/10.result";
	my $cfgFile = "$outDir/01.QC/sampleCleanFile";
	my $reportDir = "$outDir/03.$projID\_reports";
	our $faUrl;
	createDir($reportDir);
	if(-f $ReportsApp){
		open OUT,">$newReportFile" or die "can't create $newReportFile\n";
		if ($indel eq "y"){
			print OUT "$appPathHashPtr->{$PythonAppName} $appPathHashPtr->{$qtl_snpindel_pipleReportsAppName} --project $projID --raw $rawDataDir --ind $indVarDir --fa $fa --gff $gffFile --type $type --qcd $qcDir --out $reportDir --ftp_genome $faUrl ";
		}
		else{
			print OUT "$appPathHashPtr->{$PythonAppName} $appPathHashPtr->{$qtl_snp_pipleReportsAppName} --project $projID --raw $rawDataDir --ind $indVarDir --fa $fa --gff $gffFile --type $type --qcd $qcDir --out $reportDir --ftp_genome $faUrl ";
		}
			if(defined($p1)){
				print OUT "--p1 $p1 ";
        	}
        	if(defined($p2)){
				print OUT "--p2 $p2 ";
        	}
        	if(defined($s1)){
				print OUT "--s1 $s1 ";
        	}
        	if(defined($s2)){
				print OUT "--s2 $s2 ";
        	}
		print OUT "\n";
		close OUT;
	}
	
	
	
	
}
sub createDir($){
	my $dir = $_[0];	
	make_path($dir,{verbose=>1,mode=>0755});
	return;
	if(-d $dir){#略过，先不删除代码
		my $backUPDir = "$outDir/100.backup";
		if(! -d $backUPDir){
			make_path($backUPDir,{verbose=>1,mode=>0755});
		}
		my @tmp = split /\s+/,`date`;
		my $backDirName = basename $dir;
		my $backDirDate = join("_",@tmp);
		$backDirName = $backDirName."_".$backDirDate;
		`mv $dir $backUPDir/$backDirName`;
	}
	
}
sub backup($$){
	my $file = $_[0];
	my $dir = $_[1];
	my $backUPDir = "$outDir/100.backup/$dir";
	if(! -d $backUPDir){
		make_path($backUPDir,{verbose=>1,mode=>0755});
	}
	my @tmp = split /\s+/,`date`;
	my $backFileName = basename($file).join("_",@tmp);
	`mv $file $backUPDir/$backFileName`;
}

sub parserRawDataFile($){ 
	my $cfgFileInfoHashPtr;
	my $sampleNameArrPtr;
	my $cfgFile = $_[0];
	if($cfgFile =~ /\.xls/){
		$cfgFile = getInfoFromExcel($cfgFile);
	}
	($cfgFileInfoHashPtr,$sampleNameArrPtr) = parserRawDataFile_getFileInfo($cfgFile);
	parserRawDataFile_getFqDataPath($cfgFileInfoHashPtr,$sampleNameArrPtr);
	my $rawDataList = createRawDataDir($cfgFileInfoHashPtr,$sampleNameArrPtr,$cfgFile);
	print "";
	return ($cfgFileInfoHashPtr,$sampleNameArrPtr,$rawDataList);
}
sub createRawDataDir($$){
	my $cfgFileInfoHashPtr = $_[0];
	my $sampleNameArrPtr = $_[1];
	my $cfgFile = $_[2];
	my $dir = $outDir."/00.rawData";
	$dir = abs_path($dir);
	createDir($dir);
	my $arrPtr;
	my @items;
	my $rawDataListPath = "$outDir/rawDataPath.lst";
	my $rawDataListPath_hiseq = "$outDir/rawDataHiseqPath.lst";
	my $rawDataListHistory = "$outDir/rawDataPath.history";
	my %sampleHiseqFileHash;
	if(-f $rawDataListPath_hiseq){
		getSampleHiseqFilePathInfo(\%sampleHiseqFileHash,$rawDataListPath_hiseq);
	}
	open OUT,">$rawDataListPath";
	open RAWHISEQPATH,">>$rawDataListPath_hiseq";
	open HIS,">>$rawDataListHistory";#每次的记录都会以追加方式写入，并加上时间
	my $date = `date`;
	print HIS "#$date";
	open TMPINFO,">$outDir/tmpInfo";
	my $validN = 0;
	my $sample_index;#不同的库或者lane可能测的是同一个样品，因此不排除有相同的样品名称。
	
	foreach(@{$sampleNameArrPtr}){
		$arrPtr = $cfgFileInfoHashPtr->{$_};
		$validN++;
		$sample_index = 0;
		if($sampleHiseqFileHash{$_}){
			$sample_index = @{$sampleHiseqFileHash{$_}};
			$sample_index /= 2;		
		}
		my %checkSameFilePath;#同一个文库不同的lane数据放在相同的目录下		
		foreach my $ptr(@{$arrPtr}){
			next if($ptr->[-1] eq "-");#如果找不到路径，配置文件里面该处是“-”
			$sample_index++;
			next if(exists $checkSameFilePath{$ptr->[-1]});
			$checkSameFilePath{$ptr->[-1]}="";
			@items = `ls $ptr->[-1]/*`;
			chomp @items;
			my $fileName;
			my ($lib,$sample,$lane);
			my $pathForQCflag = "false";
			foreach my $file(@items){
				$fileName = basename $file;
				($lib,$sample,$lane) = ($ptr->[2],$ptr->[1],$ptr->[0]);
				
				
				my $sampleFullName=$sample."-".$sample_index;
				$fileName =~ s/$lib/$sampleFullName/;
				$fileName =~/(.*)_1.fq.gz$/;					
				if(defined($1)){
					print TMPINFO $1."\t".$lib."\t$sample"."\t$lane"."\n";
					print OUT "$dir\t$1\t$sample\n";
					print HIS "$dir\t$1\t$sample\n";
					$pathForQCflag="true";
				}else{
					if($fileName=~/.fq.gz/ and $fileName!~/_2.fq.gz/){
						warn "wrong format fq file name:$fileName\nexpected:xxxx_1.fq.gz or xxxx_2.fq.gz\n";
						#next;
					}
					
				}
				if(-d "$dir/$fileName"){
					#backup($file,"00.rawData");#已有的忽略
					warn "have same file name:$fileName\ncheck raw data dir and original cfg file\n";#应该无论如何不会有这种情况出现
				}else{
					symlink $file,$dir."/".$fileName;#没有的新添
				}
				if($fileName=~/fq.gz/){
					my $mtime = (stat $file)[9];
					print RAWHISEQPATH $sample."\t".$fileName."\t".$file."\t".$mtime."\n";
				}				
				
			}
			warn "can not get file pattern for qc piple:$sample\t$fileName\n" if($pathForQCflag eq "false");
			print "";
		}
	}
	die "can't find raw data path for all sample, check path in file $cfgFile \nand original config file if it is parsered from excel file\n" if($validN eq 0);
	close OUT;
	close TMPINFO;
	close HIS;
	close RAWHISEQPATH;
	checkRawDataList();
	return $rawDataListPath ;
}
sub getSampleHiseqFilePathInfo($$){
	my $rawDataListPath_hiseq=$_[1];
	my $hashPtr = $_[0];
	open FILE,"<$rawDataListPath_hiseq" or die $!;
	my $line;
	my @tmp;
	while($line=<FILE>){
		@tmp = split /\t/,$line;
		if(not exists $hashPtr->{$tmp[0]}){#sample name
			my @path;
			$hashPtr->{$tmp[0]}=\@path;
		}
		push @{$hashPtr->{$tmp[0]}},$tmp[2];#hiseq path
	}
	close FILE;
}

sub checkExistsFile(){
	
}

sub checkSample($){
	my $sampleNameArrPtr =$_[0];
	our ($p1,$p2,$s1,$s2);
	my %sample_flag;
	my $sample_count=0;
	unless((defined($p1) || defined($s1))&&(defined($p2) || defined($s2))&&(defined($s1) || defined($s2))){
		print "Error input parameters: -p1 $p1 -p2 $p2 -s1 $s1 -s2 $s2!\n";
	}
	if((defined $p1)&&($p1 ne "REF")){
		$sample_flag{$p1}="false";
		$sample_count++;
	}
	if((defined $p2)&&($p2 ne "REF")){
		$sample_flag{$p2}="false";
		$sample_count++;
	}
	if(defined $s1){
		$sample_flag{$s1}="false";
		$sample_count++;
	}
	if(defined $s2){
		$sample_flag{$s2}="false";
		$sample_count++;
	}
	if($sample_count<2){
		print "Error input parameters: -p1 $p1 -p2 $p2 -s1 $s1 -s2 $s2!\n";
	}
	foreach my $sample(@{$sampleNameArrPtr}){
		if(defined($sample_flag{$p1}) && ($p1 eq $sample)){
			$sample_flag{$p1}="true";
		}elsif(defined($sample_flag{$p2}) && ($p2 eq $sample)){
			$sample_flag{$p2}="true";
		}elsif(defined($sample_flag{$s1}) && ($s1 eq $sample)){
			$sample_flag{$s1}="true";
		}elsif(defined($sample_flag{$s2}) && ($s2 eq $sample)){
			$sample_flag{$s2}="true";
		}else{
			print "Error input sample $sample!\n";
		}
	}
	foreach my $tmp(keys %sample_flag){
		if($sample_flag{$tmp} eq "false"){ 
			print "Error input samples $tmp or parameters: -p1 $p1 -p2 $p2 -s1 $s1 -s2 $s2!\n";
		}
	}
}	

sub fillSampleHash($$$){
	my $sampleNamePrefix = $_[0];
	my $sampleArrPtr = $_[1];
	my $samplesHashPtr = $_[2];
	my @tmp;
	for (my $i=0; $i<@{$sampleArrPtr}; $i++){
		@tmp = split(/-\d+_L\d+$/,$sampleNamePrefix);
		if($sampleArrPtr->[$i] eq $tmp[0]){
			$samplesHashPtr->{$tmp[0]}=$sampleNamePrefix;
			splice(@{$sampleArrPtr},$i,1);
			return;
		}
	}
	#die("can't find proper sampleNamePrefix for list file of QCPipeline.");	
}
sub getSampleFileNamePrefix($){
	#从rawDataPath.lst读取第二列信息
	my $raw_data_list = $_[0];
	open FILE_IN,"<$raw_data_list" or die "$!";
	my $line="";
	my $sampleNamePrefix="";
	my %samples;
	my @tmp;
	our ($p1,$p2,$s1,$s2);
	my @sampleArrPtr =  ($p1,$p2,$s1,$s2);
	while($line=<FILE_IN>){
		chomp $line;
		@tmp=split(/\t+/,$line);
		$sampleNamePrefix=$tmp[-1];
		fillSampleHash($sampleNamePrefix,\@sampleArrPtr,\%samples);
	}
	close FILE_IN;
	return \%samples;
}
sub createQCPipleShell($$){
	our ($QCPipleAppName,$QCReports,$perlAppName,$CPrawdataName);
	my ($rawDataList,$appPathHashPtr)=($_[0],$_[1]);
	my $QCdir="$outDir/01.QC";
	createDir($QCdir);
	my $QCShellfile_run="$QCdir/createQCPipleShell.sh";
	#my $QCShellfile_reports="$QCdir/work_createReports.sh";
	if(-f $QCShellfile_run){
		backup($QCShellfile_run,"01.QC");
	}
	
	my $samplesHashPtr = getSampleFileNamePrefix($rawDataList);
	open OUT,">$QCShellfile_run" or return "can't create $QCShellfile_run";
	my $app = $appPathHashPtr->{$QCPipleAppName};
	my $samples_str = join(",", values %{$samplesHashPtr});
	if(!defined($app)){
		warn "no app for QC_Pipleline\n";
	}else{
		print OUT "$appPathHashPtr->{$perlAppName} $app -out $QCdir -speci $projID -lst $rawDataList -que ccx.q,crop.q -qop \"-r 150 -N 0.1 -q 33 -L 5 -p 0.5\"\n";
	}
	close OUT;
	`sh $QCShellfile_run`;
	my $QCstart="$QCdir/start.sh";
	my $cpraw = $appPathHashPtr->{$CPrawdataName};
	open START,">$QCstart";
	print START "sh $QCdir/shell/run_QC_Pipeline.sh\n";
	print START "wait\nsh byebye.sh\n";
	print START "$appPathHashPtr->{$perlAppName} $cpraw -list $rawDataList -outdir $QCdir\n";
	print START "qsub -cwd -l vf=1g $QCdir/00.rawMerged/cp.sh\n";
}
sub createQTLseqMutMapPipleShell($$){
	my $projID = $_[0];
	my $appPathHashPtr = $_[1];
	our($qtl_mutmap_pipleAppName,$Ref_calAppName,$perlAppName,$num,$popty);
	our ($p1,$p2,$s1,$s2);
	my $qtl_mutmap_pipleApp = $appPathHashPtr->{$qtl_mutmap_pipleAppName};
	my $ref_calApp = $appPathHashPtr->{$Ref_calAppName};
	my @faInfo = `$appPathHashPtr->{$perlAppName} $ref_calApp -i $fa`;
	$faInfo[4]=~/Contig\s+(\d+)\s+/;
	my $lenNoN = $1;
	$faInfo[5]=~/Scaffold\s+(\d+)\s+/;
	my $lenWithN=$1;
	if(!defined($lenNoN) or !defined($lenWithN)){
		warn "has no length info with fa file\n";
	}
	`$appPathHashPtr->{$perlAppName} $qtl_mutmap_pipleApp`;
	my $varDecDir="$outDir/02.varDetect";
	createDir($varDecDir); 
	my $cfg_example = `ls ./ExampleConfigFile.cfg`;
	my $cfg_spe = "$outDir/02.varDetect/$projID.cfg";
	open CFG_EX,"<$cfg_example" or die "can't open ExampleConfigFile.cfg from $cfg";
	my $old_info_arrPtr="";
	if(-f $cfg_spe){
		$old_info_arrPtr=getIndVarCfgFqInfo($cfg_spe);
	}
	open CFG_SPE,">$cfg_spe" or die "can't create $cfg_spe";
	my @CFG_EX= <CFG_EX>;
	my $lineNum=0;
	chomp(@CFG_EX);
	$CFG_EX[1] .= $fa;
	$CFG_EX[2] .= $gffFile;
	$CFG_EX[3] .= $gene_anno if(defined $gene_anno);
	$CFG_EX[5] .= $projID;
	$CFG_EX[6] .= $lenWithN;
	$CFG_EX[7] .= $lenNoN;
	my $dataListPtr = getDataListForIndividualPiple();
	foreach(@CFG_EX){
		$lineNum++;
		if($lineNum eq 13){
			print CFG_SPE @{$dataListPtr};
			print CFG_SPE @{$old_info_arrPtr} if($old_info_arrPtr ne "");
		}
				
		print CFG_SPE $_."\n";
	}
	close CFG_EX;
	close CFG_SPE;
	
	my $runWorkShell = "$outDir/02.varDetect/01.runPiple.sh";
    our $split;
    if (!defined $split){
        $split = "NO";
    }
	open OUT,">$runWorkShell" or die "can't create $runWorkShell\n";
	print OUT "$appPathHashPtr->{$perlAppName} $qtl_mutmap_pipleApp -cfg $cfg_spe -out $varDecDir -num $num -popty $popty -split $split -appCfg $appPathFile ";
        if(defined($p1)){
                print OUT "-p1 $p1 ";
        }
        if(defined($p2)){
                print OUT "-p2 $p2 ";
        }
        if(defined($s1)){
                print OUT "-s1 $s1 ";
        }
        if(defined($s2)){
                print OUT "-s2 $s2 ";
        }
	print OUT "\n";
	close OUT; 
}
sub getIndVarCfgFqInfo($){
	my $file = $_[0];
	my @fqInfoArr;
	open FILE,"<$file" or return "";
	my $line;
	while($line=<FILE>){
		last if ($line=~/<\/DATA>/);
		next unless ($line=~/^FASTQ/);
		$line ="#".$line;
		push @fqInfoArr,$line;
		
	}
	close FILE;
	return \@fqInfoArr if(@fqInfoArr > 0);
}
sub getDataListForIndividualPiple(){
	my $cleanDataDir = "$outDir/01.QC/01.DataFilter";
	my $sampleToLibInfoFile = "$outDir/tmpInfo";
	my @rawDataList;
	getRawDataList($sampleToLibInfoFile,\@rawDataList);
	
	if(@rawDataList == 0){
		die "can't find rawData in $outDir/00.rawData/";
	}
	my $infoFileExistFlag = "true";
	my %sampleNameHash;
	my %getLibHash;
	my %getSampleHash;
	my %getLaneHash;
	open IN,"<$sampleToLibInfoFile" or $infoFileExistFlag="false";
	if($infoFileExistFlag eq "false"){
		warn "can't find $sampleToLibInfoFile for $projID.cfg file\n";
	}else{
		while(<IN>){
			my @tmp=split;
			$sampleNameHash{$tmp[2]}="";
			$getLibHash{$tmp[0]}=$tmp[1];
			$getSampleHash{$tmp[0]} = $tmp[2];
			$getLaneHash{$tmp[0]} = $tmp[3];
		}
		close IN;
	}
	my ($sampleNumHashPtr,$sameSampleCleanFqArrPtr) = getExistsSampleInfo("$outDir/01.QC/sampleCleanFile",\%sampleNameHash);
	my @cleanDataList;
	my %tmpHash;
	my $fqPath;
	my ($lib,$sample,$fqFilePrefix);
	for(my $i=1;$i<@rawDataList;$i+=2){
		my @fq1 = split /_\d.fq.gz/,$rawDataList[$i-1];
		my @fq2 = split /_\d.fq.gz/,$rawDataList[$i];
		my @tmp = split /\//,$fq1[0];
		$fqFilePrefix = $tmp[-1];
		$fq1[0] =~ s/00.rawData/01.QC\/01.DataFilter/;
		$lib = $getLibHash{$fqFilePrefix};
		$sample = $getSampleHash{$fqFilePrefix};
		next if(exists $tmpHash{$sample});
		if(!defined($lib) || !defined($sample)){
			warn "can't find library name or sample name for $fq1[0] when getting info for $projID.cfg file\n";
			warn "use $fq1[0] instead.\n";	
			$lib = $fq1[0];
			$sample=$fq1[0];		
		}
		$sampleNumHashPtr->{$sample}++;
		my $lane = $sample."_".$sampleNumHashPtr->{$sample}."_".$getLaneHash{$fqFilePrefix};
		my $line = "FASTQ\t$sample\t$lib\t$lane\t1000";
		$fqPath = dirname $fq1[0];
		$fq2[0] ="$fqPath/${sample}_2_clean.fq.gz";
		$fq1[0] ="$fqPath/${sample}_1_clean.fq.gz";
		
		$line = $line."\t".$fq1[0]."\t".$fq2[0]."\n";
		push @cleanDataList,$line;
		$tmpHash{$sample}="";
	}
	backSampleCleanToFile(\@cleanDataList);
	push @cleanDataList,@{$sameSampleCleanFqArrPtr} if(defined($sameSampleCleanFqArrPtr));
	return \@cleanDataList;
}
sub getExistsSampleInfo($$){
	my $sampleCleanFile=$_[0];
	my $sampleNameHashPtr = $_[1];
	if(-f "$outDir/02.varDetect/01.runPiple.sh"){#不是第一次建立项目流程路径，如果缺少核心文件，则报错退出
		open IN,"<$sampleCleanFile" or die "can't find $outDir/01.QC/sampleCleanFile\nWe need it for alculate QC result.\n";
	}else{
		open IN,"<$sampleCleanFile" or return;
	}
	
	my $line;
	my %sampleNumHash;
	my @sameSampleCleanFqInfo;
	my $sampleName;
	my @tmp;
	while($line=<IN>){
		@tmp = split/\t/,$line;
		$sampleName = $tmp[1];
		if(exists $sampleNameHashPtr->{$sampleName}){
			push @sameSampleCleanFqInfo,$line;			
		}
		$sampleNumHash{$sampleName}++;		
	}
	
	close IN;
	return (\%sampleNumHash,\@sameSampleCleanFqInfo);
}
sub backSampleCleanToFile($){
	my $arrPtr = $_[0];
	open OUT,">>$outDir/01.QC/sampleCleanFile" or die "can't write $outDir/01.QC/sampleCleanFile\nWe need it for calculate QC result.\n";
	
	foreach(@{$arrPtr}){
		print OUT $_;
	}
	close OUT;
}

sub getRawDataList($){
	my $sampleToLibInfoFile=$_[0];
	my $dataListArrPtr=$_[1];
	open FILE,"<$sampleToLibInfoFile" or die "can't find $sampleToLibInfoFile\nWe need it to know which data used in 00.rawData Dir\n";
	my $line;
	my @items;
	while($line = <FILE>){
		@items=split /\t/,$line;
		push @{$dataListArrPtr},`ls $outDir/00.rawData/$items[0]*.fq.gz`;
	}
	close FILE;
	chomp @{$dataListArrPtr};
}
#检查rawData是否成对
sub checkRawDataList(){
	my @fqFiles = `ls $outDir/00.rawData/*.fq.gz`;
	if(@fqFiles%2==0){
		for(my $i=1;$i<@fqFiles;$i +=2){
			my @fq1 = split /\d.fq.gz/,$fqFiles[$i-1];
			my @fq2 = split /\d.fq.gz/,$fqFiles[$i];
			die "fq1 data is not the same sample with fq2 data:
			        $fqFiles[$i-1]
			        $fqFiles[$i]\n" if($fq1[0] ne $fq2[0]);
		}
	}else{
		die "raw data is not paired, check it in $outDir/00.rawData/\n";
	}
}
sub parserRawDataFile_getFileInfo($){
	open CFG,"<$_[0]" or die $!;
	my @sampleName;
	my %cfgFileInfoHash;
	my $line;
	my %sampleNameHash;
	my ($lane,$sample,$lib,$index,$insertSize,$readsLen,$dataProductionForOMS,$dataProduction,$path);
	while(chomp($line=<CFG>)){
		#该配置文件格式:lane(lane名) sample(样品名) Library(文库名) index(index名称) insertSize(插入片段) SequenceLength(测序长度) ProductionInConstruction(合同数据量) ProductionInThisLane(上机下单)
		#$line =~ s/\s+/\t/;
		next if($line=~/#/);
		($lane,$sample,$lib,$index,$insertSize,
		       $readsLen,$dataProductionForOMS,$dataProduction,$path) = split /\t/,$line;
		if(not exists $sampleNameHash{$sample}){
			$sampleNameHash{$sample}="";
			push @sampleName,$sample;
		}
		my @items =($lane,$sample,$lib,$index,$insertSize,
		                  $readsLen,$dataProductionForOMS,$dataProduction,$path);
		
		if(not exists $cfgFileInfoHash{$sample}){#防止同一个样品有不同的库，导致一个样品多行数据
			my @arr;
			$cfgFileInfoHash{$sample}=\@arr;
		}
		push @{$cfgFileInfoHash{$sample}},\@items;
	}
	close CFG;
	return (\%cfgFileInfoHash,\@sampleName);
}

sub parserRawDataFile_getFqDataPath($$){#通过将下机数据路径、文库名和index数字进行拼接，得到完整的fq数据路径
	my $fileInfoHashPtr = $_[0];
	my $sampleArrPtr = $_[1];
	my $arrPtr;
	my ($path,$lib,$index);
	#my ($fullPath,$tmpFullPath);
	foreach(@{$sampleArrPtr}){
		$arrPtr = $fileInfoHashPtr->{$_};
		foreach my $ptr(@{$arrPtr}){ 
			next if($ptr->[0]=~/#/);		
			($path,$lib,$index) = ($ptr->[-1],$ptr->[2],$ptr->[3]);
			my ($fullPath,$tmpFullPath)=("","");
			if(-d $path){
				$fullPath = $path."/$lib";
				if(!-d $fullPath){
					$tmpFullPath = $fullPath."/$index";
					if(!-d $tmpFullPath){
						$index =~ /(\d+)/;
						if(!defined($1)){
							warn "can not get index num, rawData has changed rule already?\n";
							push @{$ptr},"-";
							next;
						}
						$tmpFullPath = $fullPath."-$1";
						if(-d $tmpFullPath){
							$fullPath = $tmpFullPath;
							push @{$ptr},$fullPath;
						}else{
							warn "can not get data path of $_\n";
							push @{$ptr},"-";
						}
					}
				}else{
					push @{$ptr},$fullPath;
				}				
			}else{
				if($path !~ /\w+/){
					warn "path info is empty for $_:($path)\n";
				}else{
					warn "it is not an available path for $_:($path)\n";
				}
				
				push @{$ptr},"-";
			}
		}
		
	}
	
}
sub getInfoFromExcel($){
	my $file = $_[0];
	my $dataArr = parserExcelToArr($file);
	die "\n###########\ncan not get info from Excel file:$file\nYou should use \"read_\" as prefix for excel sheet name and give a red color for sample name.\n###########\n" if(!defined($dataArr) or $dataArr eq "");
	my $rawDataFile = $outDir."/rawData.info";
	open OUT,">$rawDataFile";
	my $n=0;
	foreach my $row(@{$dataArr}){
		if($n eq 0){
			print OUT "#";
			$n++;
		}
		print  OUT $row->[1]."\t".$row->[3]."\t".$row->[4]."\t".$row->[5]."\t".$row->[6]."\t".$row->[8]."\t".$row->[9]."\t".$row->[10]."\t".$row->[7]."\n";
#		print  OUT $row->[0]."\t".$row->[5]."\t".$row->[9]."\t".$row->[10]."\t".$row->[8]."\t".$row->[16]."\t".$row->[17]."\t".$row->[18]."\t".$row->[19]."\n";
		print OUT "#Lane\tSample\tLib\tIndex\tInsertSize\tSequenceLength\tProductionInConstruction\tProductionInThisLane\tPath\n" if($n eq 0);
	}
	close OUT;
	`echo -ne "#" >>$rawDataFile.history`;
	`date>>$rawDataFile.history`;
	`grep -v "#" $rawDataFile >>$rawDataFile.history`;
	return $rawDataFile;
}
sub parserExcelToArr(){
	my $file = $_[0];
	my @data;
	my $parser   = Spreadsheet::ParseExcel->new();
	my $formatter = Spreadsheet::ParseExcel::FmtUnicode->new(Unicode_Map=>"CP936");
	my $workbook = $parser->parse($file, $formatter);
	if ( !defined $workbook ) {
    	die $parser->error(), ".\n";
	}
	my $sheetNum = 0;
	for my $worksheet ( $workbook->worksheets() ) {
		my $sheet_name = $worksheet->get_name();
		next if $sheet_name !~/^read/;
		$sheetNum++;
    	my ( $row_min, $row_max ) = $worksheet->row_range();
    	my ( $col_min, $col_max ) = $worksheet->col_range();
    	my @tmpDataForDealMerge;

    	ROW:for my $row ( $row_min .. $row_max ) {
    		my $ypcell = $worksheet->get_cell($row,3);
    		next ROW unless $ypcell; 
    		my @tmpData;  		
        	for my $col ( $col_min .. $col_max ) {     
        		if($sheetNum>1){
        			next if($row eq $row_min);
        		}   		
           		my $cell = $worksheet->get_cell( $row, $col );           		
           		next unless $cell;
           		#my $format = $cell->get_format();
           		#my $cellfont = $cell->{Format}->{Font};
           		
          		my $value=encode_utf8(decode("gbk",$cell->value()));
          		if($cell eq $ypcell){
        			$value =~ s/\s+//g;
        		}
          		if($cell->{Merged} and $value eq ''){
          			push @tmpData,$tmpDataForDealMerge[$col];
          			#print $tmpDataForDealMerge[$col]."\t";
          		}else{
          			push @tmpData,$value;
          			#print $value."\t";
          			$tmpDataForDealMerge[$col]=$value;
          		}
           		

        		}
        		#目标数据所在行和非目标数据行存在合并单元格情况，因此还是需要遍历全部数据以获得相关数据
        		if($row ne 0 and $ypcell->{Format}->{Font}->{Color} ne 10){
        		#	next;
        		}
        		
        		if(@tmpData){
        			push @data,\@tmpData;
        		}
        		
    		}
		}
	return \@data if(@data ne 0);
}
sub initAppNameHash($){
	my $appPathHashPtr = $_[0];
	our $Ref_calAppName="Ref_cal";
	our $QCPipleAppName="QC_pipleApp";
	our $QCReports="QCCreatReport";
	our $CPrawdataName="CPrawdata";
	our $qtl_mutmap_pipleAppName="qtl_mutmap_pipleApp";
	our $perlAppName="perlApp";
	our $PythonAppName="pythonApp";
	our $qtl_snpindel_pipleReportsAppName="qtl_snpindel_pipleReportsApp";
	our $qtl_snp_pipleReportsAppName="qtl_snp_pipleReportsApp";
	$appPathHashPtr->{$QCPipleAppName}="";
	$appPathHashPtr->{$QCReports}="";
	$appPathHashPtr->{$CPrawdataName}="";
	$appPathHashPtr->{$qtl_mutmap_pipleAppName}="";
	$appPathHashPtr->{$Ref_calAppName}="";
	$appPathHashPtr->{$qtl_snpindel_pipleReportsAppName}="";
	$appPathHashPtr->{$qtl_snp_pipleReportsAppName}="";
	$appPathHashPtr->{$perlAppName}="";
	$appPathHashPtr->{$PythonAppName}="";
}
sub getApp($){	
	my $appPathFile = $_[0];
	if(!defined($appPathFile)){
		$appPathFile = getAppPathFile();
		die "no avaible appPath File\n" if !defined($appPathFile);
	}
	my %appPathHash;
	initAppNameHash(\%appPathHash);
	open APPFILE,"<$appPathFile" or die $!;
	my $line;
	my ($appName,$appPath);
	my @appNameCannotRealize;
	while(chomp($line=<APPFILE>)){
		($appName,$appPath) = ($line =~ /(\S+)\s*=(\S+)\s*/);
		
		if(not exists $appPathHash{$appName}){
			#warn "warnning:can not realize app name :$appName in config file:$appPathFile\n";
			push @appNameCannotRealize,$appName;
		}else{
			$appPathHash{$appName}=$appPath;
		}		
	}
	close APPFILE;
	my @appNameForAll = keys %appPathHash;
	my $noAppPipleNum = 0;
	foreach(@appNameForAll){
		$noAppPipleNum++ if ($appPathHash{$_} eq "");
	}
	if ($noAppPipleNum >0){
		warn "warnning:\"$noAppPipleNum\" piple has no app. Check your app info file.If you don't want to see this message, you can move this piple name from initAppNameHash function.\n";
		open OUT,">$outDir/noAppPiple.info" or warn "can not create warn message file noAppPiple.info.\n";
		print OUT "$0 saying:\n\t\t1).the app name I can realize but can't find in appPath config file:\n";
		foreach(@appNameForAll){
			print OUT "\t\t\t".$_."\n" if ($appPathHash{$_} eq "");
		}
		print OUT "\t\t2).some app name can't be realized and give them a spelling checking in appPath config file.\n" if (@appNameCannotRealize>0);
		foreach(@appNameCannotRealize){
			print OUT "\t\t\t".$_."\n";
		}
		print OUT "\nappPath config file is $appPathFile\n";
		close OUT;
	}
	return \%appPathHash;
}
sub checkCfg($){
	return "yes" if(-f $_[0]);
	return "no";
}
sub getAppPathFile(){
	my $path_curf = File::Spec->rel2abs(__FILE__);
	my ($vol, $dirs, $file) = File::Spec->splitpath($path_curf);
	my $appPathFile = $dirs."appPath.info";
	return $appPathFile if (-e $appPathFile);
}
main();
