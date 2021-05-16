#!/bin/bash

check_exit_status () {
  if [ ! "$2" -eq "0" ]; then
    echo "Step $1 failed with exit status $2"
    exit $2
  fi
  echo "Checked step $1"
}

read1="test_file_1"
read2="test_file_2"
sample="default_sample_name"
output="."
threads="1"
max_cov=0.95
min_cov=10
theta=0.005
D=0.01
priortype="informative"
fold="unfolded"
maxsnape=0.9
nflies=5
base_quality_threshold=25
illumina_quality_coding=1.8
minIndel=5
do_prep=1
do_snape=0
do_poolsnp=0

# Credit: https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -dps|--do_poolsnp)
    do_poolsnp=1
    shift # past argument
    ;;
    -dp|--dont-prep)
    do_prep=0
    shift # past argument
    ;;
    -bq|--base-quality-threshold)
    base_quality_threshold="$2"
    shift # past argument
    shift # past value
    ;;
    -ill|--illumina-quality-coding)
    illumina_quality_coding="$2"
    shift # past argument
    shift # past value
    ;;
    -mindel|--min-indel)
    minIndel="$2"
    shift # past argument
    shift # past value
    ;;
    -c|--cores)
    threads="$2"
    shift # past argument
    shift # past value
    ;;
    -x|--max-cov)
    max_cov="$2"
    shift # past argument
    shift # past value
    ;;
    -n|--min-cov)
    min_cov="$2"
    shift # past argument
    shift # past value
    ;;
    -h|--help)
    echo "Usage:"
    echo "  singularity run [options] <image> -h          Display this help message."
    echo "  singularity run [options] <image> <fastq_file_1_path> <fastq_file_2_path> <sample_name> <output_dir> <num_cores>"
    exit 0
    shift # past argument
    ;;
    -t|--theta)
    theta=$2
    shift # past argument
    shift # past value
    ;;
    -D)
    D=$2
    shift # past argument
    shift # past value
    ;;
    -p|--priortype)
    theta=$2
    shift # past argument
    shift # past value
    ;;
    -f|--fold)
    fold=$2
    shift # past argument
    shift # past value
    ;;
    -ms|--maxsnape)
    maxsnape=$2
    shift # past argument
    shift # past value
    ;;
    -nf|--num-flies)
    nflies=$2
    shift # past argument
    shift # past value
    ;;
    -ds|--do-snape)
    do_snape=1
    shift # past argument
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") #save it to an array
    shift
    ;;
esac
done

set -- "${POSITIONAL[@]}"

if [ $# != 4 ]
  then
    echo "ERROR: Need to supply <fastq_file_1_path> <fastq_file_2_path> <sample_name> <output_dir> as positional arguments, and no others"
    exit 1
fi

read1=$1; shift
read2=$1; shift
sample=$1; shift
output=$1; shift

if [ ! -f "$read1" ]; then
  echo "ERROR: $read1 does not exist"
  exit 1
fi

if [ ! -f "$read2" ]; then
  echo "ERROR: $read2 does not exist"
  exit 1
fi

if [ ! -d $output/$sample/ ]; then
  mkdir -p $output/$sample/
fi

if [ ! -d $output/$sample/${sample}_fastqc ]; then
  mkdir $output/$sample/${sample}_fastqc
fi

if [ ! -d $output/$sample/${sample}_fastqc/trimmed ]; then
  mkdir $output/$sample/${sample}_fastqc/trimmed
fi
if [ $do_prep -eq "1" ]; then

  fastqc $read1 $read2 -o $output/$sample/${sample}_fastqc

  cutadapt \
  -q 18 \
  --minimum-length 75 \
  -o $output/$sample/${sample}.trimmed1.fq.gz \
  -p $output/$sample/${sample}.trimmed2.fq.gz \
  -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  -O 15 \
  -n 3 \
  --cores=$threads \
  $read1 $read2

  check_exit_status "cutadapt" $?

  fastqc $output/$sample/${sample}.trimmed1.fq.gz $output/$sample/${sample}.trimmed2.fq.gz -o $output/$sample/${sample}_fastqc/trimmed

  check_exit_status "fastqc" $?

  #Automatically uses all available cores
  bbmerge.sh in1=$output/$sample/${sample}.trimmed1.fq.gz in2=$output/$sample/${sample}.trimmed2.fq.gz out=$output/$sample/${sample}.merged.fq.gz outu1=$output/$sample/${sample}.1_un.fq.gz outu2=$output/$sample/${sample}.2_un.fq.gz ihist=$output/$sample/${sample}.ihist.txt outa=$output/$sample/${sample}.adapters.fa

  check_exit_status "bbmerge" $?

  rm $output/$sample/${sample}.trimmed*

  bwa mem -t $threads -M -R "@RG\tID:$sample\tSM:sample_name\tPL:illumina\tLB:lib1" /opt/hologenome/holo_agly_b4_r2.1.fa $output/$sample/${sample}.1_un.fq.gz $output/$sample/${sample}.2_un.fq.gz | samtools view -@ $threads -Sbh - > $output/$sample/${sample}.q0_merged_un.bam
  
  samtools flagstat -@ $threads $output/$sample/${sample}.q0_merged_un.bam > $output/$sample/${sample}.flagstat_q0_merged_un.txt
  
  samtools view -@ $threads -Sbh -q 20 -F 0x100 $output/$sample/${sample}.q0_merged_un.bam > $output/$sample/${sample}.merged_un.bam
  
  #bwa mem -t $threads -M -R "@RG\tID:$sample\tSM:sample_name\tPL:illumina\tLB:lib1" /opt/hologenome/holo_agly_b4_r2.1.fa $output/$sample/${sample}.1_un.fq.gz $output/$sample/${sample}.2_un.fq.gz | samtools view -@ $threads -Sbh -q 20 -F 0x100 - > $output/$sample/${sample}.merged_un.bam
  
  rm $output/$sample/${sample}.1_un.fq.gz
  rm $output/$sample/${sample}.2_un.fq.gz

  bwa mem -t $threads -M -R "@RG\tID:$sample\tSM:sample_name\tPL:illumina\tLB:lib1" /opt/hologenome/holo_agly_b4_r2.1.fa $output/$sample/${sample}.merged.fq.gz | samtools view -@ $threads -Sbh - > $output/$sample/${sample}.q0_merged.bam
  
  samtools flagstat -@ $threads $output/$sample/${sample}.q0_merged.bam > $output/$sample/${sample}.flagstat_q0_merged.txt
  
  samtools view -@ $threads -Sbh -q 20 -F 0x100 $output/$sample/${sample}.q0_merged.bam > $output/$sample/${sample}.merged.bam
  
 #bwa mem -t $threads -M -R "@RG\tID:$sample\tSM:sample_name\tPL:illumina\tLB:lib1" /opt/hologenome/holo_agly_b4_r2.1.fa $output/$sample/${sample}.merged.fq.gz | samtools view -@ $threads -Sbh -q 20 -F 0x100 - > $output/$sample/${sample}.merged.bam

  check_exit_status "bwa_mem" $?

  rm $output/$sample/${sample}.merged.fq.gz

  java -jar $PICARD MergeSamFiles I=$output/$sample/${sample}.merged.bam I=$output/$sample/${sample}.merged_un.bam SO=coordinate USE_THREADING=true O=$output/$sample/${sample}.sorted_merged.bam

  check_exit_status "Picard_MergeSamFiles" $?

  rm $output/$sample/${sample}.q0_merged.bam
  rm $output/$sample/${sample}.q0_merged_un.bam
  
  rm $output/$sample/${sample}.merged.bam
  rm $output/$sample/${sample}.merged_un.bam

  java -jar $PICARD MarkDuplicates \
  REMOVE_DUPLICATES=true \
  I=$output/$sample/${sample}.sorted_merged.bam \
  O=$output/$sample/${sample}.dedup.bam \
  M=$output/$sample/${sample}.mark_duplicates_report.txt \
  VALIDATION_STRINGENCY=SILENT

  check_exit_status "Picard_MarkDuplicates" $?

  # From this file combine multiple technical replicates
  #rm $output/$sample/${sample}.sorted_merged.bam

  samtools index $output/$sample/${sample}.dedup.bam
  
  java -jar $GATK -T RealignerTargetCreator \
  -nt $threads \
  -R /opt/hologenome/holo_agly_b4_r2.1.fa \
  -I $output/$sample/${sample}.dedup.bam \
  -o $output/$sample/${sample}.hologenome.intervals

  check_exit_status "GATK_RealignerTargetCreator" $?

  java -jar $GATK \
  -T IndelRealigner \
  -R /opt/hologenome/holo_agly_b4_r2.1.fa \
  -I $output/$sample/${sample}.dedup.bam \
  -targetIntervals $output/$sample/${sample}.hologenome.intervals \
  -o $output/$sample/${sample}.contaminated_realigned.bam

  check_exit_status "GATK_IndelRealigner" $?

  rm $output/$sample/${sample}.dedup.bam*

  #Parse the scaffolds names of Agly_B4_v2.1 genome
  agly_scaffolds="scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 scaffold_6 scaffold_7 scaffold_8 scaffold_9 scaffold_10 scaffold_11 scaffold_12 scaffold_13 scaffold_14 scaffold_15 scaffold_16 scaffold_17 scaffold_18 scaffold_19 scaffold_20 scaffold_21 scaffold_22 scaffold_23 scaffold_24 scaffold_25 scaffold_26 scaffold_27 scaffold_28 scaffold_29 scaffold_30 scaffold_31 scaffold_32 scaffold_33 scaffold_34 scaffold_35 scaffold_36 scaffold_37 scaffold_38 scaffold_39 scaffold_40 scaffold_41 scaffold_42 scaffold_43 scaffold_44 scaffold_45 scaffold_46 scaffold_47 scaffold_48 scaffold_49 scaffold_50 scaffold_51 scaffold_52 scaffold_53 scaffold_54 scaffold_55 scaffold_56 scaffold_57 scaffold_58 scaffold_59 scaffold_60 scaffold_61 scaffold_62 scaffold_63 scaffold_64 scaffold_65 scaffold_66 scaffold_67 scaffold_68 scaffold_69 scaffold_70 scaffold_71 scaffold_72 scaffold_73 scaffold_74 scaffold_75 scaffold_76 scaffold_77 scaffold_78 scaffold_79 scaffold_80 scaffold_81 scaffold_82 scaffold_83 scaffold_84 scaffold_85 scaffold_86 scaffold_87 scaffold_88 scaffold_89 scaffold_90 scaffold_91 scaffold_92 scaffold_93 scaffold_94 scaffold_95 scaffold_96 scaffold_97 scaffold_98 scaffold_99 scaffold_100 scaffold_101 scaffold_102 scaffold_103 scaffold_104 scaffold_105 scaffold_106 scaffold_107 scaffold_108 scaffold_109 scaffold_110 scaffold_111 scaffold_112 scaffold_113 scaffold_114 scaffold_115 scaffold_116 scaffold_117 scaffold_118 scaffold_119 scaffold_120 scaffold_121 scaffold_122 scaffold_123 scaffold_124 scaffold_125 scaffold_126 scaffold_127 scaffold_128 scaffold_129 scaffold_130 scaffold_131 scaffold_132 scaffold_133 scaffold_134 scaffold_135 scaffold_136 scaffold_137 scaffold_138 scaffold_139 scaffold_140 scaffold_141 scaffold_142 scaffold_143 scaffold_144 scaffold_145 scaffold_146 scaffold_147 scaffold_148 scaffold_149 scaffold_150 scaffold_151 scaffold_152 scaffold_153 scaffold_154 scaffold_155 scaffold_156 scaffold_157 scaffold_158 scaffold_159 scaffold_160 scaffold_161 scaffold_162 scaffold_163 scaffold_164 scaffold_165 scaffold_166 scaffold_167 scaffold_168 scaffold_169 scaffold_170 scaffold_171 scaffold_172 scaffold_173 scaffold_174 scaffold_175 scaffold_176 scaffold_177 scaffold_178 scaffold_179 scaffold_180 scaffold_181 scaffold_182 scaffold_183 scaffold_184 scaffold_185 scaffold_186 scaffold_187 scaffold_188 scaffold_189 scaffold_190 scaffold_191 scaffold_192 scaffold_193 scaffold_194 scaffold_195 scaffold_196 scaffold_197 scaffold_198 scaffold_199 scaffold_200 scaffold_201 scaffold_202 scaffold_203 scaffold_204 scaffold_205 scaffold_206 scaffold_207 scaffold_208 scaffold_209 scaffold_210 scaffold_211 scaffold_212 scaffold_213 scaffold_214 scaffold_215 scaffold_216 scaffold_217 scaffold_218 scaffold_219 scaffold_220 scaffold_221 scaffold_222 scaffold_223 scaffold_224 scaffold_225 scaffold_226 scaffold_227 scaffold_228 scaffold_229 scaffold_230 scaffold_231 scaffold_232 scaffold_233 scaffold_234 scaffold_235 scaffold_236 scaffold_237 scaffold_238 scaffold_239 scaffold_240 scaffold_241 scaffold_242 scaffold_243 scaffold_244 scaffold_245 scaffold_246 scaffold_247 scaffold_248 scaffold_249 scaffold_250 scaffold_251 scaffold_252 scaffold_253 scaffold_254 scaffold_255 scaffold_256 scaffold_257 scaffold_258 scaffold_259 scaffold_260 scaffold_261 scaffold_262 scaffold_263 scaffold_264 scaffold_265 scaffold_266 scaffold_267 scaffold_268 scaffold_269 scaffold_270 scaffold_271 scaffold_272 scaffold_273 scaffold_274 scaffold_275 scaffold_276 scaffold_277 scaffold_278 scaffold_279 scaffold_280 scaffold_281 scaffold_282 scaffold_283 scaffold_284 scaffold_285 scaffold_286 scaffold_287 scaffold_288 scaffold_289 scaffold_290 scaffold_291 scaffold_292 scaffold_293 scaffold_294 scaffold_295 scaffold_296 scaffold_297 scaffold_298 scaffold_299 scaffold_300 scaffold_301 scaffold_302 scaffold_303 scaffold_304 scaffold_305 scaffold_306 scaffold_307 scaffold_308 scaffold_309 scaffold_310 scaffold_311 scaffold_312 scaffold_313 scaffold_314 scaffold_315 scaffold_316 scaffold_317 scaffold_318 scaffold_319 scaffold_320 scaffold_321 scaffold_322 scaffold_323 scaffold_324 scaffold_325 scaffold_326 scaffold_327 scaffold_328 scaffold_329 scaffold_330 scaffold_331 scaffold_332 scaffold_333 scaffold_334 scaffold_335 scaffold_336 scaffold_337 scaffold_338 scaffold_339 scaffold_340 scaffold_341 scaffold_342 scaffold_343 scaffold_344 scaffold_345 scaffold_346 scaffold_347 scaffold_348 scaffold_349 scaffold_350 scaffold_351 scaffold_352 scaffold_353 scaffold_354 scaffold_355 scaffold_356 scaffold_357 scaffold_358 scaffold_359 scaffold_360 scaffold_361 scaffold_362 scaffold_363 scaffold_364 scaffold_365 scaffold_366 scaffold_367 scaffold_368 scaffold_369 scaffold_370 scaffold_371 scaffold_372 scaffold_373 scaffold_374 scaffold_375 scaffold_376 scaffold_377 scaffold_378 scaffold_379 scaffold_380 scaffold_381 scaffold_382 scaffold_383 scaffold_384 scaffold_385 scaffold_386 scaffold_387 scaffold_388 scaffold_389 scaffold_390 scaffold_391 scaffold_392 scaffold_393 scaffold_394 scaffold_395 scaffold_396 scaffold_397 scaffold_398 scaffold_399 scaffold_400 scaffold_401 scaffold_402 scaffold_403 scaffold_404 scaffold_405 scaffold_406 scaffold_407 scaffold_408 scaffold_409 scaffold_410 scaffold_411 scaffold_412 scaffold_413 scaffold_414 scaffold_415 scaffold_416 scaffold_417 scaffold_418 scaffold_419 scaffold_420 scaffold_421 scaffold_422 scaffold_423 scaffold_424 scaffold_425 scaffold_426 scaffold_427 scaffold_428 scaffold_429 scaffold_430 scaffold_431 scaffold_432 scaffold_433 scaffold_434 scaffold_435 scaffold_436 scaffold_437 scaffold_438 scaffold_439 scaffold_440 scaffold_441 scaffold_442 scaffold_443 scaffold_444 scaffold_445 scaffold_446 scaffold_447 scaffold_448 scaffold_449 scaffold_450 scaffold_451 scaffold_452 scaffold_453 scaffold_454 scaffold_455 scaffold_456 scaffold_457 scaffold_458 scaffold_459 scaffold_460 scaffold_461 scaffold_462 scaffold_463 scaffold_464 scaffold_465 scaffold_466 scaffold_467 scaffold_468 scaffold_469 scaffold_470 scaffold_471 scaffold_472 scaffold_473 scaffold_474 scaffold_475 scaffold_476 scaffold_477 scaffold_478 scaffold_479 scaffold_480 scaffold_481 scaffold_482 scaffold_483 scaffold_484 scaffold_485 scaffold_486 scaffold_487 scaffold_488 scaffold_489 scaffold_490 scaffold_491 scaffold_492 scaffold_493 scaffold_494 scaffold_495 scaffold_496 scaffold_497 scaffold_498 scaffold_499 scaffold_500 scaffold_501 scaffold_502 scaffold_503 scaffold_504 scaffold_505 scaffold_506 scaffold_507 scaffold_508 scaffold_509 scaffold_510 scaffold_511 scaffold_512 scaffold_513 scaffold_514 scaffold_515 scaffold_516 scaffold_517 scaffold_518 scaffold_519 scaffold_520 scaffold_521 scaffold_522 scaffold_523 scaffold_524 scaffold_525 scaffold_526 scaffold_527 scaffold_528 scaffold_529 scaffold_530 scaffold_531 scaffold_532 scaffold_533 scaffold_534 scaffold_535 scaffold_536 scaffold_537 scaffold_538 scaffold_539 scaffold_540 scaffold_541 scaffold_542 scaffold_543 scaffold_544 scaffold_545 scaffold_546 scaffold_547 scaffold_548 scaffold_549 scaffold_550 scaffold_551 scaffold_552 scaffold_553 scaffold_554 scaffold_555 scaffold_556 scaffold_557 scaffold_558 scaffold_559 scaffold_560 scaffold_561 scaffold_562 scaffold_563 scaffold_564 scaffold_565 scaffold_566 scaffold_567 scaffold_568 scaffold_569 scaffold_570 scaffold_571 scaffold_572 scaffold_573 scaffold_574 scaffold_575 scaffold_576 scaffold_577 scaffold_578 scaffold_579 scaffold_580 scaffold_581 scaffold_582 scaffold_583 scaffold_584 scaffold_585 scaffold_586 scaffold_587 scaffold_588 scaffold_589 scaffold_590 scaffold_591 scaffold_592 scaffold_593 scaffold_594 scaffold_595 scaffold_596 scaffold_597 scaffold_598 scaffold_599 scaffold_600 scaffold_601 scaffold_602 scaffold_603 scaffold_604 scaffold_605 scaffold_606 scaffold_607 scaffold_608 scaffold_609 scaffold_610 scaffold_611 scaffold_612 scaffold_613 scaffold_614 scaffold_615 scaffold_616 scaffold_617 scaffold_618 scaffold_619 scaffold_620 scaffold_621 scaffold_622 scaffold_623 scaffold_624 scaffold_625 scaffold_626 scaffold_627 scaffold_628 scaffold_629 scaffold_630 scaffold_631 scaffold_632 scaffold_633 scaffold_634 scaffold_635 scaffold_636 scaffold_637 scaffold_638 scaffold_639 scaffold_640 scaffold_641 scaffold_642 scaffold_643 scaffold_644 scaffold_645 scaffold_646 scaffold_647 scaffold_648 scaffold_649 scaffold_650 scaffold_651 scaffold_652 scaffold_653 scaffold_654 scaffold_655 scaffold_656 scaffold_657 scaffold_658 scaffold_659 scaffold_660 scaffold_661 scaffold_662 scaffold_663 scaffold_664 scaffold_665 scaffold_666 scaffold_667 scaffold_668 scaffold_669 scaffold_670 scaffold_671 scaffold_672 scaffold_673 scaffold_674 scaffold_675 scaffold_676 scaffold_677 scaffold_678 scaffold_679 scaffold_680 scaffold_681 scaffold_682 scaffold_683 scaffold_684 scaffold_685 scaffold_686 scaffold_687 scaffold_688 scaffold_689 scaffold_690 scaffold_691 scaffold_692 scaffold_693 scaffold_694 scaffold_695 scaffold_696 scaffold_697 scaffold_698 scaffold_699 scaffold_700 scaffold_701 scaffold_702 scaffold_703 scaffold_704 scaffold_705 scaffold_706 scaffold_707 scaffold_708 scaffold_709 scaffold_710 scaffold_711 scaffold_712 scaffold_713 scaffold_714 scaffold_715 scaffold_716 scaffold_717 scaffold_718 scaffold_719 scaffold_720 scaffold_721 scaffold_722 scaffold_723 scaffold_724 scaffold_725 scaffold_726 scaffold_727 scaffold_728 scaffold_729 scaffold_730 scaffold_731 scaffold_732 scaffold_733 scaffold_734 scaffold_735 scaffold_736 scaffold_737 scaffold_738 scaffold_739 scaffold_740 scaffold_741 scaffold_742 scaffold_743 scaffold_744 scaffold_745 scaffold_746 scaffold_747 scaffold_748 scaffold_749 scaffold_750 scaffold_751 scaffold_752 scaffold_753 scaffold_754 scaffold_755 scaffold_756 scaffold_757 scaffold_758 scaffold_759 scaffold_760 scaffold_761 scaffold_762 scaffold_763 scaffold_764 scaffold_765 scaffold_766 scaffold_767 scaffold_768 scaffold_769 scaffold_770 scaffold_771 scaffold_772 scaffold_773 scaffold_774 scaffold_775 scaffold_776 scaffold_777 scaffold_778 scaffold_779 scaffold_780 scaffold_781 scaffold_782 scaffold_783 scaffold_784 scaffold_785 scaffold_786 scaffold_787 scaffold_788 scaffold_789 scaffold_790 scaffold_791 scaffold_792 scaffold_793 scaffold_794 scaffold_795 scaffold_796 scaffold_797 scaffold_798 scaffold_799 scaffold_800 scaffold_801 scaffold_802 scaffold_803 scaffold_804 scaffold_805 scaffold_806 scaffold_807 scaffold_808 scaffold_809 scaffold_810 scaffold_811 scaffold_812 scaffold_813 scaffold_814 scaffold_815 scaffold_816 scaffold_817 scaffold_818 scaffold_819 scaffold_820 scaffold_821 scaffold_822 scaffold_823 scaffold_824 scaffold_825 scaffold_826 scaffold_827 scaffold_828 scaffold_829 scaffold_830 scaffold_831 scaffold_832 scaffold_833 scaffold_834 scaffold_835 scaffold_836 scaffold_837 scaffold_838 scaffold_839 scaffold_840 scaffold_841 scaffold_842 scaffold_843 scaffold_844 scaffold_845 scaffold_846 scaffold_847 scaffold_848 scaffold_849 scaffold_850 scaffold_851 scaffold_852 scaffold_853 scaffold_854 scaffold_855 scaffold_856 scaffold_857 scaffold_858 scaffold_859 scaffold_860 scaffold_861 scaffold_862 scaffold_863 scaffold_864 scaffold_865 scaffold_866 scaffold_867 scaffold_868 scaffold_869 scaffold_870 scaffold_871 scaffold_872 scaffold_873 scaffold_874 scaffold_875 scaffold_876 scaffold_877 scaffold_878 scaffold_879 scaffold_880 scaffold_881 scaffold_882 scaffold_883 scaffold_884 scaffold_885 scaffold_886 scaffold_887 scaffold_888 scaffold_889 scaffold_890 scaffold_891 scaffold_892 scaffold_893 scaffold_894 scaffold_895 scaffold_896 scaffold_897 scaffold_898 scaffold_899 scaffold_900 scaffold_901 scaffold_902 scaffold_903 scaffold_904 scaffold_905 scaffold_906 scaffold_907 scaffold_908 scaffold_909 scaffold_910 scaffold_911 scaffold_912 scaffold_913 scaffold_914 scaffold_915 scaffold_916 scaffold_917 scaffold_918 scaffold_919 scaffold_920 scaffold_921 scaffold_922 scaffold_923 scaffold_924 scaffold_925 scaffold_926 scaffold_927 scaffold_928 scaffold_929 scaffold_930 scaffold_931 scaffold_932 scaffold_933 scaffold_934 scaffold_935 scaffold_936 scaffold_937 scaffold_938 scaffold_939 scaffold_940 scaffold_941 mitochondrion_genome"
  
  # Only A. glycines metters here since I didn't look for other Aphids species as potential contaminants
  # In the output file, only reads mapped to Agly_B4_v2.1 genome are kept
  samtools view -@ $threads $output/$sample/${sample}.contaminated_realigned.bam $agly_scaffolds -b > $output/$sample/${sample}.agly.bam
  
  # Here I keep a bam with all reads in case we want to analyze the symbionts data
  mv $output/$sample/${sample}.contaminated_realigned.bam  $output/$sample/${sample}.original.bam
  rm $output/$sample/${sample}.contaminated_realigned.bai

  check_exit_status "mpileup" $?

fi

if [ $do_poolsnp -eq "1" ]; then

  samtools mpileup $output/$sample/${sample}.agly.bam \
  -B \
  -Q ${base_quality_threshold} \
  -f /opt/hologenome/raw/A_glycines_b4_r2.1.fasta > $output/$sample/${sample}.agly_mpileup.txt


  python3 /opt/DEST-AglyPoolseq/mappingPipeline/scripts/Mpileup2Sync.py \
  --mpileup $output/$sample/${sample}.agly_mpileup.txt \
  --ref /opt/hologenome/raw/A_glycines_b4_r2.1.fasta.pickled.ref \
  --output $output/$sample/${sample} \
  --base-quality-threshold $base_quality_threshold \
  --coding $illumina_quality_coding \
  --minIndel $minIndel

  check_exit_status "Mpileup2Sync" $?

  #For the PoolSNP output
  python3 /opt/DEST-AglyPoolseq/mappingPipeline/scripts/MaskSYNC_snape_complete.py \
  --sync $output/$sample/${sample}.sync.gz \
  --output $output/$sample/${sample} \
  --indel $output/$sample/${sample}.indel \
  --coverage $output/$sample/${sample}.cov \
  --mincov $min_cov \
  --maxcov $max_cov \
  --te /opt/DEST-AglyPoolseq/RepeatMasker/ref/Aphis_glycines_4.v2.1.scaffolds.fa.out.gff\
  --maxsnape $maxsnape

  check_exit_status "MaskSYNC" $?

  # gzip $output/$sample/${sample}.cov
  # gzip $output/$sample/${sample}.indel

  mv $output/$sample/${sample}_masked.sync.gz $output/$sample/${sample}.masked.sync.gz
  gunzip $output/$sample/${sample}.masked.sync.gz
  bgzip $output/$sample/${sample}.masked.sync
  tabix -s 1 -b 2 -e 2 $output/$sample/${sample}.masked.sync.gz

  check_exit_status "tabix" $?

  echo "Read 1: $read1" >> $output/$sample/${sample}.parameters.txt
  echo "Read 2: $read2" >> $output/$sample/${sample}.parameters.txt
  echo "Sample name: $sample" >> $output/$sample/${sample}.parameters.txt
  echo "Output directory: $output" >> $output/$sample/${sample}.parameters.txt
  echo "Number of cores used: $threads" >> $output/$sample/${sample}.parameters.txt
  echo "Max cov: $max_cov" >> $output/$sample/${sample}.parameters.txt
  echo "Min cov $min_cov" >> $output/$sample/${sample}.parameters.txt
  echo "base-quality-threshold $base_quality_threshold" >> $output/$sample/${sample}.parameters.txt
  echo "illumina-quality-coding $illumina_quality_coding" >> $output/$sample/${sample}.parameters.txt
  echo "min-indel $minIndel" >> $output/$sample/${sample}.parameters.txt

fi

#Generate the SNAPE SYNC files
if [ $do_snape -eq "1" ]; then

  /opt/DEST-AglyPoolseq/mappingPipeline/scripts/Mpileup2Snape.sh \
  ${sample}.agly_mpileup.txt \
  $output \
  $sample \
  $theta \
  $D \
  $priortype \
  $fold \
  $nflies

  check_exit_status "Mpileup2SNAPE" $?

  gzip -f $output/$sample/${sample}.SNAPE.output.txt

  python3 /opt/DEST-AglyPoolseq/mappingPipeline/scripts/SNAPE2SYNC.py \
  --input $output/$sample/${sample}.SNAPE.output.txt.gz \
  --ref /opt/hologenome/raw/A_glycines_b4_r2.1.fasta.pickled.ref \
  --output $output/$sample/${sample}.SNAPE

  check_exit_status "SNAPE2SYNC" $?

  python3 /opt/DEST-AglyPoolseq/mappingPipeline/scripts/MaskSYNC_snape_complete.py \
  --sync $output/$sample/${sample}.SNAPE.sync.gz \
  --output $output/$sample/${sample}.SNAPE.complete \
  --indel $output/$sample/${sample}.indel \
  --coverage $output/$sample/${sample}.cov \
  --mincov $min_cov \
  --maxcov $max_cov \
  --te /opt/DEST-AglyPoolseq/RepeatMasker/ref/Aphis_glycines_4.v2.1.scaffolds.fa.out.gff \
  --maxsnape $maxsnape \
  --SNAPE

  check_exit_status "MaskSYNC_SNAPE_Complete" $?

  mv $output/$sample/${sample}.SNAPE.complete_masked.sync.gz $output/$sample/${sample}.SNAPE.complete.masked.sync.gz

  python3 /opt/DEST-AglyPoolseq/mappingPipeline/scripts/MaskSYNC_snape_monomorphic_filter.py \
  --sync $output/$sample/${sample}.SNAPE.sync.gz \
  --output $output/$sample/${sample}.SNAPE.monomorphic \
  --indel $output/$sample/${sample}.indel \
  --coverage $output/$sample/${sample}.cov \
  --mincov $min_cov \
  --maxcov $max_cov \
  --te /opt/DEST-AglyPoolseq/RepeatMasker/ref/Aphis_glycines_4.v2.1.scaffolds.fa.out.gff \
  --maxsnape $maxsnape \
  --SNAPE

  check_exit_status "MaskSYNC_SNAPE_Monomporphic_Filter" $?

  mv $output/$sample/${sample}.SNAPE.monomorphic_masked.sync.gz $output/$sample/${sample}.SNAPE.monomorphic.masked.sync.gz

  gunzip $output/$sample/${sample}.SNAPE.complete.masked.sync.gz
  bgzip $output/$sample/${sample}.SNAPE.complete.masked.sync
  tabix -s 1 -b 2 -e 2 $output/$sample/${sample}.SNAPE.complete.masked.sync.gz

  gunzip $output/$sample/${sample}.SNAPE.monomorphic.masked.sync.gz
  bgzip $output/$sample/${sample}.SNAPE.monomorphic.masked.sync
  tabix -s 1 -b 2 -e 2 $output/$sample/${sample}.SNAPE.monomorphic.masked.sync.gz

  check_exit_status "tabix" $?

  gzip $output/$sample/${sample}.agly_mpileup.txt

  echo "Maxsnape $maxsnape" >> $output/$sample/${sample}.parameters.txt
  echo "theta:  $theta" >> $output/$sample/${sample}.parameters.txt
  echo "D:  $D" >> $output/$sample/${sample}.parameters.txt
  echo "priortype: $priortype" >> $output/$sample/${sample}.parameters.txt

fi
