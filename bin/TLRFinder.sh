#!/bin/bash
set -euo pipefail

version="v1.0.1"


##############################################################
##################### SETTING VARIABLES ######################
##############################################################

#default values, unless denoted when running PAQman
assembly=""
assemblylist=""
tipsize="50000"
telomererepeat="TTAGGG"
covthreshold="0.75"
sizemin="2000"
threads="1"
window="10"
slide="5"
bootstraps="1000"


prefix="TLRFinder"
output="TLRFinder_output"
help="nohelp"

## to clean up a bunch of output from the tools in order to reduce all the unnecessary output
cleanup="yes"

while [[ $# -gt 0 ]]
do
key="$1"

case "$key" in
    -a|--assembly)
    assembly="$2"
    shift
    shift
    ;;
    -al|--assemblylist)
    assemblylist="$2"
    shift
    shift
    ;;
    -ts|--tipsize)
    tipsize="$2"
    shift
    shift
    ;;
    -tr|--telomererepeat)
    telomererepeat="$2"
    shift
    shift
    ;;
    -ct|--covthreshold)
    covthreshold="$2"
    shift
    shift
    ;;
    -sm|--sizemin)
    sizemin="$2"
    shift
    shift
    ;;
    -t|--threads)
    threads="$2"
    shift
    shift
    ;;
    -w|--window)
    window="$2"
    shift
    shift
    ;;
    -s|--slide)
    slide="$2"
    shift
    shift
    ;;
    -b|--bootstraps)
    bootstraps="$2"
    shift
    shift
    ;;
    -p|--prefix)
    prefix="$2"
    shift
    shift
    ;;
    -o|--output)
    output="$2"
    shift
    shift
    ;;
    -h|--help)
    echo "
 
    TLRFinder (version: ${version})

    TLRFinder.sh -a assembly.fa
    or
    TLRFinder.sh -al list.tsv
    
    Required inputs:
    -a | --assembly     Genome assembly in fasta format (*.fa / *.fasta / *.fna) and can be gzipped (*.gz) with bgzip
    or
    -al | --assemblylist     A tsv file containing sample names in the first column and assembly paths in the second column

    Recommended inputs:
    -ts | --tipsize     Length of contig ends to be extracted for TLR detection (Default: 50000)
    -tr | --telomererepeat       Telomeric repeat pattern (Default: TTAGGG)
    -ct | --covthreshold    The amount of coverage required for a region to be considered for TLR clustering relative to the number fo telomeres (Default = 0.75)
    -sm | --sizemin     The minimum size of a region passing the coverage threshold to be considered as a potential TLR (Default: 2000)
    
    Multiple assembly specific parameters (if using --al)
    -b | --bootstraps   Number of bootstrap tests to be performed by mashtree (Default: 1000)

    Optional parameters:
    -w | --window       Number of basepairs for window averaging for coverage (Default: 10)
    -s | --slide        Number of basepairs for the window to slide for coverage (Default: 5)
    -p | --prefix       Prefix for output (Default: TLRFinder)
    -o | --output       Name of output folder for all results (default: TLRFinder_output)
    -h | --help         Print this help message
    "
    exit
    ;;
    *)  # catch invalid args
    echo "ERROR: Unknown option: '$1'"
    echo "Run 'TLRFinder.sh -h' to see valid options"
    exit 1
    ;;
    esac
done


##check the assembly is given and that I can find it
[[ $assembly == "" && $assemblylist == "" ]] && echo "ERROR: No path to assembly/assemblylist was provided, assign one using -a or -al" && exit
[[ $assembly != "" && $assemblylist != "" ]] && echo "ERROR: Both an assembly AND list were provided, provide only one or the other using -a or -al" && exit

##warning if the telomeric repeat and/or number of threads was not modified
[[ $telomererepeat == "TTAGGG" ]] && echo "WARNING: Using default telomeric repeat sequence (TTAGGG)"
[[ $threads == "1" ]] && echo "WARNING: Using 1 thread for analysis"

##get a kb relative tipsize for output file tags
tipsize2=$( echo ${tipsize} | awk '{print $1/1000}' )

##begin run if only an assembly is provided
if [[ $assembly != ""  ]]
then

echo "An assembly was provided using -a; therefore running TLRFinder on a single input assembly"

assemblypath=$( realpath ${assembly} )
[ ! -f "${assemblypath}" ] && echo "ERROR: Cannot find path to assembly file provided by -a; check path is correct and file exists" && exit

##check if assembly is compressed or not
compressed=$( echo ${assembly} | awk -F "." '{if($NF == "gz") {print "yes"} else {print "no"}}' )

#######################################################################################################
##############################BEGIN FINDING THE TELOMERE-LINKED-REPEATS################################
#######################################################################################################


echo "######### Starting TLRFinder"

[ -d "${output}" ] && echo "ERROR: output folder already exists" && exit

echo "######### Output in ${output}"


##make an output directory and move into it
mkdir ${output}
cd ${output}

mkdir assemblies
mkdir contig_ends
mkdir contig_ends_alignments
mkdir telomeric_repeats
mkdir contig_ends_coverage
mkdir subtelomeric_repeats
mkdir plotting_Rscripts

##geneate header for the summary file 
echo "sample;assembly;contigs;telomeric_repeats;telomeric_repeats_contig_ends;TLR_regions;TLR_regions_contig_ends;contig_ends_TLR;TLR_representative_coords;TLR_representative_size;TLR_average_size;TLR_representative" | tr ';' '\t' > summary_stats.tsv


##create link to assembly in the output folder
ln -sf ${assemblypath} assemblies/
assembly=$( echo ${assembly} | awk -F "/" '{print "assemblies/"$NF}' )

##get the sametools index as a quick way to generate a file with contig sizes etc
samtools faidx ${assembly}

##getting the tips of contigs larger than 2 times the tip size
cat ${assembly}.fai  | awk -v tipsize="$tipsize" '{if($2 > (2*tipsize)) {print $1"\t"$2} }' | while read line
do
contig=$( echo "${line}" | awk -F "\t" '{print $1}' )
size=$( echo "${line}" | awk -F "\t" '{print $2}' )
if [ $size -gt ${tipsize} ]
then
end=$( grep ^${contig} ${assembly}.fai | awk -F "\t" -v contig="$contig" '{if($1 == contig) print $2}' | awk -v tipsize="$tipsize" '{if(($1-tipsize)> 0) {print ($1-tipsize)"-"$1} else {print "1-"$1}}' )
samtools faidx ${assembly} ${contig}:1-${tipsize}
samtools faidx ${assembly} ${contig}:${end}
fi
done > contig_ends/${prefix}.${tipsize2}kb_ends.fa

##can also label all regions with the canonical telomeric repeat
##use seqkit to locate the position of the canonical repeat then merge all those locations with a buffer of 7bp incase one repeat is off and export a bed file
echo "contig;start;end;sense" | tr ';' '\t' > telomeric_repeats/${prefix}.telomeres.bed
seqkit locate -j ${threads} --ignore-case -r -p "\"${telomererepeat}\""  ${assembly} | tail -n+2 | awk '{print $1"\t"$5"\t"$6"\t"$4}' | sort -k1,1 -k2,2n | bedtools merge -d 7 -c 4 -o distinct -i - | awk '{if($3-$2 > 11) print}' | awk -F "," '{print $1}' >> telomeric_repeats/${prefix}.telomeres.bed

echo "contig;start;end;sense" | tr ';' '\t' > telomeric_repeats/${prefix}.${tipsize2}kb_ends.telomeres.bed
seqkit locate -j ${threads} --ignore-case -r -p "\"${telomererepeat}\""  contig_ends/${prefix}.${tipsize2}kb_ends.fa | tail -n+2 | awk '{print $1"\t"$5"\t"$6"\t"$4}' | sort -k1,1 -k2,2n | bedtools merge -d 7 -c 4 -o distinct -i - | awk '{if($3-$2 > 11) print}' | awk -F "," '{print $1}' >> telomeric_repeats/${prefix}.${tipsize2}kb_ends.telomeres.bed

##subset the ends to only those with a good telomere sequence (>50bp)
#cat telomeric_repeats/${prefix}.${tipsize2}kb_ends.telomeres.bed  | awk '{if($3-$2 > 50) print $1}' > contig_ends/${prefix}.${tipsize2}kb_ends.with_tel.txt
#seqkit grep -f contig_ends/${prefix}.${tipsize2}kb_ends.with_tel.txt contig_ends/${prefix}.${tipsize2}kb_ends.fa > contig_ends/${prefix}.${tipsize2}kb_ends.with_tel.fa

##get simple bed file of just the ends to use for visualisation
samtools faidx contig_ends/${prefix}.${tipsize2}kb_ends.fa
echo "contig;start;end" | tr ';' '\t' > contig_ends/${prefix}.${tipsize2}kb_ends.bed
cat contig_ends/${prefix}.${tipsize2}kb_ends.fa.fai | awk -F "\t" '{print $1"\t1\t"$2}' >> contig_ends/${prefix}.${tipsize2}kb_ends.bed

##align the whole regions to themselves (once with all ends and another with just telomeric ends)
nucmer --maxmatch --threads ${threads} --delta contig_ends_alignments/${prefix}.${tipsize2}kb_ends.nucmer.delta contig_ends/${prefix}.${tipsize2}kb_ends.fa contig_ends/${prefix}.${tipsize2}kb_ends.fa
##convert to paf format and remove self matches and the reciprical alignments
paftools.js delta2paf contig_ends_alignments/${prefix}.${tipsize2}kb_ends.nucmer.delta | awk '{
    q = $1; r = $6;
    qs = $3; qe = $4;
    rs = $8; re = $9;

    if (q == r) next;  # Skip self-alignments (same query and reference)

    # Canonical key: always sort (query, reference) alphabetically
    if (q < r) {
        key = q "|" qs "|" qe "|" r "|" rs "|" re;
    } else {
        key = r "|" rs "|" re "|" q "|" qs "|" qe;
    }

    if (!seen[key]++) print;
}' > contig_ends_alignments/${prefix}.${tipsize2}kb_ends.nucmer.paf


##make a bedfile in order to use for coverage analysis of the contig ends alignments
cat contig_ends/${prefix}.${tipsize2}kb_ends.fa.fai | awk '{print $1"\t1\t"$2}' | bedtools makewindows -w ${window} -s ${slide} -b - > contig_ends/${prefix}.${tipsize2}kb_ends.${window}bpwindow.bed

##create a bed file from the alignments
cat contig_ends_alignments/${prefix}.${tipsize2}kb_ends.nucmer.paf  | awk '{if($3 > $4) {print $1"\t"$4"\t"$3} else {print $1"\t"$3"\t"$4}}' > contig_ends_alignments/${prefix}.${tipsize2}kb_ends.nucmer.paf.bed

##get all regions with at least XX coverage and larger than XXkb
##to get the minimum coverage we could estimate this by taking the number of contigs with +50bp telomeric sequences found and want at least 75% that coverage (unless 0 then make it 4)
##a 'good' telomere here is hard set to within 100bp of the end of the contigs (rough but needed for some filtering..)
goodtel=$(awk -F '\t' -v tipsize="$tipsize" \
    '($3-$2 > 50) && (($2 < 100 && $3 > 0) || ($2 < tipsize && $3 > tipsize-100)) {print $1}' \
    telomeric_repeats/${prefix}.${tipsize2}kb_ends.telomeres.bed | wc -l)
covmin=$( echo "${goodtel}" | awk -v covthreshold="$covthreshold" '{if($1 != "0") {print $1*covthreshold} else {print "4"}}' )
##size min set to 2kb as default 
bedtools coverage -a contig_ends/${prefix}.${tipsize2}kb_ends.${window}bpwindow.bed -b contig_ends_alignments/${prefix}.${tipsize2}kb_ends.nucmer.paf.bed > contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.cov.bed
echo "contig;start;end" | tr ';' '\t' > contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.bed
cat contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.cov.bed | awk -v covmin="$covmin" '{if($4 > covmin) print}' | bedtools merge -d 10 | awk -v sizemin="$sizemin"  '{if($3-$2 > sizemin) print}' >> contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.bed

repeatcount=$( wc -l contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.bed | awk '{print $1}' )
##create warning if no repeats were found
if [ ${repeatcount} -gt 2 ]
then


##we shall now have identified any repeats that a common across these contig ends
##now we want to extract a reference sequence for the repeat based on all the coordinates identifed
tail -n+2 contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.bed | bedtools getfasta -fi contig_ends/${prefix}.${tipsize2}kb_ends.fa -bed - -fo contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.fa

##then generate clusters based on a 80/80/80 threshold
##and selecting a representative per cluster (the longest sequence meeting the threshold)
cd-hit-est -i contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.fa -o contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.clustered.fa -d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500 > cd-hit-est.log
##now we have the representative and clusters, we want to select the cluster with the great number of alignments (removing some smaller regions)
##then we can BLAST for our representative across the regions again and get a better indicator of its appearance.
cluster=$( awk '/Cluster/{c++} />/{print "Cluster"c-1}' contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.clustered.fa.clstr | sort | uniq -c | sort -nr  | head -n1 | awk '{print $2}' )
awk '/Cluster/{c++} />/{print "Cluster"c-1"\t"$0}' contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.clustered.fa.clstr | grep "^${cluster}" | awk '{print $4}' | sed 's/\.\.\.//' | sed 's/>//' > subtelomeric_repeats/${prefix}.repeat_cluster.txt
rep=$( awk '/Cluster/{c++} />/{print "Cluster"c-1"\t"$0}' contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.clustered.fa.clstr | grep "^${cluster}" | grep "\*$" | awk '{print $4}' | sed 's/\.\.\.//' | sed 's/>//' )
##extract just the rep
seqkit grep --quiet -p "${rep}" contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.clustered.fa > subtelomeric_repeats/${prefix}.repeat_rep.fa
seqkit grep --quiet -f subtelomeric_repeats/${prefix}.repeat_cluster.txt contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.fa > subtelomeric_repeats/${prefix}.repeat_cluster.fa

##now use the representative to BLAST against the ends again and create a new bed file for the coords
##just need some filters...
##an alignment of at least 80% identity and 500bb
blastn  -subject contig_ends/${prefix}.${tipsize2}kb_ends.fa -query subtelomeric_repeats/${prefix}.repeat_rep.fa -outfmt 6 > subtelomeric_repeats/${prefix}.repeat_rep.ends_blast.tsv
##create bed from nonredundant positions
echo "contig;start;end" | tr ';' '\t' > subtelomeric_repeats/${prefix}.repeat_rep.ends_blast.bed
cat subtelomeric_repeats/${prefix}.repeat_rep.ends_blast.tsv | awk '{if($3 > 80 && $4 > 500) print}' | awk '{print $2"\t"$9"\t"$10}' | awk '{if($2>$3) {print $1"\t"$3"\t"$2} else {print}}' | sort -k1,1 -k2,2n | bedtools merge -d 10  >> subtelomeric_repeats/${prefix}.repeat_rep.ends_blast.bed

##and the whole genome
if [[ $compressed == "yes" ]]
then
gzip -dc ${assembly} | blastn -subject - -query subtelomeric_repeats/${prefix}.repeat_rep.fa -outfmt 6 > subtelomeric_repeats/${prefix}.repeat_rep.WG_blast.tsv
else
blastn -subject ${assembly} -query subtelomeric_repeats/${prefix}.repeat_rep.fa -outfmt 6 > subtelomeric_repeats/${prefix}.repeat_rep.WG_blast.tsv
fi
echo "contig;start;end" | tr ';' '\t' > subtelomeric_repeats/${prefix}.repeat_rep.WG_blast.bed
cat subtelomeric_repeats/${prefix}.repeat_rep.WG_blast.tsv | awk '{if($3 > 80 && $4 > 500) print}' | awk '{print $2"\t"$9"\t"$10}' | awk '{if($2>$3) {print $1"\t"$3"\t"$2} else {print}}' | sort -k1,1 -k2,2n | bedtools merge -d 10 >> subtelomeric_repeats/${prefix}.repeat_rep.WG_blast.bed
##get bed file for the whole genome so we can plot the positions of the repeats alongside it
echo "contig;start;end" | tr ';' '\t' > assemblies/${prefix}.genome.bed
cat ${assembly}.fai | awk '{print $1"\t1\t"$2}'  >> assemblies/${prefix}.genome.bed

echo "######### Plotting contig end alignments and genome-wide distribution"


Rscriptpath=$( which Finderplots_TLRs.R )
cat ${Rscriptpath} | sed "s/SAMPLE/${prefix}/g" | sed "s/TIPSIZE2/${tipsize2}/g" | sed "s/TIPSIZE/${tipsize}/g" > plotting_Rscripts/${prefix}.R

Rscript plotting_Rscripts/${prefix}.R

else
echo "######### WARNING: No TLRs were found in ${prefix}:${assembly}"
fi

echo "######### Adding results to summary file"

##sample = $prefix
assembly2=$( echo ${assembly}  | awk -F "/" '{print $NF}' )
assemblytype=$( echo ${assembly2} | awk -F "." '{if($NF == "gz" || $NF == "bgzip" || $NF == "gzip" ) {print "compressed"} else {print "uncompressed"}}' )
contigs=$( if [[ $assemblytype == "uncompressed"  ]]; then grep '>' ${assemblypath} | wc -l ; else zgrep '>' ${assemblypath} | wc -l ; fi )
telomericrepeats=$( if [ -e "telomeric_repeats/${prefix}.${tipsize2}kb_ends.telomeres.bed" ]; then cat telomeric_repeats/${prefix}.${tipsize2}kb_ends.telomeres.bed  |  awk '{if($3-$2 > 50) print $1}'  | wc -l ; else echo "0" ; fi )
telomericends=$( if [ -e "telomeric_repeats/${prefix}.${tipsize2}kb_ends.telomeres.bed" ]; then cat telomeric_repeats/${prefix}.${tipsize2}kb_ends.telomeres.bed  |  awk '{if($3-$2 > 50) print $1}'  | sort -u | wc -l ; else echo "0" ; fi  )
TLRs=$( if [ -e "subtelomeric_repeats/${prefix}.repeat_rep.WG_blast.bed" ]; then cat subtelomeric_repeats/${prefix}.repeat_rep.WG_blast.bed | wc -l ; else echo "0" ; fi  )
TLRends=$( if [ -e "subtelomeric_repeats/${prefix}.repeat_rep.ends_blast.bed" ]; then cat subtelomeric_repeats/${prefix}.repeat_rep.ends_blast.bed | wc -l ; else echo "0" ; fi  )
endsTLR=$( if [ -e "subtelomeric_repeats/${prefix}.repeat_rep.ends_blast.bed" ]; then cat subtelomeric_repeats/${prefix}.repeat_rep.ends_blast.bed | awk -F "\t" '{print $1}' | sort -u | wc -l ; else echo "0" ; fi  )
repeatrepcoords=$( if [ -e "subtelomeric_repeats/${prefix}.repeat_rep.fa" ]; then grep '>' subtelomeric_repeats/${prefix}.repeat_rep.fa | awk '{if($0 ~ ">") {print $0} else {print "NA"}}' | sed 's/>//g' | tr '-' '\t' | tr ':' '\t' | awk '{if($0 == "NA") {print $0} else if($4 == "") { print  $1":"$2"-"$3} else if($2 == "1") {print $1":"$4"-"$5} else {print $1":"($2+$4)"-"(($2+$4)+($5-$4))}}' ; else echo "NA" ; fi  )
repeatrepsize=$(if [ -e "subtelomeric_repeats/${prefix}.repeat_rep.fa" ]; then grep -v '>' subtelomeric_repeats/${prefix}.repeat_rep.fa | tr '\n' 'X' | sed 's/X//g' | wc -c ; else echo "NA" ; fi )
repeatavgsize=$( if [ -e "subtelomeric_repeats/${prefix}.repeat_cluster.fa" ]; then seqkit stat -T subtelomeric_repeats/${prefix}.repeat_cluster.fa  | awk 'NR>1 {print $7}' | awk -F "." '{print $1}' ; else echo "NA" ; fi )
repeatrep=$( if [ -e "subtelomeric_repeats/${prefix}.repeat_rep.fa" ]; then grep -v '>' subtelomeric_repeats/${prefix}.repeat_rep.fa | tr '\n' 'X' | sed 's/X//g' ; else echo "NA" ; fi  )

echo "${prefix};${assembly2};${contigs};${telomericrepeats};${telomericends};${TLRs};${TLRends};${endsTLR};${repeatrepcoords};${repeatrepsize};${repeatavgsize};${repeatrep}" | tr ';' '\t' >> summary_stats.tsv


else
##RUNNING NOW IF A LIST OF ASSEMBLIES WERE PROVIDED

echo "A list of assemblies was provided using -al; therefore running TLRFinder on a set of assemblies"


assemblylistpath=$( realpath ${assemblylist} )
[ ! -f "${assemblylistpath}" ] && echo "ERROR: Cannot find path to assembly list provided by -al; check path is correct and file exists" && exit

#######################################################################################################
###############################BEGIN FINDING THE TELOMERE-LINKED-REPEATS############################### 
#######################################################################################################


echo "######### Starting TLRFinder"
echo "######### Output in ${output}"


##make an output directory and move into it
mkdir ${output}

mkdir ${output}/assemblies
mkdir ${output}/contig_ends
mkdir ${output}/contig_ends_alignments
mkdir ${output}/telomeric_repeats
mkdir ${output}/contig_ends_coverage
mkdir ${output}/subtelomeric_repeats
mkdir ${output}/plotting_Rscripts

##geneate header for the summary file 
echo "sample;assembly;contigs;telomeric_repeats;telomeric_repeats_contig_ends;TLR_regions;TLR_regions_contig_ends;contig_ends_TLR;TLR_representative_coords;TLR_representative_size;TLR_average_size;TLR_representative" | tr ';' '\t' > ${output}/summary_stats.tsv


cat ${assemblylistpath} | while read line
do
prefix=$( echo "${line}" | awk -F "\t" '{print $1}' )
assembly=$( echo "${line}" | awk -F "\t" '{print $2}' )
assemblypath=$( realpath ${assembly} )
[ ! -f "${assemblypath}" ] && echo "ERROR: Cannot find path to assembly from sample ${prefix} provided by -al; check path to assembly is correct and file exists" && exit

cd ${output}

echo "######### Running TLRFinder on sample ${prefix}"

##check if assembly is compressed or not
compressed=$( echo ${assembly} | awk -F "." '{if($NF == "gz") {print "yes"} else {print "no"}}' )

##create link to assembly in the output folder
if [[ $compressed == "yes" ]]
then
ln -sf ${assemblypath} assemblies/${prefix}.fa.gz
assembly=$( echo "assemblies/${prefix}.fa.gz" )
else 
ln -sf ${assemblypath} assemblies/${prefix}.fa
assembly=$( echo "assemblies/${prefix}.fa" )
fi

##get the sametools index as a quick way to generate a file with contig sizes etc
samtools faidx ${assembly}

##getting the tips of contigs larger than 2 times the tip size
cat ${assembly}.fai  | awk -v tipsize="$tipsize" '{if($2 > (2*tipsize)) {print $1"\t"$2} }' | while read line
do
contig=$( echo "${line}" | awk -F "\t" '{print $1}' )
size=$( echo "${line}" | awk -F "\t" '{print $2}' )
if [ $size -gt ${tipsize} ]
then
end=$( grep ^${contig} ${assembly}.fai | awk -F "\t" -v contig="$contig" '{if($1 == contig) print $2}' | awk -v tipsize="$tipsize" '{if(($1-tipsize)> 0) {print ($1-tipsize)"-"$1} else {print "1-"$1}}' )
samtools faidx ${assembly} ${contig}:1-${tipsize}
samtools faidx ${assembly} ${contig}:${end}
fi
done > contig_ends/${prefix}.${tipsize2}kb_ends.fa

##can also label all regions with the canonical telomeric repeat
##use seqkit to locate the position of the canonical repeat then merge all those locations with a buffer of 7bp incase one repeat is off and export a bed file
echo "contig;start;end;sense" | tr ';' '\t' > telomeric_repeats/${prefix}.telomeres.bed
seqkit locate -j ${threads} --ignore-case -r -p "\"${telomererepeat}\""  ${assembly} | tail -n+2 | awk '{print $1"\t"$5"\t"$6"\t"$4}' | sort -k1,1 -k2,2n | bedtools merge -d 7 -c 4 -o distinct -i - | awk '{if($3-$2 > 11) print}' | awk -F "," '{print $1}' >> telomeric_repeats/${prefix}.telomeres.bed

echo "contig;start;end;sense" | tr ';' '\t' > telomeric_repeats/${prefix}.${tipsize2}kb_ends.telomeres.bed
seqkit locate -j ${threads} --ignore-case -r -p "\"${telomererepeat}\""  contig_ends/${prefix}.${tipsize2}kb_ends.fa | tail -n+2 | awk '{print $1"\t"$5"\t"$6"\t"$4}' | sort -k1,1 -k2,2n | bedtools merge -d 7 -c 4 -o distinct -i - | awk '{if($3-$2 > 11) print}' | awk -F "," '{print $1}' >> telomeric_repeats/${prefix}.${tipsize2}kb_ends.telomeres.bed

##subset the ends to only those with a good telomere sequence (>50bp)
#cat telomeric_repeats/${prefix}.${tipsize2}kb_ends.telomeres.bed  | awk '{if($3-$2 > 50) print $1}' > contig_ends/${prefix}.${tipsize2}kb_ends.with_tel.txt
#seqkit grep -f contig_ends/${prefix}.${tipsize2}kb_ends.with_tel.txt contig_ends/${prefix}.${tipsize2}kb_ends.fa > contig_ends/${prefix}.${tipsize2}kb_ends.with_tel.fa

##get simple bed file of just the ends to use for visualisation
samtools faidx contig_ends/${prefix}.${tipsize2}kb_ends.fa
echo "contig;start;end" | tr ';' '\t' > contig_ends/${prefix}.${tipsize2}kb_ends.bed
cat contig_ends/${prefix}.${tipsize2}kb_ends.fa.fai | awk -F "\t" '{print $1"\t1\t"$2}' >> contig_ends/${prefix}.${tipsize2}kb_ends.bed

##align the whole regions to themselves (once with all ends and another with just telomeric ends)
nucmer --maxmatch --threads ${threads} --delta contig_ends_alignments/${prefix}.${tipsize2}kb_ends.nucmer.delta contig_ends/${prefix}.${tipsize2}kb_ends.fa contig_ends/${prefix}.${tipsize2}kb_ends.fa
##convert to paf format and remove self matches and the reciprical alignments
paftools.js delta2paf contig_ends_alignments/${prefix}.${tipsize2}kb_ends.nucmer.delta | awk '{
    q = $1; r = $6;
    qs = $3; qe = $4;
    rs = $8; re = $9;

    if (q == r) next;  # Skip self-alignments (same query and reference)

    # Canonical key: always sort (query, reference) alphabetically
    if (q < r) {
        key = q "|" qs "|" qe "|" r "|" rs "|" re;
    } else {
        key = r "|" rs "|" re "|" q "|" qs "|" qe;
    }

    if (!seen[key]++) print;
}' > contig_ends_alignments/${prefix}.${tipsize2}kb_ends.nucmer.paf


##make a bedfile in order to use for coverage analysis of the contig ends alignments
cat contig_ends/${prefix}.${tipsize2}kb_ends.fa.fai | awk '{print $1"\t1\t"$2}' | bedtools makewindows -w ${window} -s ${slide} -b - > contig_ends/${prefix}.${tipsize2}kb_ends.${window}bpwindow.bed

##create a bed file from the alignments
cat contig_ends_alignments/${prefix}.${tipsize2}kb_ends.nucmer.paf  | awk '{if($3 > $4) {print $1"\t"$4"\t"$3} else {print $1"\t"$3"\t"$4}}' > contig_ends_alignments/${prefix}.${tipsize2}kb_ends.nucmer.paf.bed

##get all regions with at least XX coverage and larger than XXkb
##to get the minimum coverage we could estimate this by taking the number of contigs with +50bp telomeric sequences found and want at least 75% that coverage (unless 0 then make it 4)
##a 'good' telomere here is hard set to within 100bp of the end of the contigs (rough but needed for some filtering..)
goodtel=$(awk -F '\t' -v tipsize="$tipsize" \
    '($3-$2 > 50) && (($2 < 100 && $3 > 0) || ($2 < tipsize && $3 > tipsize-100)) {print $1}' \
    telomeric_repeats/${prefix}.${tipsize2}kb_ends.telomeres.bed | wc -l)
covmin=$( echo "${goodtel}" | awk -v covthreshold="$covthreshold" '{if($1 != "0") {print $1*covthreshold} else {print "4"}}' )
##size min set to 2kb as default 
bedtools coverage -a contig_ends/${prefix}.${tipsize2}kb_ends.${window}bpwindow.bed -b contig_ends_alignments/${prefix}.${tipsize2}kb_ends.nucmer.paf.bed > contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.cov.bed
echo "contig;start;end" | tr ';' '\t' > contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.bed
cat contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.cov.bed | awk -v covmin="$covmin" '{if($4 > covmin) print}' | bedtools merge -d 10 | awk -v sizemin="$sizemin"  '{if($3-$2 > sizemin) print}' >> contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.bed

repeatcount=$( wc -l contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.bed | awk '{print $1}' )
##create warning if no repeats were found
if [ ${repeatcount} -gt 2 ]
then

##we shall now have identified any repeats that a common across these contig ends
##now we want to extract a reference sequence for the repeat based on all the coordinates identifed
tail -n+2 contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.bed | bedtools getfasta -fi contig_ends/${prefix}.${tipsize2}kb_ends.fa -bed - -fo contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.fa

##then generate clusters based on a 80/80/80 threshold
##and selecting a representative per cluster (the longest sequence meeting the threshold)
cd-hit-est -i contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.fa -o contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.clustered.fa -d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500 > cd-hit-est.log
##now we have the representative and clusters, we want to select the cluster with the great number of alignments (removing some smaller regions)
##then we can BLAST for our representative across the regions again and get a better indicator of its appearance.
cluster=$( awk '/Cluster/{c++} />/{print "Cluster"c-1}' contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.clustered.fa.clstr | sort | uniq -c | sort -nr  | head -n1 | awk '{print $2}' )
awk '/Cluster/{c++} />/{print "Cluster"c-1"\t"$0}' contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.clustered.fa.clstr | grep "^${cluster}" | awk '{print $4}' | sed 's/\.\.\.//' | sed 's/>//' > subtelomeric_repeats/${prefix}.repeat_cluster.txt
rep=$( awk '/Cluster/{c++} />/{print "Cluster"c-1"\t"$0}' contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.clustered.fa.clstr | grep "^${cluster}" | grep "\*$" | awk '{print $4}' | sed 's/\.\.\.//' | sed 's/>//' )
##extract just the rep
seqkit grep --quiet -p "${rep}" contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.clustered.fa > subtelomeric_repeats/${prefix}.repeat_rep.fa
seqkit grep --quiet -f subtelomeric_repeats/${prefix}.repeat_cluster.txt contig_ends_coverage/${prefix}.${tipsize2}kb_ends.nucmer.paf.repeats.fa > subtelomeric_repeats/${prefix}.repeat_cluster.fa

##now use the representative to BLAST against the ends again and create a new bed file for the coords
##just need some filters...
##an alignment of at least 80% identity and 500bb
blastn  -subject contig_ends/${prefix}.${tipsize2}kb_ends.fa -query subtelomeric_repeats/${prefix}.repeat_rep.fa -outfmt 6 > subtelomeric_repeats/${prefix}.repeat_rep.ends_blast.tsv
##create bed from nonredundant positions
echo "contig;start;end" | tr ';' '\t' > subtelomeric_repeats/${prefix}.repeat_rep.ends_blast.bed
cat subtelomeric_repeats/${prefix}.repeat_rep.ends_blast.tsv | awk '{if($3 > 80 && $4 > 500) print}' | awk '{print $2"\t"$9"\t"$10}' | awk '{if($2>$3) {print $1"\t"$3"\t"$2} else {print}}' | sort -k1,1 -k2,2n | bedtools merge -d 10  >> subtelomeric_repeats/${prefix}.repeat_rep.ends_blast.bed

##and the whole genome
if [[ $compressed == "yes" ]]
then
gzip -dc ${assembly} | blastn -subject - -query subtelomeric_repeats/${prefix}.repeat_rep.fa -outfmt 6 > subtelomeric_repeats/${prefix}.repeat_rep.WG_blast.tsv
else
blastn -subject ${assembly} -query subtelomeric_repeats/${prefix}.repeat_rep.fa -outfmt 6 > subtelomeric_repeats/${prefix}.repeat_rep.WG_blast.tsv
fi
echo "contig;start;end" | tr ';' '\t' > subtelomeric_repeats/${prefix}.repeat_rep.WG_blast.bed
cat subtelomeric_repeats/${prefix}.repeat_rep.WG_blast.tsv | awk '{if($3 > 80 && $4 > 500) print}' | awk '{print $2"\t"$9"\t"$10}' | awk '{if($2>$3) {print $1"\t"$3"\t"$2} else {print}}' | sort -k1,1 -k2,2n | bedtools merge -d 10 >> subtelomeric_repeats/${prefix}.repeat_rep.WG_blast.bed
##get bed file for the whole genome so we can plot the positions of the repeats alongside it
echo "contig;start;end" | tr ';' '\t' > assemblies/${prefix}.genome.bed
cat ${assembly}.fai | awk '{print $1"\t1\t"$2}'  >> assemblies/${prefix}.genome.bed

echo "######### Plotting contig end alignments and genome-wide distribution"


Rscriptpath=$( which Finderplots_TLRs.R )
cat ${Rscriptpath} | sed "s/SAMPLE/${prefix}/g" | sed "s/TIPSIZE2/${tipsize2}/g" | sed "s/TIPSIZE/${tipsize}/g" > plotting_Rscripts/${prefix}.R

Rscript plotting_Rscripts/${prefix}.R


else
echo "######### WARNING: No TLRs were found in ${prefix}:${assembly}"
fi


echo "######### Adding results to summary file"

##sample = $prefix
assembly2=$( cat ${assemblylistpath}  | awk -v prefix="$prefix" -F "\t" '{if($1== prefix) print $2}' )
assemblytype=$( echo ${assembly2} | awk -F "." '{if($NF == "gz" || $NF == "bgzip" || $NF == "gzip" ) {print "compressed"} else {print "uncompressed"}}' )
contigs=$( if [[ $assemblytype == "uncompressed"  ]]; then grep '>' ${assemblypath} | wc -l ; else zgrep '>' ${assemblypath} | wc -l ; fi )
telomericrepeats=$( if [ -e "telomeric_repeats/${prefix}.${tipsize2}kb_ends.telomeres.bed" ]; then cat telomeric_repeats/${prefix}.${tipsize2}kb_ends.telomeres.bed  |  awk '{if($3-$2 > 50) print $1}'  | wc -l ; else echo "0" ; fi )
telomericends=$( if [ -e "telomeric_repeats/${prefix}.${tipsize2}kb_ends.telomeres.bed" ]; then cat telomeric_repeats/${prefix}.${tipsize2}kb_ends.telomeres.bed  |  awk '{if($3-$2 > 50) print $1}'  | sort -u | wc -l ; else echo "0" ; fi  )
TLRs=$( if [ -e "subtelomeric_repeats/${prefix}.repeat_rep.WG_blast.bed" ]; then cat subtelomeric_repeats/${prefix}.repeat_rep.WG_blast.bed | wc -l ; else echo "0" ; fi  )
TLRends=$( if [ -e "subtelomeric_repeats/${prefix}.repeat_rep.ends_blast.bed" ]; then cat subtelomeric_repeats/${prefix}.repeat_rep.ends_blast.bed | wc -l ; else echo "0" ; fi  )
endsTLR=$( if [ -e "subtelomeric_repeats/${prefix}.repeat_rep.ends_blast.bed" ]; then cat subtelomeric_repeats/${prefix}.repeat_rep.ends_blast.bed | awk -F "\t" '{print $1}' | sort -u | wc -l ; else echo "0" ; fi  )
repeatrepcoords=$( if [ -e "subtelomeric_repeats/${prefix}.repeat_rep.fa" ]; then grep '>' subtelomeric_repeats/${prefix}.repeat_rep.fa | awk '{if($0 ~ ">") {print $0} else {print "NA"}}' | sed 's/>//g' | tr '-' '\t' | tr ':' '\t' | awk '{if($0 == "NA") {print $0} else if($4 == "") { print  $1":"$2"-"$3} else if($2 == "1") {print $1":"$4"-"$5} else {print $1":"($2+$4)"-"(($2+$4)+($5-$4))}}' ; else echo "NA" ; fi  )
repeatrepsize=$(if [ -e "subtelomeric_repeats/${prefix}.repeat_rep.fa" ]; then grep -v '>' subtelomeric_repeats/${prefix}.repeat_rep.fa | tr '\n' 'X' | sed 's/X//g' | wc -c ; else echo "NA" ; fi )
repeatavgsize=$( if [ -e "subtelomeric_repeats/${prefix}.repeat_cluster.fa" ]; then seqkit stat -T subtelomeric_repeats/${prefix}.repeat_cluster.fa  | awk 'NR>1 {print $7}' | awk -F "." '{print $1}' ; else echo "NA" ; fi )
repeatrep=$( if [ -e "subtelomeric_repeats/${prefix}.repeat_rep.fa" ]; then grep -v '>' subtelomeric_repeats/${prefix}.repeat_rep.fa | tr '\n' 'X' | sed 's/X//g' ; else echo "NA" ; fi  )

echo "${prefix};${assembly2};${contigs};${telomericrepeats};${telomericends};${TLRs};${TLRends};${endsTLR};${repeatrepcoords};${repeatrepsize};${repeatavgsize};${repeatrep}" | tr ';' '\t' >> summary_stats.tsv




cd ../

done


cd ${output}

echo "######### Running comparisons on TLRs found across the set of assemblies provided"

###because several assemblies were provided we can now compare their TLRs


echo "######### Generating a k-mer based NJ tree using mashtree"

##generate a mashtree of the assemblies
##first just check if the assemblies are all compressed/uncompressed or a mix
assemblytype=$( cat ${assemblylistpath} | awk -F "\t" 'BEGIN{compressed="no"; uncompressed="no"; mixed="no"} {if($2 ~ ".gz" && uncompressed == "no") {compressed="yes"} else if($2 ~ ".gz" && uncompressed == "yes") {compressed="yes"; mixed="yes"} else if($2 !~ ".gz" && compressed == "no") {uncompressed="yes"} else if($2 !~ ".gz" && compressed == "yes") {uncompressed="yes"; mixed="yes"}} END{if(mixed=="yes") {print "mixed"} else if(compressed=="yes") {print "compressed"} else {print "uncompressed"}}' )

if [[ $assemblytype == "mixed"  ]]
then
mashtree_bootstrap.pl --reps ${bootstraps} --numcpus ${threads} assemblies/*.fa assemblies/*.fa.gz -- --mindepth 0 --sort-order random > assemblies.mashtree.bootstrap.dnd 2> mashtree.log
else
if [[ $assemblytype == "compressed"  ]]
then
mashtree_bootstrap.pl --reps ${bootstraps} --numcpus ${threads} assemblies/*.fa.gz -- --mindepth 0 --sort-order random > assemblies.mashtree.bootstrap.dnd 2> mashtree.log
else
if [[ $assemblytype == "uncompressed"  ]]
then
mashtree_bootstrap.pl --reps ${bootstraps}  --numcpus ${threads} assemblies/*.fa -- --mindepth 0 --sort-order random > assemblies.mashtree.bootstrap.dnd 2> mashtree.log
fi
fi
fi


echo "######### Comparing global-ANI stats between TLRs within and between assemblies"

##generate lz-ani similarity comparisons
##within an assembly
mkdir subtelomeric_repeats_comparisons/
echo "sample;repeat_representative;average_gANI;count" | tr ';' '\t' > gANI.within_repeats.tsv
ls subtelomeric_repeats/*.WG_blast.bed | while read file
do
sample=$( echo "${file}" | awk -F "/" '{print $NF}' |  sed 's/.repeat_rep.WG_blast.bed//g' )
assembly=$( ls assemblies/ | grep "${sample}\.fa" | grep -v "\.fai"$ | grep -v "\.bed"$ | grep -v "\.gzi"$ )
replength=$( grep -v '>' subtelomeric_repeats/${sample}.repeat_rep.fa | tr '\n' 'X' | sed 's/X//g' | wc -c  )
rep=$( grep '>' subtelomeric_repeats/${sample}.repeat_rep.fa | sed 's/>//g' | tr '-' '\t' | tr ':' '\t' | awk -v tipsize="$tipsize" '{if($4 == "") { print  $1":"$2"-"$3} else if($2 == "1") {print $1":"$4"-"$5} else {print $1":"($2+$4)"-"(($2+$4)+($5-$4))}}'   )
tail -n+2 $file | awk -v replength="$replength" '{if($3-$2 > (0.25*replength)) {print}}' > subtelomeric_repeats_comparisons/${sample}.repeat_rep.WG_blast.bed
bedtools getfasta -fi assemblies/${assembly} -bed subtelomeric_repeats_comparisons/${sample}.repeat_rep.WG_blast.bed -fo subtelomeric_repeats_comparisons/${sample}.repeat_rep.WG_blast.fa
##number of repeats used
count=$( grep '>' subtelomeric_repeats_comparisons/${sample}.repeat_rep.WG_blast.fa | wc -l )
lz-ani all2all --in-fasta subtelomeric_repeats_comparisons/${sample}.repeat_rep.WG_blast.fa --out subtelomeric_repeats_comparisons/${sample}.repeat_rep.WG_blast.ani.tsv 2> lz-ani.${sample}.log

##taking the global ANI = "The number of identical bases across local alignments divided by the length of the query/reference genome"
##otherwise small alignments get found to have good ANI but could cover just a few bases
tail -n+2 subtelomeric_repeats_comparisons/${sample}.repeat_rep.WG_blast.ani.tsv | awk '{print $6}' | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }' | awk -v sample="$sample" -v rep="$rep" -v count="$count" '{print sample"\t"rep"\t"$1"\t"count}' >> gANI.within_repeats.tsv

done

##between assemblies using the representative

##make a useful bed file of these representatives
echo "contig;start;end;sample" | tr ';' '\t' > repeat_representatives.bed

##get a refined list of the representatives per repeat
tail -n+2 gANI.within_repeats.tsv | awk -F "\t" '{print $1"\t"$2}' | while read bed
do
sample=$( echo "${bed}" | awk -F "\t" '{print $1}' )
coords=$( echo "${bed}" | awk -F "\t" '{print $2}' )
assembly=$( ls assemblies/ | grep "${sample}\.fa" | grep -v "\.fai"$ | grep -v "\.bed"$ | grep -v "\.gzi"$ )
echo "${coords};${sample}" | tr ';' '\t' | tr '-' '\t' | tr ':' '\t' >> repeat_representatives.bed
samtools faidx assemblies/${assembly} "${coords}"
done > repeat_representatives.fa

lz-ani all2all --in-fasta repeat_representatives.fa --out subtelomeric_repeats_comparisons/repeat_representatives.ani.tsv 2> lz-ani.repeat_representatives.log
echo "ref_sample;ref_rep;query_sample;query_rep;gANI" | tr ';' '\t' > gANI.between_repeat_representatives.tsv 
#tail -n+2 subtelomeric_repeats_comparisons/repeat_representatives.ani.tsv | awk -F "\t" '{print $3"\t"$4"\t"$6}' | awk -F "_" '{print $1"\t"$2"\t"$3}' | awk -F "\t" '{print $1"\t"$1"_"$2"\t"$3"\t"$3"_"$4"\t"$5}' >> gANI.between_repeat_representatives.tsv 

tail -n+2 subtelomeric_repeats_comparisons/repeat_representatives.ani.tsv | awk -F "\t" '{print $3"\t"$4"\t"$6}' | while read line
do
refrep=$( echo "${line}"  | awk -F "\t" '{print $1}' )
queryrep=$( echo "${line}"  | awk -F "\t" '{print $2}' )
gANI=$( echo "${line}"  | awk -F "\t" '{print $3}' )
refsample=$( cat repeat_representatives.bed | awk -F "\t" '{print $1":"$2"-"$3"\t"$4}' | awk -F "\t" -v refrep="$refrep" '{if($1 == refrep) {print $2}}' )
querysample=$( cat repeat_representatives.bed | awk -F "\t" '{print $1":"$2"-"$3"\t"$4}' | awk -F "\t" -v queryrep="$queryrep" '{if($1 == queryrep) {print $2}}' )
echo "${refsample};${refrep};${querysample};${queryrep};${gANI}" | tr ';' '\t' >> gANI.between_repeat_representatives.tsv 
done

tail -n+2 subtelomeric_repeats_comparisons/repeat_representatives.ani.tsv | awk -F "\t" '{print $3"\n"$4}' | sort -u | while read rep
do
repsample=$( cat repeat_representatives.bed | awk -F "\t" '{print $1":"$2"-"$3"\t"$4}' | awk -F "\t" -v rep="$rep" '{if($1 == rep) {print $2}}' )
echo "${repsample};${rep};${repsample};${rep};1" | tr ';' '\t'  >> gANI.between_repeat_representatives.tsv
done 


##generate the tree-heatmap plot

Rscriptpath2=$( which Compreplots_TLRs.R )
cat ${Rscriptpath2} > phylogeny_plus_gANI_heatmap.R

Rscript phylogeny_plus_gANI_heatmap.R


fi

rm Rplots.pdf


echo "######### TLRFinder has finished; E noho rā"
