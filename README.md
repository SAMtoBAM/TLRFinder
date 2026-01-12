<p align="center" >
    <img src="https://github.com/SAMtoBAM/TLRFinder/blob/main/logo/TLRFinder_logo.png" width=100%>
</p>

[![Zenodo DOI](https://zenodo.org/badge/1052979411.svg)](https://doi.org/10.5281/zenodo.17093898)
[![Anaconda_version](https://anaconda.org/samtobam/TLRFinder/badges/version.svg)](https://anaconda.org/samtobam/TLRFinder)
[![Anaconda_platforms](https://anaconda.org/samtobam/TLRFinder/badges/platforms.svg)](https://anaconda.org/samtobam/TLRFinder)
[![Anaconda_downloads](https://anaconda.org/samtobam/TLRFinder/badges/downloads.svg)](https://anaconda.org/samtobam/TLRFinder)
[![Anaconda-Server Badge](https://anaconda.org/samtobam/TLRFinder/badges/latest_release_date.svg)](https://anaconda.org/samtobam/TLRFinder)


# TLRFinder
TLRFinder is a tools desgined to detect and analyse **T**elomere-**L**inked-**R**epeats (TLR)

TLRs are, as the name describes, regions adjacent to telomeres that are conserved across chromosome ends <br/>
The most commonly described version contains Helicases, such as the Y' element in _S. cerevisiae_, and Telomere-Linked-Helicases (TLHs) <br/>
However other coding sequences also exist within TLRs and are also conserved across diverse fungi <br/>
With the advent of Long-read sequencing we now have an increasing number of well assembled genomes with the subtelomeres intact <br/>
Therefore we are now primed for looking at the TLRs and their evolution.

All that is required to run TLRFinder is an assembly in fasta format (can be bgzip compressed) or a tsv file containing a list of assemblies (first column: sample name; second column: assembly) <br/>
NOTE: Do not use hyphens '-' or other special characters in the sample/file names

This tool was developed for [O'Donnell et al. 2025]() (please cite this publication if you find this tool useful)

# Conda installation

    conda install samtobam::tlrfinder

# How to use

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
    -sm | --sizemin     The minimum size of a region passing the coverage threshold to be considered as a potential TLR region (Default: 2000)
    
    Multiple assembly specific parameters (if using --al)
    -b | --bootstraps   Number of bootstrap tests to be performed by mashtree (Default: 1000)

    Optional parameters:
    -w | --window       Number of basepairs for window averaging for coverage (Default: 10)
    -s | --slide        Number of basepairs for the window to slide for coverage (Default: 5)
    -p | --prefix       Prefix for output (Default: TLRFinder)
    -o | --output       Name of output folder for all results (default: TLRFinder_output)
    -h | --help         Print this help message



## How does TLRFinder work

TLRFinder works in 10 main steps (each step is run on all assemblies provided -al)

1. Extract 50kb from the ends of contigs >100kb
2. Identify Telomeric sequences genomes wide and determine the number within the contig ends
3. Align all contig ends to one another
4. Calculate the average coverage of sliding windows within the 50kb ends
5. Extract contiguous regions with a coverage > (0.75*'the number of telomeric sequences')
6. Remove regions < 2kb
7. Cluster the nucleotide sequence of the remaining regions using the 80/80 principal
8. Consider the largest cluster as the TLR and take the largest version as the representative
9. Search the whole genome using BLASTn using the representative TLR
10. Plot the alignment and whole genome positions of TLRs for manual verification/scrutiny <br/>

TLRFinder runs an additional 3 steps if provided multiple assemblies using -al <br/>
&nbsp; 11. Rapidly generate a k-mer NJ phylogeny using mashtree <br/>
&nbsp; 12. Calculate global-ANI (g-ANI) statistics for TLRs within an assembly and between TLR representatives <br/>
&nbsp; 13. Plot the tree and g-ANI stats side by side <br/>


## Important output files
-The primary output is the **summary_stats.tsv** file with some summarising data per assembly <br/>
    Column 1 : 'sample' : name given to samples in -al file or prefix if a single assembly (-a) <br/>
    Column 2 : 'assembly' : name of assembly file provided <br/>
    Column 3 : 'contigs' : Number of contigs in assembly <br/>
    Column 4 : 'telomeric_repeats' : Number of good telomeric repeat regions found (minimum of 50bp) <br/>
    Column 5 : 'telomeric_repeats_contig_ends' : Number of the above telomeric repeats found within the contig ends (Default: 50kb from end of contigs larger than 100kb) <br/>
    Column 6 : 'TLR_regions' : Number of good TLRs detected genome wide <br/>
    Column 7 : 'TLR_regions_contig_ends' : Number of the above TLRs found within the contig ends (Default: 50kb from end of contigs larger than 100kb) <br/>
    Column 8 : 'TLR_representative_coords' : The coordinates in the assembly to the TLR representative used in the analyses <br/>
    Column 9 : 'TLR_representative_size' : The size in basepairs of the TLR representative <br/>
    Column 10 : 'TLR_representative_size' : The average size of the regions detected that clustered with the TLR representative <br/>
    Column 12 : 'TLR_representative' : The nucleotide sequence of the TLR representative <br/>

-Per sample plots are all placed in **plotting_Rscripts** <br/>
    plotting_Rscripts/\*.end_alignments.svg : Contigs with ends aligned with TLRs highlighted (red) and telomeric repeats showing (black dots) <br/>
    plotting_Rscripts/\*.end_alignments.filtered.svg : Same as above but only with contig ends containing at least one telomeric repeat or TLR <br/>
    plotting_Rscripts/\*.whole_genome.svg : Positions of telomeric_repeats (dots) and TLR_regions (triangles) in the whole genome <br/>

-Comparisons between samples are placed in output folder
    phylogeny_plus_gANI_heatmap.svg : a midrooted mashtree of all assemblies horizontally aligned with lz-ani calulated global ANI (gANI) scores comparing the representatives of each assemblies TLR. Far right is the same gANI scores but calculated between repeats for the same assembly.

## The coverage threshold
This value is essentially is used to find repeats that are at least in this many copies compared to telomeric sequences <br/>
This tools was designed primarily on Fusarium and Pyricularia assemblies which contain (if present) TLRs on most chromosome ends <br/>
However they may be less frequent in other species <br/>
The default of 0.75 is too make sure that other potentially repetitive regions within 50kb of the contig ends are not kept

## Manually identified a secondary repeat?
When looking a the alignments of the contig ends (and maybe manually reordering the ends etc), did you notice another region next to telomeres with several alignments? <br/>
Here is a few lines you can run in order to extract and detail another repeat <br/>
All you need to know is the coordinates for one of the repeats (will be considered the representative copy) <br/>
For Example: I extracted a secondary repeat in the assembly 'GCA030345115.fa' from contig 'GCA030345115_CP128282.1' between position 2 and 9881 (GCA030345115_CP128282.1:2-9881) <br/>
In my run of TLRFinder on 'GCA030345115.fa', the prefix was GCA030345115 and therefore in the TLRFinder output folder I ran the below to add the secondary repeat:

        ##my run of TLRFinder commented out
        #TLRFinder -a GCA030345115.fa -p GCA030345115 -o GCA030345115_TLRFinder_output
        
        prefix="GCA030345115"
        assembly="GCA030345115.fa"
        coords="GCA030345115_CP128282.1:2-9881"

        ##move into the TLRFinder output for 'GCA030345115.fa'
        cd GCA030345115_TLRFinder_output/
        
        ##manually extract the second repeat for ${prefix} based on using the synteny alignment script (only the primary repeat was found previously)
        samtools faidx ${assembly} ${coords} > subtelomeric_repeats/${prefix}.manual_second.repeat_rep.fa
        ##get the positions of this second repeat in the contig ends
        blastn  -subject contig_ends/${prefix}.${tipsize2}kb_ends.fa -query subtelomeric_repeats/${prefix}.manual_second.repeat_rep.fa  -outfmt 6 > subtelomeric_repeats/${prefix}.manual_second.repeat_rep.ends_blast.tsv
        ##create bed from nonredundant positions (used for plotting)
        echo "contig;start;end" | tr ';' '\t' > subtelomeric_repeats/${prefix}.manual_second.repeat_rep.ends_blast.bed
        cat subtelomeric_repeats/${prefix}.manual_second.repeat_rep.ends_blast.tsv | awk '{if($3 > 80 && $4 > 500) print}' | awk '{print $2"\t"$9"\t"$10}' | awk '{if($2>$3) {print $1"\t"$3"\t"$2} else {print}}' | sort -k1,1 -k2,2n | bedtools merge -d 10  >> subtelomeric_repeats/${prefix}.manual_second.repeat_rep.ends_blast.bed
        ##and the same now for the whole genome
        blastn  -subject assemblies/${prefix}.fa -query subtelomeric_repeats/${prefix}.manual_second.repeat_rep.fa  -outfmt 6 > subtelomeric_repeats/${prefix}.manual_second.repeat_rep.WG_blast.tsv
        echo "contig;start;end" | tr ';' '\t' > subtelomeric_repeats/${prefix}.manual_second.repeat_rep.WG_blast.bed
        cat subtelomeric_repeats/${prefix}.manual_second.repeat_rep.WG_blast.tsv | awk '{if($3 > 80 && $4 > 500) print}' | awk '{print $2"\t"$9"\t"$10}' | awk '{if($2>$3) {print $1"\t"$3"\t"$2} else {print}}' | sort -k1,1 -k2,2n | bedtools merge -d 10 >> subtelomeric_repeats/${prefix}.manual_second.repeat_rep.WG_blast.bed

        ##Now you can plot both side by side


## How to look for known functional domains
TLRFinder looks for the conserved region however the actual expressed portion of the repeat (generally a helicase) is only contained within the repeat <br/>
Therefore you may be interested to find domains within the repeat as an indication of the function <br/>
Notably, in [O'Donnell et al. 2025]() using TLRFinder we found several examples of repeats containing genes that were not helicases.

This step requires you to download a large database of conserved domains in order to search against <br/>

        ##first download the appropriate conserved domain library (https://pmc.ncbi.nlm.nih.gov/articles/PMC7378889/)
        wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz
        mkdir Cdd_LE
        mv Cdd_LE.tar.gz Cdd_LE
        cd Cdd_LE
        tar -zxvf Cdd_LE.tar.gz

        ##also download a library to read the cdd index numbers
        wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz
        gunzip cddid.tbl.gz
        cd ../

Then use rpstblastn to search against these domains with your representative TLR OR all TLRs found in your assembly dataset <br/>
Here we are using the TLR representative for GCA030345115 <br/>
We also just grab the top 5 for the blast results for looking at the function using the cdd index table (this can definitely be increased if necessary) 

        sample="GCA030345115"
        input="GCA030345115_TLRFinder_output/subtelomeric_repeats/GCA030345115.repeat_rep.fa"
        ##alternatively if having used multiple assemblies with -al you could use the combined 'TLRFinder_output/repeat_representatives.fa' file which contains all TLR representatives from assemblies where it was found
        tophits="5"

        ##get the outfmt 6 plus the actual sequence aigned in the query (our repeats)
        rpstblastn -query ${input} -db Cdd_LE/Cdd -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq" > ${sample}.repeat_rep.rpsblast_output.txt 
 
        ##now take the best 5 hits for each representative
        echo "repeat_representative;CDD;domain;title;description" | tr ';' '\t' > ${sample}.repeat_rep.rpsblast_output.besthit_description.tsv
        grep '>' ${input} | sed 's/>//g' | while read rep
        do
        grep "${rep}" ${sample}.repeat_rep.rpsblast_output.txt | head -n${tophits}
        done | while read line
        do
        cdd=$( echo "${line}" | cut -f2 | awk -F "|" '{print $3}' )
        cdd2=$( grep ^${cdd} Cdd_LE/cddid.tbl | awk -F "\t" '{print "CDD:"$1"\t"$2"\t"$3"\t"$4}' )
        echo "${line}" | awk -v cdd2="$cdd2" '{print $1"\t"cdd2}'
        done  >> ${sample}.repeat_rep.rpsblast_output.besthit_description.tsv



