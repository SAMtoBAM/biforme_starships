
##environment
#conda config --add channels bioconda
#conda config --add channels conda-forge
#conda create -n stargraph samtobam::stargraph
source ~/.bashrc
conda activate stargraph

threads=32

project="biforme"
cd /home/samuel/projects/${project}



############################################################################################################
#################################### 1. Assembly set up
############################################################################################################

####was given absolute paths to all camemberti, caseifulvum and birforme assemblies
##removed the duplicates i.e. the short read assemblies that have a long-read equivalent
##saved the paths as assembly_paths.txt

##however these assemblies will need to be renamed with the semi-pan_spec contig naming convention
##this means a sample identifier, a common separator, then the contig name

##therefore I need to copy the assemblies here with the renaming
##also going to filter for a minimum size of contig

contigmin=1000

mkdir genomes
mkdir genomes/SR
mkdir genomes/SR/Pfus
mkdir genomes/SR/Pbif
mkdir genomes/SR/Pcam
mkdir genomes/SR/Pcas
mkdir genomes/LR

## first for the short read assemblies
##renaming the contigs from "NODE_#" to "contig#"
cat assembly_paths.txt | grep -v "genome_graph" | while read assembly
do
species=$( echo $assembly | awk -F "/" '{print $NF}' | awk -F "_" '{print $1}' | sed 's/Pcamcam/Pcam/g' | sed 's/Pcamcas/Pcas/g' | sed 's/Pbifc/Pbif/' | sed 's/Pbifn/Pbif/' )
sample=$( echo $assembly | awk -F "/" '{print $NF}' | awk -F "_" '{print $2}' | awk -F "." '{print $1}' )
cat $assembly | awk -F "_" -v sample="$sample" '{if($0 ~ ">"){print ">"sample"_contig"$4} else {print}}' | seqkit seq -m ${contigmin} - > genomes/SR/${species}/${sample}.renamed_filtered.fa
done

##then the long-read assemblies
cat assembly_paths.txt | grep "genome_graph" | while read assembly
do
sample=$( echo $assembly | awk -F "/" '{print $NF}' | awk -F "." '{print $1}' )
cat $assembly | awk -F "_" -v sample="$sample" '{if($0 ~ ">"){print ">"sample"_"$3} else {print}}' > genomes/LR/${sample}.renamed.fa
done


############################################################################################################
#################################### 2. Stargraph
############################################################################################################



##first need to create a file with the paths to all the assemblies
ls genomes/*/*.fa genomes/*/*/*.fa  > assemblies_panSN.txt

##run the stargraph starfish wrapper
starfish_wrapper.sh -a assemblies_panSN.txt -t ${threads}


##maually curated the starships by looking at the pairviz and elementviz plots
##made a list of all elements to be removed 'starfish_elements_to_be_removed.txt'
cp starfish_output/elementFinder/starfish.elements.ann.feat starfish.elements.ann.feat.temp
cat starfish_elements_to_be_removed.txt  | awk -F "\t" '{print $1}' | while read element
do
sed -i "s/^${element}/TOBEREMOVED/g" starfish.elements.ann.feat.temp
done

grep -v "TOBEREMOVED" starfish.elements.ann.feat.temp > starfish.elements.ann.feat.manual_curation.tsv
rm starfish.elements.ann.feat.temp


###now to run stargraph
##need a list of only the long-read assemblies
ls genomes/LR/*.fa  > assemblies_panSN.LR.txt
##and need to reduce the input data to only the elements found in these assemblies (otherwise it'll screw with some downstream analyses where it tries to extract the sequences)
head -n1 starfish.elements.ann.feat.manual_curation.tsv > starfish.elements.ann.feat.manual_curation.LR.tsv
cat assemblies_panSN.LR.txt | awk -F "/" '{print $NF}' | awk -F "." '{print $1}' | while read sample
do
grep "^${sample}_" starfish.elements.ann.feat.manual_curation.tsv
done >> starfish.elements.ann.feat.manual_curation.LR.tsv

##do the opposite as well (only SR)
head -n1 starfish.elements.ann.feat.manual_curation.tsv > starfish.elements.ann.feat.manual_curation.SR.tsv
grep -v -F -x -f starfish.elements.ann.feat.manual_curation.LR.tsv starfish.elements.ann.feat.manual_curation.tsv >> starfish.elements.ann.feat.manual_curation.SR.tsv




##now run stargraph with everything
stargraph.sh -a assemblies_panSN.LR.txt -r starfish_output/starfish.filt.SRGs_combined.gff -e starfish.elements.ann.feat.manual_curation.LR.tsv -t ${threads}


##use the combined dataset to look for HGT using cargobay.sh

##can use all assemblies from the concatenated file in starfish (-a) and all SRGs from starfish (-g)

##have to create a metadata file with the species assignments (-m)
##can use the genome folders they were placed into with the species names, just modifying them to the full name
cat assemblies_panSN.txt | awk -F "/" '{if($(NF-2) == "SR") {print $NF"."$(NF-1)} else if($NF ~ "ESE00090") {print $NF".Penicillium fuscoglaucum"} else if($NF ~ "ESE00019") {print $NF".Penicillium caseifulvum"} else if($NF ~ "LCP06093") {print $NF".Penicillium camemberti"} else {print $NF".Penicillium biforme"} }' | awk -F "." '{print $1"\t"$NF}' | sed 's/Pfus/Penicillium fuscoglaucum/g' | sed 's/Pbif/Penicillium biforme/g' | sed 's/Pcas/Penicillium caseifulvum/g' | sed 's/Pcam/Penicillium camemberti/g' > species_metadata.all_assemblies.tsv


##need to combine the stargraph output with the short-read output from starfish in a fasta (-e) and bed (-b) file
##first create a bed file just for the SR cals
tail -n+2 starfish.elements.ann.feat.manual_curation.SR.tsv | awk -F "\t" '{print $4"\t"$6"\t"$7"\t"$1}' > starfish.elements.ann.feat.manual_curation.SR.bed
##now combine that with the stargraph final bed file 'stargraph_output/4.SLR_starship_combination/*.starships_SLRs.bed'
cat stargraph_output/4.SLR_starship_combination/*.starships_SLRs.bed starfish.elements.ann.feat.manual_curation.SR.bed > starfish_stargraph.LR_SR_combined.bed
##now extract all those regions from the full assembly set
bedtools getfasta -name -bed starfish_stargraph.LR_SR_combined.bed -fi starfish_output/blastdb/starfish.assemblies.fa | sed 's/::/ /g' > starfish_stargraph.LR_SR_combined.fa


##running cargobay to detect HGT and using a strict containment threshold of 90% for the initial candidate identification
cargobay.sh -g starfish_output/starfish.filt.SRGs_combined.gff -a starfish_output/blastdb/starfish.assemblies.fa -m species_metadata.all_assemblies.tsv -e starfish_stargraph.LR_SR_combined.fa -b starfish_stargraph.LR_SR_combined.bed -t ${threads} -c 0.9


###can we calculate the prevelance of elements to transfer to other species
##here we need to consider the number of genomes per species
echo "element;species;present;absent" | tr ';' '\t' > cargobay.species_enrichment.tsv
tail -n+2 cargobay_output/1.database_search/cargobay.sourmash_multisearch.candidates.final.tsv  | awk -F "\t" '{print $1}' | sort -u | while read element
do
tail -n+2 cargobay_output/1.database_search/cargobay.sourmash_multisearch.candidates.final.tsv  | awk -F "\t" -v element="$element" '{if($1 == element) print $3}' | sort -u | while read species
do
##get count of genomes with this element
present=$( tail -n+2 cargobay_output/1.database_search/cargobay.sourmash_multisearch.candidates.final.tsv  | awk -F "\t" -v element="$element" -v species="$species" '{if($1 = element && $3 == species) {print} }' | wc -l )
##rename the species to be found in the cargobay databse and then get the total count
species2=$( echo $species | sed 's/_/ /g' )
##and make sure it is just for your species
total=$( grep ",${species2}" cargobay_output/0.database/lineages.fungi.csv | awk -F "," -v species2="${species2}" '{if($NF ~ species2) {print $NF} else {print $(NF-1)}}' | sed "s/${species2}/${species}/g" | awk -F " " '{print $1}' | grep ${species}$ | wc -l )
echo "${element};${species};${present};${total}" | tr ';' '\t' | awk '{if($3 != "") print $1"\t"$2"\t"$3"\t"$4-$3}'
done
done | grep -v "Penicillium_sp." >> cargobay.species_enrichment.tsv



###########
######USE THE FULL DATASET TO DEFINE CLUSTERS FOR THE FULL DATASET
##using the same scheme used by stargraph but now with all elements
## sourmash k-mer containment using k-31, 0.3 minimum similarity, and mcl clustering
kmerthreshold="0.3"
mkdir SLR_starship_network_clustering
cd SLR_starship_network_clustering

##generate a folder to place the sourmash singatures etc
mkdir sourmash_signatures

##for sourmash sketching there is one important parameter to think about: k-mer size
##here we are using a k-mer size of 31, this gives stringency and resolution
awk -F " " '{print $1}' ../starfish_stargraph.LR_SR_combined.fa > temp.fa
sourmash sketch dna -p k=31,noabund --singleton -o sourmash_signatures/ temp.fa
##first we can get the jaccard similarity
sourmash compare -k 31 sourmash_signatures/*.sig.gz --csv sourmash_signatures.compare_k31.jaccard.csv
##and second we can get the max pairwise containment score
sourmash compare -k 31 sourmash_signatures/*.sig.gz --max-containment --csv sourmash_signatures.compare_k31.containment.csv

rm temp.fa

##generate header for the file that will be used to build the network
echo "to;from;weight" | tr ';' '\t' > starships_SLRs.pairwise_jaccard.tsv
##now get a list of nonredundant pairwise jaccard similarities 
cat sourmash_signatures.compare_k31.jaccard.csv | tr -d '\r'  | awk -F',' 'NR==1{for(i=1;i<=NF;i++)samples[i]=$i;next}{row=NR-1;for(i=row+1;i<=NF;i++)print samples[row],samples[i],$i}' OFS='\t' | awk -F "\t" -v kmerthreshold="$kmerthreshold" '{if($3 >= kmerthreshold) {print}}' >> starships_SLRs.pairwise_jaccard.tsv

##same but for the containment scores (we used max containment so the pairwise values are symmetric making this easy)
echo "to;from;weight" | tr ';' '\t' > starships_SLRs.pairwise_containment.tsv
##now get a list of nonredundant pairwise jaccard similarities 
cat sourmash_signatures.compare_k31.containment.csv | tr -d '\r'  | awk -F',' 'NR==1{for(i=1;i<=NF;i++)samples[i]=$i;next}{row=NR-1;for(i=row+1;i<=NF;i++)print samples[row],samples[i],$i}' OFS='\t' | awk -F "\t" -v kmerthreshold="$kmerthreshold" '{if($3 >= kmerthreshold) {print}}' >> starships_SLRs.pairwise_containment.tsv


##also want a simplified metadata file used for plotting the networks
##first create a template for all the metadata as done in stargraph but now for all elements
echo "starship_SLR;navis-haplotype;contig;start;end;size;captain;captain_start;captain_end;captain_sense" | tr ';' '\t' > starships_SLRs.tsv
tail -n+2 ../starfish.elements.ann.feat.manual_curation.tsv | awk '{print $1"\t"$3"\t"$4"\t"$6"\t"$7"\t"$7-$6"\t"$5}' | while read line
do
tyr=$( echo "${line}" | awk -F "\t" '{print $NF}' )
cat ../starfish_output/starfish.filt.SRGs_combined.gff | grep "${tyr}"$ | awk -v line="$line" '{print line"\t"$4"\t"$5"\t"$7}'
done >> starships_SLRs.tsv
tail -n+2 ../stargraph_output/4.SLR_starship_combination/stargraph.SLRs.starships_subtracted.tyrRs_agg.tsv | cut -f1-9,11 >> starships_SLRs.tsv

##now add predefined navis information
#echo "name;type;family;navis_haplotype;navis;navis_slim" | tr ';' '\t' > ${prefix}.starships_SLRs.metadata.tsv
echo "name;type;navis_haplotype;navis;navis_slim" | tr ';' '\t' > starships_SLRs.metadata.tsv
tail -n+2 starships_SLRs.tsv | awk '{print $1}' | sort -u | while read element
do
navhap=$( tail -n+2 starships_SLRs.tsv  | awk -v element="$element" '{if($1 == element) print $2}' )
navis=$( echo "${navhap}" | awk -F "-" '{print $1}' )
navisslim=$( echo "${navis}" | awk '{if($1 ~ "nav"){print "NA"} else {print}}' )
#family=$( tail -n+2 ../${prefix}.starships_SLRs.tsv  | awk -v element="$element" '{if($1 == element) print $2}' )
type=$( echo "${element}" | awk '{if($1 ~ "SLR"){print "SLR"} else {print "Starship"}}' )
#echo "${element};${type};${family};${navhap};${navis};${navisslim}" | tr ';' '\t'
echo "${element};${type};${navhap};${navis};${navisslim}" | tr ';' '\t'
done >> starships_SLRs.metadata.tsv

##get rid of header before feeding the pairwise similarities to mcl
tail -n+2 starships_SLRs.pairwise_jaccard.tsv > temp.tsv
##cluster all the Starships and SLRs using the same kmer similarities (containment)
##now use mcl to quickly find the clusters
mcl temp.tsv --abc -o starships_SLRs.pairwise_jaccard.mcl.txt
rm temp.tsv
##now name the clusters and then append to the summary files
awk -F '\t' '{for (i=1; i <= NF; i++) {print "cluster"NR "\t" $i}}' starships_SLRs.pairwise_jaccard.mcl.txt > starships_SLRs.pairwise_jaccard.mcl.clusters.txt

##get rid of header before feeding the pairwise similarities to mcl
tail -n+2 starships_SLRs.pairwise_containment.tsv > temp.tsv
##cluster all the Starships and SLRs using the same kmer similarities (containment)
##now use mcl to quickly find the clusters
mcl temp.tsv --abc -o starships_SLRs.pairwise_containment.mcl.txt
rm temp.tsv
##now name the clusters and then append to the summary files
awk -F '\t' '{for (i=1; i <= NF; i++) {print "cluster"NR "\t" $i}}' starships_SLRs.pairwise_containment.mcl.txt > starships_SLRs.pairwise_containment.mcl.clusters.txt

echo "element;contig;start;end;cluster_contaiment;cluster_jaccard" | tr ';' '\t' > ${prefix}.starships_SLRs.plus_clusters.tsv
tail -n+2 starships_SLRs.tsv | awk '{print $1"\t"$3"\t"$4"\t"$5}' | while read line
do
SLR=$( echo "${line}" | awk '{print $1}' )
contclust=$( cat starships_SLRs.pairwise_containment.mcl.clusters.txt | awk -v SLR="$SLR" -v line="$line" '{if($2==SLR) print $1}' )
jaccclust=$( cat starships_SLRs.pairwise_jaccard.mcl.clusters.txt | awk -v SLR="$SLR" -v line="$line" '{if($2==SLR) print $1}' )
echo "${line};${contclust};${jaccclust}" | tr ';' '\t'
done | awk -F'\t' -v OFS='\t' '{ if ($5 == "") $5 = $1; if ($6 == "") $6 = $1; print }' >> starships_SLRs.plus_clusters.tsv

cd ..



###use starfish covereage to genotype al the short read datasets
##need a file with the sampleID in the first column and each paired end read dataset in the next two (-r)
##we can use a file with all the rwa paired end illumina paths in a txt file
cat /scratch/jp/P_biforme_fusco_camemberti/data/path_to_data_short_reads.txt | awk -F "/" '{print $NF}' | awk -F "_" '{print $2}' | sort -u | while read sample
do
R1=$( grep "${sample}_R1" /scratch/jp/P_biforme_fusco_camemberti/data/path_to_data_short_reads.txt )
R2=$( grep "${sample}_R2" /scratch/jp/P_biforme_fusco_camemberti/data/path_to_data_short_reads.txt )
echo "${sample};${R1};${R2}" | tr ';' '\t'
done > coverage_SR_paths.tsv

###UPDATED ~miniconda3/envs/stargraph/main/coverage script with the version from https://github.com/egluckthaler/starfish/pull/23

##now run the coverqge analysis
mkdir starfish_coverage
starfish coverage -m paired --aligner minimap2 -r coverage_SR_paths.tsv -l starfish_stargraph.LR_SR_combined.fa -x starfish_coverage -o starfish_coverage --clean --threads ${threads}


##now we need to extract the coverage output for each strain and each starship/SLR
##place it all into a single matrix with samples as the rows and elements as the columns
##use the covered amount as the variable to keep

awk -F'\t' '
# Skip header lines starting with #
$1 ~ /^#/ { next }

{
    sample=$1
    split($2,a," "); elem=a[1]
    cov=$4
    data[sample,elem]=cov
    samples[sample]=1
    elements[elem]=1
}
END {
    # Build sorted element list
    n=0
    for (e in elements) element_list[++n]=e
    asort(element_list)

    # Print header
    printf "sampleID"
    for (i=1;i<=n;i++) printf "\t%s", element_list[i]
    print ""

    # Print each sample row (sorted)
    m=0
    for (s in samples) sample_list[++m]=s
    asort(sample_list)

    for (i=1;i<=m;i++) {
        s=sample_list[i]
        printf "%s", s
        for (j=1;j<=n;j++) {
            e=element_list[j]
            printf "\t%s", ((s,e) in data ? data[s,e] : 0)
        }
        print ""
    }
}
' starfish_coverage/*.cov > starfish_coverage.summary_matrix.tsv



###now take the above matrix and summarise by the k-mer containment, mcl clusters
MATRIX="starfish_coverage.summary_matrix.tsv"
CLUSTERS="SLR_starship_network_clustering/starships_SLRs.plus_clusters.tsv"

awk -F'\t' -v CLUSTERS="$CLUSTERS" '
BEGIN {
    # Load cluster mapping
    lineNo = 0
    while ((getline line < CLUSTERS) > 0) {
        lineNo++
        if (lineNo == 1) continue  # skip header in clusters.tsv
        split(line, f, "\t")
        elem = f[1]
        cluster = f[5]
        elem2cluster[elem] = cluster
        clusters[cluster] = 1
    }
    close(CLUSTERS)
}

NR==1 {
    # Read header from matrix
    for (i=2; i<=NF; i++) header[i] = $i
    next
}

{
    sample = $1
    for (i=2; i<=NF; i++) {
        elem = header[i]
        val = $i
        if (elem in elem2cluster) {
            c = elem2cluster[elem]
            sum[sample, c] += val
            count[sample, c]++
            samples[sample] = 1
        }
    }
}

END {
    # Sort clusters
    n = 0
    for (c in clusters) clust_list[++n] = c
    asort(clust_list)

    # Header
    printf "sampleID"
    for (i=1; i<=n; i++) printf "\t%s", clust_list[i]
    print ""

    # Sort samples
    m = 0
    for (s in samples) sample_list[++m] = s
    asort(sample_list)

    # Output averaged matrix
    for (i=1; i<=m; i++) {
        s = sample_list[i]
        printf "%s", s
        for (j=1; j<=n; j++) {
            c = clust_list[j]
            avg = (count[s, c] ? sum[s, c] / count[s, c] : 0)
            printf "\t%.6f", avg
        }
        print ""
    }
}
' "$MATRIX" > starfish_coverage.summary_matrix.clusters.tsv




##path to short reads for all strains
/scratch/jp/P_biforme_fusco_camemberti/data/path_to_data_short_reads.txt
