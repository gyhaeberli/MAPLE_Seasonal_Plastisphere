## Post DADA2 and LULU Curation -> Swarm clustering algorithm

############# Set Up ###################
#$ conda activate MetabMaple #Metabarcoding environment
#$ conda create -n swarm python=2.7 
#$ conda activate swarm
#$ conda install bioconda::swarm
#$ conda install bioconda::vsearch
#######################################
cd data/glennsdata/MAPLE/18S/EUK_data/X204SC25074012-Z01-F001/01.RawData/swarm/

#Normal swarm
swarm nP_swarm.fasta > nP_post_swarm_clusters.txt
swarm psP_swarm.fasta > psP_post_swarm_clusters.txt  
swarm fP_swarm.fasta > fP_post_swarm_clusters.txt

#Swarn to get representative sequences + clusters
## USE THIS ONE

# Use Swarm to get representative sequences directly
swarm -d 1 -w nP_otu_repseq.fasta -o nP_post_swarm_clusters.txt -s nP_swarm_stats.txt nP_swarm.fasta
swarm -d 1 -w psP_otu_repseq.fasta -o psP_post_swarm_clusters.txt -s psP_swarm_stats.txt psP_swarm.fasta  
swarm -d 1 -w fP_otu_repseq.fasta -o fP_post_swarm_clusters.txt -s fP_swarm_stats.txt fP_swarm.fasta

