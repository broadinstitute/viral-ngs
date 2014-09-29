#-------- MULTIPLE-SEQUENCE ALIGNMENT & TREE BUILDING @ BROAD --------#
# Make sure all file-names are unique when cut down to 10 characters - e.g. if analysing Lassa 'LASV-' identifier needs to be removed from the input sequence file.

# MAKE REQUIRED SUB-DIRECTORIES
for directory in
do
bsub -o ~/log.txt -P sabeti_trees "mkdir $directory/_msa $directory/_sequences $directory/_trees $directory/_logs $directory/_temp $directory/_trees/mr-bayes"
done

#-------- SEQUENCE ALIGNMENT --------#
# CREATE MSA USING MAFFT & TRIM USING TRIMAL
for sequences in
do
for directory in
do
bsub -R "rusage[mem=4]" -n 4 -R "span[hosts=1]" -W 4:00 -q hour -o $directory/_logs/$sequences.log.bsub.txt -P sabeti_trees -J $sequences.al1 "/idi/sabeti-scratch/kandersen/bin/mafft/core/mafft --localpair --maxiterate 1000 --reorder --ep 0.123 --preservecase --thread 4 $directory/_sequences/$sequences.fasta > $directory/_msa/$sequences.mafft.fasta && /idi/sabeti-scratch/kandersen/bin/trimal/trimal -phylip -automated1 -in $directory/_msa/$sequences.mafft.fasta -out $directory/_msa/$sequences.pruned.phy -htmlout $directory/_logs/$sequences.log.trimal.html -colnumbering"
done
done

#-------- TREE BUILDING USING MAXIMUM LIKELIHOOD - RAXML --------#
# CREATE TREE
for sequences in
do
for directory in
do
for substitution_model in GTRGAMMA
do
for bootstraps in 500
do
bsub -sp 100 -R "rusage[mem=2]" -n 2 -R "span[hosts=1]" -W 4:00 -q hour -o $directory/_logs/$sequences.log.bsub.txt -P sabeti_trees -J $sequences.tr1 "/idi/sabeti-scratch/kandersen/bin/raxml/raxmlHPC-PTHREADS-SSE3 -f d -T 2 -p 123421 -m $substitution_model -N 20 -n $sequences.tree1 -w $directory/_trees/ -s $directory/_msa/$sequences.pruned.phy"
bsub -sp 90 -n 4 -R "span[hosts=1]" -q week -o $directory/_logs/$sequences.log.bsub.txt -P sabeti_trees -J $sequences.tr2 "/idi/sabeti-scratch/kandersen/bin/raxml/raxmlHPC-PTHREADS-SSE3 -f d -T 4 -p 12438 -m $substitution_model -b 12438 -N $bootstraps -k -n $sequences.tree2 -w $directory/_trees/ -s $directory/_msa/$sequences.pruned.phy && /idi/sabeti-scratch/kandersen/bin/raxml/raxmlHPC-SSE3 -T 1 -m $substitution_model -n $sequences.tree3 -f b -t $directory/_trees/RAxML_bestTree.$sequences.tree1 -z $directory/_trees/RAxML_bootstrap.$sequences.tree2 -w $directory/_trees/ && mv $directory/_trees/RAxML_bipartitions.$sequences.tree3 $directory/_trees/$sequences.raxml.tree"
done
done
done
done

#-------- CREATE TREES USING PHYML --------#
# CREATE TREE
for sequences in
do
for directory in
do
for substitution_model in GTR # HKY85, JC69, K80, F81, F84, TN93
do
for bootstraps in 500
do
bsub -n 1 -R "span[hosts=1]" -q week -o $directory/_logs/$sequences.log.bsub.txt -P sabeti_trees -J $sequences.pm "/idi/sabeti-scratch/kandersen/bin/phyml/PhyML-3.1_linux64 -i $directory/_msa/$sequences.pruned.phy -d nt -b $bootstraps -m $substitution_model --pinv 0 --nclasses 4 -s BEST --rand_start --n_rand_starts 10 --r_seed 1553 -f m --no_memory_check && mv $directory/_msa/$sequences.pruned.phy_phyml* $directory/_trees/ && /idi/sabeti-scratch/kandersen/bin/raxml/raxmlHPC-SSE3 -f b -t $directory/_trees/$sequences.pruned.phy_phyml_tree.txt -z $directory/_trees/$sequences.pruned.phy_phyml_boot_trees.txt -m GTRGAMMA -s $directory/_msa/$sequences.pruned.phy -n $sequences.phyml.tree -w $directory/_trees/ && mv $directory/_trees/RAxML_bipartitions.$sequences.phyml.tree $directory/_trees/$sequences.phyml.tree"
done
done
done
done

#-------- CREATE TREES USING MR BAYES --------#
# OUTPUT NEXUS FILE
for sequences in
do
for directory in
do
bsub -R "rusage[mem=2]" -n 1 -R "span[hosts=1]" -W 4:00 -q hour -o $directory/_logs/$sequences.log.bsub.txt -P sabeti_trees -J $sequences.pr "/idi/sabeti-scratch/kandersen/bin/trimal/trimal -automated1 -in $directory/_msa/$sequences.pruned.phy -nexus -out $directory/_msa/$sequences.pruned.nex"
done
done

# RUN MR BAYES
for sequences in
do
for directory in
do
bsub -R "rusage[mem=2]" -n 1 -R "span[hosts=1]" -q week -o $directory/_logs/$sequences.log.bsub.txt -P sabeti_trees -J $sequences.mb "cp $directory/_msa/$sequences.pruned.nex $directory/_trees/mr-bayes/$sequences.nex && /idi/sabeti-scratch/kandersen/bin/mrbayes/run_mrbayes.sh $directory/_trees/mr-bayes/$sequences.nex"
done
done