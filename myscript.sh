#!/bin/sh

set -x
# . env.sh
if test -z $TREEMACHINE_ROOT 
then
	echo "Treemachine root must be defined"
	exit 1
fi
treemachine_jar="$TREEMACHINE_ROOT/target/treemachine-0.0.1-SNAPSHOT-jar-with-dependencies.jar"
if ! test -f "$treemachine_jar"
then
	echo "Treemachine jar not found"
	exit 1
fi
#declare variables
taxonomy_file="$1"
synonyms_file="$2"
newick_taxonomy="primates_ott_newick.txt"
primates_big_tree="primates_big_tree.txt" 
primates_true_tree_ids="primates_true_tree_ids.txt"
primates_true_tree_wo_ids="primates_true_tree_no_ids.txt"
primates_input_trees="primates_input_trees.txt"
big_tree_input_info="big_tree_leaf_overlap.txt"
primates_MRP="primates_MRP.txt"
primates_input_trees_nexon="primates_input_trees_nexon.txt"
primates_bogus_ott="primates_bogus_ott.txt"
primates_SAS_tree="primates_SAS_tree.tre"
primates_MRP_tree="primates_MRP_tree.tre"
distance_results="distances.txt"
original_dir="$PWD"

#load taxonomy and convert to newick format. 
taxonomyToNewick.py "${taxonomy_file}" > "$original_dir/${newick_taxonomy}" || exit
#head -n1 "$original_dir/${newick_taxonomy}" > newick_tax_ids.txt
#tail -n1 "$original_dir/${newick_taxonomy}" > newick_tax_noids.txt
#delete test directory if already present.
if test -d myoutput
then
	echo "myoutput is in the way!"
	exit 1
fi
mkdir myoutput
cd myoutput

for i in $(seq 1 1 10) # Parameter swoop begin.
do
	mkdir run$i
	#Run taxonomy through big-tree.
	big-tree-sim-0.0.2a -c "$original_dir/basic-commands.txt" "$original_dir/${newick_taxonomy}" > run$i/"${primates_big_tree}" || exit

	head -n1 run$i/"${primates_big_tree}" > run$i/"${primates_true_tree_ids}" || exit
	nl=$(wc -l < run$i/"${primates_big_tree}") || exit
	ninp=$(expr $nl - 1) || exit
	tail -n$ninp run$i/"${primates_big_tree}" > run$i/"${primates_input_trees}" || exit #; split -d -l 1 run$i/"${primates_input_trees}" run$i/input$i.txt || exit

	tree_leaf_overlap.py run$i/"${primates_input_trees}" > run$i/"${big_tree_input_info}"	

	name_split.py run$i/"${primates_true_tree_ids}"> run$i/"${primates_true_tree_wo_ids}" || exit
	#Run MRP algorithm (MRP output seems to be off still, no '?' appear... 
	MRP_Matrix_Generator.py "run$i/${primates_input_trees}" > run$i/"${primates_MRP}" || exit
	#Eventually run MRP file through PAUP to return a nexus file readable by dendropy..
	
	mrp_start=$(date +"%s")
	echo "Set MaxTrees = 100 Increase = No ;" >> "run$i/${primates_MRP}"
	echo "HSearch MulTrees = Yes ;" >> "run$i/${primates_MRP}"
	echo "SaveTrees file = run$i/${primates_MRP_tree} ;" >> "run$i/${primates_MRP}"
	echo "ConTree /strict save treefile = run$i/mrp_con_tree.txt ;" >> "run$i/${primates_MRP}"
	time paup -n "run$i/${primates_MRP}"
	mrp_end=$(date +"%s")
	
	#Convert to Nexon for TAG algorithm.
	#Set the field seperator to a newline


	#Create bogus ott for TAG.
	bogus_ott_gen.py "$original_dir/${newick_taxonomy}" > run${i}/"${primates_bogus_ott}" || exit
	#Run TAG algorithm (as Nexson) compare to bogus ott generated.

	#set primates[1-10].db to automatically delete each time to avoid confusion.
	if test -d run$i/primates$i.db
	then
		echo "primates$i.db is in the way! removing now.."
		rm -r run$i/primates$i.db
	fi

		
	#SAS Algorithm.
	sas_start=$(date +"%s")
	time java -jar "$treemachine_jar" inittax run${i}/"${primates_bogus_ott}" "$synonyms_file" run$i/primates$i.db || exit
	time java -jar "$treemachine_jar" addtaxonomymetadatanodetoindex 3 run$i/primates$i.db || exit
	time java -jar "$treemachine_jar" listsources run$i/primates$i.db || exit
	for j in $(seq 1 1 $ninp)
	do
		head -n$j run$i/"${primates_input_trees}" | tail -n1 > run$i/inp${j}.tre
		newick_to_nexon.py run$i/inp${j}.tre > run$i/inp${j}.nexson || exit
		time java -jar "$treemachine_jar" pgloadind run$i/primates$i.db run$i/inp${j}.nexson || exit
	done
	time java -jar "$treemachine_jar" synthesizedrafttree 805080 run$i/primates$i.db || exit
	time java -jar "$treemachine_jar" extractdrafttree 805080 run$i/"${primates_SAS_tree}" run$i/primates$i.db || exit
	sas_end=$(date +"%s")
	#Compare MRP and TAG distance.
	distance.py run$i/"${primates_true_tree_wo_ids}" run$i/"${primates_MRP_tree}" run$i/mrp_con_tree.txt run$i/"${primates_SAS_tree}" $i false_positives_and_negatives$ninp.txt > run$i/"${distance_results}" || exit
	mrp_diff=$(($mrp_end-$mrp_start))
	sas_diff=$(($sas_end-$sas_start))
	echo "total time for MRP: " >> run$i/"${distance_results}"
	date -u -d @"$mrp_diff" +'%-Hh %-Mm %-Ss' >> run$i/"${distance_results}"
	echo "total time for SAS: "  >> run$i/"${distance_results}"
	date -u -d @"$sas_diff" +'%-Hh %-Mm %-Ss' >> run$i/"${distance_results}"
done #Parameter swoop end.
#false_pos_and_neg.py false_positives_and_negatives$ninp.txt "$original_dir"/false_averages.txt $ninp
#mv "$original_dir"/myoutput "$original_dir"/myoutput$ninp