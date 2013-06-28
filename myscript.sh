#!/bin/sh
#load taxonomy and convert to newick format. 
set -x
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
taxonomy_file="$1"
synonyms_file="$2"
newick_taxonomy="primates_ott_newick.txt"
primates_bigtree="primates_big_tree.txt" # fake 'true' tree. has ottolids attached.
primates_MRP="primates_MRP.txt"
primates_bigtree_nexon="primates_bt_nexon.txt"
primates_bogus_ott="primates_bogus_ott.txt"
primates_SAS_tree="primates_SAS_tree.tre"
primates_MRP_tree="primates_MRP_tree.tre"
distance_results="distances.txt"
original_dir="$PWD"

taxonomyToNewick.py "${taxonomy_file}" > "$original_dir/${newick_taxonomy}" || exit

if test -d myoutput
then
	echo "myoutput is in the way!"
	exit 1
fi
mkdir myoutput
cd myoutput

for i in $(seq 1 1 10) # Parameter swoop begin
do
	mkdir run$i
	#for j in $(seq 1 1 1000) # probably will have to replace i with j in sequence's later..
	#do
		

		#Run taxonomy through big-tree.
		big-tree-sim-0.0.2a -c "$original_dir/basic-commands.txt" "$original_dir/${newick_taxonomy}" > run$i/"${primates_bigtree}" || exit # fake 'true' tree output here. 
		#Run MRP algorithm (MRP output seems to be off still, no '?' appear... Is it neccessary to get rid of ottolids in names?.
		MRP_Matrix_Generator.py "run$i/${primates_bigtree}" > run$i/"${primates_MRP}" || exit
		#Eventually run MRP file through PAUP to return a nexus file readable by dendropy..
		echo "HSearch;" >> run$i/"${primates_MRP}"
		echo "SaveTrees file = run$i/${primates_MRP_tree} ;" >> "run$i/${primates_MRP}"
		time paup -n "run$i/${primates_MRP}"
		#Convert to Nexon for TAG algorithm.
		newick_to_nexon.py run$i/"${primates_bigtree}" > run$i/"${primates_bigtree_nexon}" || exit
		#Create bogus ott for TAG.
		bogus_ott_gen.py "$original_dir/${newick_taxonomy}" > run${i}/"${primates_bogus_ott}" || exit
		#Run TAG algorithm (as Nexson) compare to bogus ott generated.

		#set primates1.db to automatically delete each time to avoid confusion.
		if test -d run$i/primates$i.db
		then
			echo "primates$i.db is in the way! removing now.."
			rm -r run$i/primates$i.db
		fi

		
		#SAS Algorithm.
		java -jar "$treemachine_jar" inittax run$i/"${primates_bogus_ott}" "$synonyms_file" run$i/primates$i.db || exit
		java -jar "$treemachine_jar" addtaxonomymetadatanodetoindex 3 run$i/primates$i.db || exit
		java -jar "$treemachine_jar" listsources run$i/primates$i.db || exit
		java -jar "$treemachine_jar" pgloadind run$i/primates$i.db run$i/"${primates_bigtree_nexon}" || exit
		java -jar "$treemachine_jar" synthesizedrafttreelist 805080 _0_0 run$i/primates$i.db || exit
		java -jar "$treemachine_jar" extractdrafttree 805080 run$i/"${primates_SAS_tree}" run$i/primates$i.db || exit

		#Compare MRP and TAG distance. Have to do something about single quote..
		#distance.py run$i/"${primates_bigtree}" run$i/"${primates_MRP_tree}" run$i/"${primates_SAS_tree}" > run$i/"${distance_results}" || exit
	#done
done
