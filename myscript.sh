#!/bin/sh
#load taxonomy and convert to newick format. 
set -x

taxonomy_file="$1"
newick_taxonomy="primates_ott_newick1.txt"
taxonomyToNewick.py "${taxonomy_file}" > "${newick_taxonomy}" || exit

if test -d myoutput
then
	echo "myoutput is in the way !"
	exit 1
fi

for i in $(seq 1 1 10) # Parameter swoop begin
do
	#for j in $(seq 1 1 1000) # probably will have to replace i with j in sequence's later..
	#do
		#Set each folder to automatically delete..
		mkdir myoutput/run$i
		

		#Run taxonomy through big-tree.
		time ~/big-tree-collection-simulator/src/big-tree-sim-0.0.2a -c ~/big-tree-collection-simulator/test/basic-commands.txt "${newick_taxonomy}" > myoutput/run$i/primates_big_tree_TAG$i.txt || exit # fake 'true' tree output here. 
		#Run MRP algorithm (MRP output seems to be off still, no '?' appear... Is it neccessary to get rid of ottolids in names?.
		time python ~/supertree-study/MRP_Matrix_Generator.py myoutput/run$i/primates_big_tree_TAG$i.txt > myoutput/run$i/primates_MRP$i.txt || exit
		#Eventually run MRP file through PAUP to return a nexus file readable by dendropy..
		# Not working correctly, need to find a way to run all commands through the terminal...
		echo "HSearch;" >> myoutput/run$i/primates_MRP$i.txt
		echo "SaveTrees file = myoutput/run${i}/primates_MRP${i}_out.txt ;" >> myoutput/run$i/primates_MRP$i.txt
		

		time paup -n myoutput/run$i/primates_MRP$i.txt || exit
		#time paup -r myoutput/run$i/primates_MRP$i.txt -l myoutput/run$i/primates_MRP$i || exit	
		#Convert to Nexon for TAG algorithm.
		#time python ~/supertree-study/newick_to_nexon.py myoutput/run$i/primates_big_tree_TAG$i.txt > myoutput/run$i/primates_bt_nexon$i.txt || exit
		#Create bogus ott for TAG.
		#time python ~/supertree-study/bogus_ott_gen.py "${newick_taxonomy}" > myoutput/run$i/primates_bogus_ott$i.txt || exit
		#Run TAG algorithm (as Nexson) compare to bogus ott generated.

		#set primates1.db to automatically delete each time to avoid confusion.
		if test -d myoutput/run$i/primates$i.db
		then
			echo "primates$i.db is in the way! removing now.."
			rm -r myoutput/primates$i.db
		fi

		#if test -z "$TRACE_LEVEL"
		#then
		#	config_arg_val=debuglog4j.properties
		#else
		#	config_arg_val=tracelog4j.properties
		#fi

		#set -x
		
		#SAS Algorithm.
		java -jar ~/treemachine/target/treemachine-0.0.1-SNAPSHOT-jar-with-dependencies.jar inittax myoutput/run$i/primates_bogus_ott$i.txt ~/Dropbox/holder_lab/supertrees/primate_taxonomy_with_synonyms/primates_synonyms.txt myoutput/run$i/primates$i.db || exit
		java -jar ~/treemachine/target/treemachine-0.0.1-SNAPSHOT-jar-with-dependencies.jar addtaxonomymetadatanodetoindex 3 myoutput/run$i/primates$i.db || exit
		java -jar ~/treemachine/target/treemachine-0.0.1-SNAPSHOT-jar-with-dependencies.jar listsources myoutput/run$i/primates$i.db || exit
		java -jar ~/treemachine/target/treemachine-0.0.1-SNAPSHOT-jar-with-dependencies.jar pgloadind myoutput/run$i/primates$i.db myoutput/run$i/primates_bt_nexon$i.txt || exit
		java -jar ~/treemachine/target/treemachine-0.0.1-SNAPSHOT-jar-with-dependencies.jar synthesizedrafttreelist 805080 _0_0 myoutput/run$i/primates$i.db || exit
		java -jar ~/treemachine/target/treemachine-0.0.1-SNAPSHOT-jar-with-dependencies.jar extractdrafttree 805080 myoutput/run$i/primates_tree$i.tre myoutput/run$i/primates$i.db || exit

		#Compare MRP and TAG distance.
		time python ~/supertree-study/distance.py myoutput/run$i/primates_big_tree_TAG$i.txt myoutput/run$i/primates_MRP$i.tre myoutput/run$i/primates_tree$i.tre > myoutput/run$i/symmetric_differences$i.txt || exit
	#done
done
