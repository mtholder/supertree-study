#load taxonomy and convert to newick format. 
python ~/supertree-study/taxonomyToNewick.py ~/Dropbox/holder_lab/supertrees/primate_taxonomy_with_synonyms/primates_taxonomy.txt > ~/myfile/primates_ott_newick1.txt || exit
for i in $(seq 1 1 10) # Parameter swoop begin
do
	#for j in $(seq 1 1 1000) # probably will have to replace i with j in sequence's later..
	#do
		#Set each folder to automatically delete..
		if test -d ~/myfile/run$i
		then
			echo "the folder run$i is in the way and is being removed.."
			rm -r ~/myfile/run$i
		fi
		mkdir ~/myfile/run$i
		

		#Run taxonomy through big-tree.
		time ~/big-tree-collection-simulator/src/big-tree-sim-0.0.2a -c ~/big-tree-collection-simulator/test/basic-commands.txt ~/myfile/primates_ott_newick1.txt > ~/myfile/run$i/primates_big_tree_TAG$i.txt || exit # fake 'true' tree output here. 
		#Run MRP algorithm (MRP output seems to be off still, no '?' appear... Is it neccessary to get rid of ottolids in names?.
		time python ~/supertree-study/MRP_Matrix_Generator.py ~/myfile/run$i/primates_big_tree_TAG$i.txt > ~/myfile/run$i/primates_MRP$i.txt || exit
		#Eventually run MRP file through PAUP to return a nexus file readable by dendropy..
		# Not working correctly, need to find a way to run all commands through the terminal...
		time paup -n ~/myfile/run$i/primates_MRP$i.txt -l ~/myfile/run$i/primates_MRP$i || exit
		#time paup -r ~/myfile/run$i/primates_MRP$i.txt -l ~/myfile/run$i/primates_MRP$i || exit	
		#GenerateTrees || exit
		#RootTrees || exit
		#time SaveTrees || exit
		#quit || exit
		#Convert to Nexon for TAG algorithm.
		time python ~/supertree-study/newick_to_nexon.py ~/myfile/run$i/primates_big_tree_TAG$i.txt > ~/myfile/run$i/primates_bt_nexon$i.txt || exit
		#Create bogus ott for TAG.
		time python ~/supertree-study/bogus_ott_gen.py ~/myfile/primates_ott_newick1.txt > ~/myfile/run$i/primates_bogus_ott$i.txt || exit
		#Run TAG algorithm (as Nexson) compare to bogus ott generated.

		#set primates1.db to automatically delete each time to avoid confusion.
		if test -d ~/myfile/run$i/primates$i.db
		then
			echo "primates$i.db is in the way! removing now.."
		rm -r ~/myfile/primates$i.db
		fi

		#if test -z "$TRACE_LEVEL"
		#then
		#	config_arg_val=debuglog4j.properties
		#else
		#	config_arg_val=tracelog4j.properties
		#fi

		#set -x
		
		#SAS Algorithm.
		java -jar ~/treemachine/target/treemachine-0.0.1-SNAPSHOT-jar-with-dependencies.jar inittax ~/myfile/run$i/primates_bogus_ott$i.txt ~/Dropbox/holder_lab/supertrees/primate_taxonomy_with_synonyms/primates_synonyms.txt ~/myfile/run$i/primates$i.db || exit
java -jar ~/treemachine/target/treemachine-0.0.1-SNAPSHOT-jar-with-dependencies.jar addtaxonomymetadatanodetoindex 3 ~/myfile/run$i/primates$i.db || exit
		java -jar ~/treemachine/target/treemachine-0.0.1-SNAPSHOT-jar-with-dependencies.jar listsources ~/myfile/run$i/primates$i.db || exit
		java -jar ~/treemachine/target/treemachine-0.0.1-SNAPSHOT-jar-with-dependencies.jar pgloadind ~/myfile/run$i/primates$i.db ~/myfile/run$i/primates_bt_nexon$i.txt || exit
		java -jar ~/treemachine/target/treemachine-0.0.1-SNAPSHOT-jar-with-dependencies.jar synthesizedrafttreelist 805080 _0_0 ~/myfile/run$i/primates$i.db || exit
		java -jar ~/treemachine/target/treemachine-0.0.1-SNAPSHOT-jar-with-dependencies.jar extractdrafttree 805080 ~/myfile/run$i/primates_tree$i.tre ~/myfile/run$i/primates$i.db || exit

		#Compare MRP and TAG distance.
		time python ~/supertree-study/distance.py ~/myfile/run$i/primates_big_tree_TAG$i.txt ~/myfile/run$i/primates_MRP$i.tre ~/myfile/run$i/primates_tree$i.tre > ~/myfile/run$i/symmetric_differences$i.txt || exit
	#done
done
