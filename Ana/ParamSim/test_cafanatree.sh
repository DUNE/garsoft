#!/bin/bash

#bin/mcp_skimmer --infile anatree_100.root --outfile skimmed.root --debug 0
#echo "done with running skimming module"
#bin/cafanatree_module --infile skimmed.root --outfile caf.root
#bin/cafanatree_module --infile andytree.root --outfile caf.root

cafanatree_module --infile andytree.root --outfile caf.root

echo "done with running cafanatree module"
