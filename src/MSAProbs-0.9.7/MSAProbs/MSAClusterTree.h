/***********************************************
 * # Copyright 2009-2010. Liu Yongchao
 * # Contact: Liu Yongchao, School of Computer Engineering,
 * #			 Nanyang Technological University.
 * # Emails:	 liuy0039@ntu.edu.sg; nkcslyc@hotmail.com
 * #
 * # GPL version 3.0 applies.
 * #
 * ************************************************/

#ifndef _MSA_CLUSTER_TREE_H
#define _MSA_CLUSTER_TREE_H

#include "MSAGuideTree.h"

class MSAClusterTree: public MSAGuideTree {
public:
	MSAClusterTree(MSA* msa, VVF& distMatrix, int numSeqs);
	~MSAClusterTree();

	//construct the cluster tree
	void create();
private:
	//generate the cluster tree
	void generateClusterTree();
};
#endif
