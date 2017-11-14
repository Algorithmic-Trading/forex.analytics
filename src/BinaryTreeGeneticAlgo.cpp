#include "../include/BinaryTreeGeneticAlgo.h"

BinaryTreeGeneticAlgo::BinaryTreeGeneticAlgo(
	int selectNumber,
	double leafValueMutation,
	double leafValueSignMutation,
	double logicalNodeMutation,
	double leafValueIndicatorMutation,
	double crossover) {
	this->select = selectNumber;
	this->leafValueMutationProbability = leafValueMutation;
	this->leafValueSignMutationProbability = leafValueSignMutation;
	this->logicalNodeMutationProbability = logicalNodeMutation;
	this->crossoverProbability = crossover;
	this->leafIndicatorMutationProbability = leafValueIndicatorMutation;
}

BinaryTreeGeneticAlgo::~BinaryTreeGeneticAlgo() {
}

void BinaryTreeGeneticAlgo::Select(
	std::vector<BinaryTreeChromosome*>* newGeneration,
	std::vector<BinaryTreeChromosome*>* oldGeneration,
	unsigned size) {

	std::default_random_engine generator;
	std::uniform_int_distribution<unsigned> distribution(0, this->select);

	for (unsigned i = size - 1, y = 0; i >= size - this->select; i--, y++) {
		newGeneration->at(y)->setFitness(oldGeneration->at(i)->getFitness());
		oldGeneration->at(i)->copyTo(newGeneration->at(y));

	}

	for (unsigned i = this->select; i < size; i++) {
		unsigned index = distribution(generator);
		this->Mutate(newGeneration, index, newGeneration->at(i));
	}

	this->Crossover(newGeneration);
}

void BinaryTreeGeneticAlgo::Mutate(
	std::vector<BinaryTreeChromosome*>* generation,
	unsigned index,
	BinaryTreeChromosome * outputChromosome) {

	generation->at(index)->copyTo(outputChromosome);

	outputChromosome->Mutate(leafValueMutationProbability,
		leafValueSignMutationProbability,
		logicalNodeMutationProbability,
		crossoverProbability,
		leafIndicatorMutationProbability);
}

void BinaryTreeGeneticAlgo::Crossover(std::vector<BinaryTreeChromosome*>* generation) {
	for (unsigned long i = this->select + 1; i < generation->size(); i++) {
		double crossover = static_cast<double>(rand() % 100) / 100;

		if (crossover <= this->crossoverProbability) {
			unsigned long rnd = i;

			while (rnd == i || rnd <= this->select)
				rnd = rand() % generation->size();

			this->Crossover(generation->at(i), generation->at(rnd));
		}
	}
}

TreeNode * CutTree(
	TreeNode* node,
	unsigned int cuttingPoint,
	unsigned int i);
TreeNode * CutTree(
	TreeNode* node,
	unsigned int cuttingPoint,
	unsigned int i) {
	if (node->left != nullptr && node->right != nullptr) {
		if (i >= cuttingPoint)
			return node;
		if (rand() % 2 == 0) {
			TreeNode* tryCut = CutTree(node->left, cuttingPoint, i + 1);
			if (tryCut != nullptr)
				return tryCut;
			else
				return CutTree(node->right, cuttingPoint, i + 1);				
		}
		else {				
			TreeNode* tryCut = CutTree(node->right, cuttingPoint, i + 1);
			if (tryCut != nullptr)
				return tryCut;
			else
				return CutTree(node->left, cuttingPoint, i + 1);
		}
	}
	else
		return nullptr; //Reached end of tree without reaching cuttingPoint

}

TreeNode * CutTree(
	TreeNode* node,
	unsigned int cuttingPoint) {
	return CutTree(node, cuttingPoint, 1);
}

unsigned int Height(TreeNode* node);
unsigned int Height(TreeNode* node) {
	if (node->left != nullptr && node->right != nullptr)
		return std::max(1 + Height(node->left), 1 + Height(node->right));
	else
		return 1;
}

void CrossoverTree(
	TreeNode* left,
	TreeNode* right) {
	unsigned int leftHeight = Height(left);
	unsigned int rightHeight = Height(right);
	unsigned int maxHeight = 5;
	unsigned int rightCuttingPoint, leftCuttingPoint, min, max;

	//Keep cutting points between points that don't make the tree grow beyon maxHeight
	if (leftHeight > rightHeight) {
		leftCuttingPoint =  rand() % (leftHeight - 1) + 1;
		min = rightHeight - std::min(rightHeight, maxHeight - leftCuttingPoint);
		max = std::min(rightHeight - 1, maxHeight - (leftHeight - leftCuttingPoint));
		rightCuttingPoint = (rand() % (max - min + 1)) + min;
	}
	else {		
		rightCuttingPoint = rand() % (rightHeight - 1) + 1;
		min = leftHeight - std::min(leftHeight, maxHeight - rightCuttingPoint);
		max = std::min(leftHeight - 1, maxHeight - (rightHeight - rightCuttingPoint));
		leftCuttingPoint = (rand() % (max - min + 1)) + min;		
	}
	
	TreeNode* leftCut = CutTree(left, leftCuttingPoint);
	TreeNode* rightCut = CutTree(right, rightCuttingPoint);
	if (rand() % 2 == 0) {
		TreeNode* copy = leftCut->left->Copy();
		delete leftCut->left;
		if (rand() % 2 == 0) {
			leftCut->left = rightCut->right->Copy();
			delete rightCut->right;
			rightCut->right = copy;
		}
		else {
			leftCut->left = rightCut->left->Copy();			
			delete rightCut->left;
			rightCut->left = copy;
		}
	}
	else {			
		TreeNode* copy = leftCut->right->Copy();
		delete leftCut->right;
		if (rand() % 2 == 0) {
			leftCut->right = rightCut->right->Copy();		
			delete rightCut->right;
			rightCut->right = copy;
		}
		else {
			leftCut->right = rightCut->left->Copy();		
			delete rightCut->left;
			rightCut->left = copy;
		}
	}
}

void BinaryTreeGeneticAlgo::Crossover(
	BinaryTreeChromosome * left,
	BinaryTreeChromosome * right) {
	CrossoverTree(left->buy, right->buy);
	CrossoverTree(left->sell, right->sell);
}
