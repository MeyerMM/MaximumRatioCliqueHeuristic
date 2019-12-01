#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <math.h>

struct candidateNode {
	int nodeId;
	float potentialObjectiveF;
};

bool candidateComparator(candidateNode n1, candidateNode n2) { 
	return (n1.potentialObjectiveF > n2.potentialObjectiveF); // Sort results in ascending order but we want descending order. Because of that the comparator is inverted
}


void printArray(int size, int* array) {
	printf("\n");
	for (int i = 0; i < size; i++) {
		printf("%i ", array[i]);
	}
	printf("\n");
}

void printAdyList(int numNodes, std::vector<std::vector<int>> adyList) {
	printf("ADY LIST:\n");
	for (int i = 0; i < numNodes; i++) {
		std::cout << i << ": ";
		for (std::vector<int>::iterator it = adyList[i].begin(); it != adyList[i].end(); ++it)
			std::cout << ' ' << *it;
		std::cout << '\n';
	}
	printf("ADY LIST END\n");
}

void initData(int numNodes, int maxWeight, float density, int* pWeights, int* qWeights, int* adyMat) {
	int baseIndex;
	for (int i = 0; i < numNodes; i++) {
		pWeights[i] = rand() % maxWeight;
		qWeights[i] = rand() % maxWeight;
		baseIndex = i*numNodes;
		for (int j =  i + 1; j < numNodes; j++) {
			float value = rand() / (float)RAND_MAX;
			if (value < density) {
				adyMat[baseIndex + j] = 1;
				adyMat[j * numNodes + i] = 1;
			}
			else {
				adyMat[baseIndex + j] = 0;
				adyMat[j * numNodes + i] = 0;
			}
				
		}
	}
}

void createAdyListFromAdyMat(int numNodes, int* adyMat, std::vector<std::vector<int>>* adyList) {
	int baseIndex;
	for (int i = 0; i < numNodes; i++) {
		baseIndex = i*numNodes;
		for (int j = 0; j < numNodes; j++) {
			if (adyMat[baseIndex + j]) {
				adyList->at(i).push_back(j);
			}
		}
		sort(adyList->at(i).begin(), adyList->at(i).end());
	}
}

void printGraph(int numNodes, int maxWeight, float density, int* pWeights, int* qWeights, int* adyMat) {
	printf("numNodes: %i, maxWeight: %i, density: %g \n", numNodes, maxWeight, density);
	printf("Node\tWeight P\tWeight Q\n");
	for (int i = 0; i < numNodes; i++) {
		printf("%i\t%i\t%i \n", i, pWeights[i], qWeights[i]);
	}
	printf("Adyacency Matrix\n");
	int baseIndex;
	for (int i = 0; i < numNodes; i++) {
		baseIndex = i*numNodes;
		for (int j = 0; j < numNodes; j++) {
			printf("%i\t", adyMat[baseIndex + j]);
		}
		printf("\n");
	}
}

int* greedyConstructor(int numNodes, int* pWeights, int* qWeights, int* adyMat, std::vector<std::vector<int>> adyList) {
	// Returns pointer to int array where the 0 position indicates the size of the array

	int bestNode = -1;
	float bestRatio = 0; // Select the first node
	for (int i = 0; i < numNodes; i++) {
		float newRatio = (float)(pWeights[i]) / (float)(qWeights[i]);
		if (newRatio > bestRatio) {
			bestNode = i;
			bestRatio = newRatio;
		}
	}

	int* solution = (int*)malloc(numNodes * sizeof(int)); // The biggest clique possible is the whole graph
	solution[0] = bestNode;
	int nodesInClique = 1;
	std::vector<int> potentialNodes = adyList[bestNode];
	int pAcum = pWeights[bestNode];
	int qAcum = qWeights[bestNode];
	std::vector<int> curr_intersection;
	float newRatio = 0.0;
	while (!potentialNodes.empty()) {
		bestNode = -1;
		bestRatio = 0;  // It must start at 0 because the cliques must be maximal
		for (std::vector<int>::iterator it = potentialNodes.begin(); it != potentialNodes.end(); ++it) {
			int candidateNode = *it;
			newRatio = (float)(pAcum + pWeights[candidateNode])/ (float)(qAcum + qWeights[candidateNode]);
			if (newRatio > bestRatio) {
				bestNode = candidateNode;
				bestRatio = newRatio;
			}
		}

		solution[nodesInClique] = bestNode;
		nodesInClique++;
		pAcum += pWeights[bestNode];
		qAcum += qWeights[bestNode];

		
		std::set_intersection(potentialNodes.begin(), potentialNodes.end(), adyList.at(bestNode).begin(), adyList.at(bestNode).end(), std::back_inserter(curr_intersection));
		std::swap(potentialNodes, curr_intersection);
		curr_intersection.clear();
	}

	// Copy solution to an array with only the necessary space 
	int* finalSolution = (int*)malloc((nodesInClique + 1) * sizeof(int)); // Allocate space so the first position of the array stores its own size
	finalSolution[0] = nodesInClique;
	memcpy(&finalSolution[1], solution, nodesInClique * sizeof(int));
	free(solution);
	return finalSolution;
}
/*
int* GRASP(int numNodes, int* pWeights, int* qWeights, int* adyMat, std::vector<std::vector<int>> adyList) {
	// Returns pointer to int array where the 0 position indicates the size of the array
	float alfa = 0.25;
	int candidateNumber = numNodes *alfa;
	if (candidateNumber < 1) {
		candidateNumber = 1;
	}
	int* candidateList = (int*)malloc(candidateNumber * sizeof(int));
	float* candidateListOF = (float*)malloc(candidateNumber * sizeof(float));

	int worstCandidate = -1;
	float worstCandidateRatio = 999999;


	// Find the first node for the solution
	// Fill the candidate list
	for (int i = 0; i < candidateNumber; i++) {
		float newRatio = (float)(pWeights[i]) / (float)(qWeights[i]);
		candidateList[i] = i;
		candidateListOF[i] = newRatio;
		if (newRatio < worstCandidateRatio) {
			worstCandidate = i;
			worstCandidateRatio = newRatio;
		}
	}
	// Continue evaluating nodes
	for (int i = candidateNumber; i < numNodes; i++) {
		float newRatio = (float)(pWeights[i]) / (float)(qWeights[i]);
		if (newRatio > worstCandidateRatio) {
			candidateList[worstCandidate] = i;
			candidateListOF[worstCandidate] = newRatio;
			if (i + 1 < numNodes) { // Find the new worst candidate, but only if the node evaluation has nodes left
				worstCandidate = 0;
				worstCandidateRatio = candidateListOF[0];
				for (int j = 1; j < candidateNumber; j++) {
					if (candidateListOF[j] < worstCandidateRatio) {
						worstCandidate = j;
						worstCandidateRatio = candidateListOF[j];
					}
				}
			}
		}
	}

	// Choose the first node randomly from the candidate list
	int selectedIndex = rand() % candidateNumber;
	printf("%i \n", selectedIndex);
	int selectedCandidate = candidateList[selectedIndex];
	free(candidateList);
	free(candidateListOF);


	int* solution = (int*)malloc(numNodes * sizeof(int)); // The biggest clique possible is the whole graph
	solution[0] = selectedCandidate;
	int nodesInClique = 1;
	std::vector<int> potentialNodes = adyList[selectedCandidate];
	int pAcum = pWeights[selectedCandidate];
	int qAcum = qWeights[selectedCandidate];
	std::vector<int> curr_intersection;
	float newRatio = 0.0;
	while (!potentialNodes.empty()) {
		candidateNumber = (float)potentialNodes.size() *alfa;
		if (candidateNumber < 1) {
			candidateNumber = 1;
		}

		int* candidateList = (int*)malloc(candidateNumber * sizeof(int));
		float* candidateListOF = (float*)malloc(candidateNumber * sizeof(float));
		worstCandidate = -1;
		worstCandidateRatio = 999999;

		for (int i = 0; i < candidateNumber; i++) {
			int candidateNode = potentialNodes[i];
			newRatio = (float)(pAcum + pWeights[candidateNode]) / (float)(qAcum + qWeights[candidateNode]);
			candidateList[i] = candidateNode;
			candidateListOF[i] = newRatio;
			if (newRatio < worstCandidateRatio) {
				worstCandidate = i;
				worstCandidateRatio = newRatio;
			}
		}
		for (unsigned int i = candidateNumber; i < potentialNodes.size(); i++) {
			int candidateNode = potentialNodes[i];
			newRatio = (float)(pAcum + pWeights[candidateNode]) / (float)(qAcum + qWeights[candidateNode]);
			if (newRatio > worstCandidateRatio) {
				candidateList[worstCandidate] = candidateNode;
				candidateListOF[worstCandidate] = newRatio;
				if (i + 1 < numNodes) { // Find the new worst candidate, but only if the node evaluation has nodes left
					worstCandidate = 0;
					worstCandidateRatio = candidateListOF[0];
					for (int j = 1; j < candidateNumber; j++) {
						if (candidateListOF[j] < worstCandidateRatio) {
							worstCandidate = j;
							worstCandidateRatio = candidateListOF[j];
						}
					}
				}
			}

		}


		selectedIndex = 0;
		if (candidateNumber > 1) {
			selectedIndex = (int)(rand() % candidateNumber);
		}
		selectedCandidate = candidateList[selectedIndex];


		solution[nodesInClique] = selectedCandidate;
		nodesInClique++;
		pAcum += pWeights[selectedCandidate];
		qAcum += qWeights[selectedCandidate];

		std::set_intersection(potentialNodes.begin(), potentialNodes.end(), adyList.at(selectedCandidate).begin(), adyList.at(selectedCandidate).end(), std::back_inserter(curr_intersection));
		std::swap(potentialNodes, curr_intersection);
		curr_intersection.clear();

		std::set_intersection(potentialNodes.begin(), potentialNodes.end(), adyList.at(selectedCandidate).begin(), adyList.at(selectedCandidate).end(), std::back_inserter(curr_intersection));
		std::swap(potentialNodes, curr_intersection);
		curr_intersection.clear();
		free(candidateList);
		free(candidateListOF);
	}

	// Copy solution to an array with only the necessary space 
	int* finalSolution = (int*)malloc((nodesInClique + 1) * sizeof(int)); // Allocate space so the first position of the array stores its own size
	finalSolution[0] = nodesInClique;
	memcpy(&finalSolution[1], solution, nodesInClique * sizeof(int));
	free(solution);
	return finalSolution;
}
*/
int* GRASPV2(int numNodes, int* pWeights, int* qWeights, int* adyMat, std::vector<std::vector<int>> adyList) {
	// Returns pointer to int array where the 0 position indicates the size of the array
	float alfa = 0.25f;
	int* solution = (int*)malloc(numNodes * sizeof(int)); // The biggest clique possible is the whole graph
	std::vector<candidateNode> candidateList(numNodes);  // The biggest candidate list possible is the whole graph
	int candidateNumber = (int)(numNodes *alfa + 0.5f);  // How many canidates to consider
	if (candidateNumber < 1) {							// At least one candidate
		candidateNumber = 1;
	}

	// Find the first node for the solution
	for (int i = 0; i < numNodes; i++) { // Evaluate every node in the graph
		candidateNode cN;
		cN.nodeId = i;
		cN.potentialObjectiveF = (float)(pWeights[i]) / (float)(qWeights[i]);
		candidateList.push_back(cN);
	}
	std::sort(candidateList.begin(), candidateList.end(), candidateComparator);  // Sort so the best nodes are at the beggining

	// Choose the first node randomly from the candidate list
	int selectedIndex = rand() % candidateNumber;
	int selectedCandidate = candidateList[selectedIndex].nodeId;

	// Add the node to the solution
	solution[0] = selectedCandidate;
	int nodesInClique = 1;
	std::vector<int> potentialNodes = adyList[selectedCandidate];  // Potential nodes for the clique
	int pAcum = pWeights[selectedCandidate];
	int qAcum = qWeights[selectedCandidate];
	std::vector<int> curr_intersection;  // Auxiliar vector for intersections

	// Complete the solution
	while (!potentialNodes.empty()) {
		// Empty the candidate list from the previous iteration
		candidateList.clear();

		candidateNumber = (int)((float)potentialNodes.size() *alfa + 0.5f); // Reevaluate the number of candidates
		if (candidateNumber < 1) {	// At least one candidate
			candidateNumber = 1;
		}

		for (unsigned int i = 0; i < potentialNodes.size(); i++) { // Evaluate every potential node
			candidateNode cN;
			cN.nodeId = potentialNodes[i];
			cN.potentialObjectiveF = (float)(pAcum + pWeights[cN.nodeId]) / (float)(qAcum + qWeights[cN.nodeId]);
			candidateList.push_back(cN);
		}

		std::sort(candidateList.begin(), candidateList.end(), candidateComparator); // Sort the nodes by descending objective function
		
		selectedIndex = (int)(rand() % candidateNumber);  // Randomly select one of the candidates
		selectedCandidate = candidateList[selectedIndex].nodeId;
		// Add node to solution
		solution[nodesInClique] = selectedCandidate;
		nodesInClique++;
		pAcum += pWeights[selectedCandidate];
		qAcum += qWeights[selectedCandidate];
		// Find the new list of potential nodes for the clique
		std::set_intersection(potentialNodes.begin(), potentialNodes.end(), adyList.at(selectedCandidate).begin(), adyList.at(selectedCandidate).end(), std::back_inserter(curr_intersection));
		std::swap(potentialNodes, curr_intersection);
		curr_intersection.clear();
	}

	// Copy solution to an array with only the necessary space 
	int* finalSolution = (int*)malloc((nodesInClique + 1) * sizeof(int)); // Allocate space so the first position of the array stores its own size
	finalSolution[0] = nodesInClique;
	memcpy(&finalSolution[1], solution, nodesInClique * sizeof(int));
	free(solution);
	return finalSolution;
}

int* generateSolution(int numNodes, int* adyMat) {
	// Returns pointer to int array where the 0 position indicates the size of the array
	int nodesInClique = 1;
	int* solution = (int*)malloc(numNodes * sizeof(int)); // The biggest clique possible is the whole graph
	
	solution[0] = rand() % numNodes; // Start with a random node for a random solution
	int baseIndex;
	for (int i = 0; i < numNodes; i++) { // Check every node in the graph
		baseIndex = i*numNodes;
		bool enterClique = true;
		for (int j = 0; j < nodesInClique; j++) { // And every node in the clique
			if (i != solution[j] && adyMat[baseIndex + solution[j]]) {  // If the two nodes are not the same and they are connected
				enterClique = enterClique && true;  // The new node is a pontential new member for the clique
			}
			else {
				enterClique = enterClique && false; // Else, discard the new node
				break;
			}
		}
		if (enterClique) {
			solution[nodesInClique] = i;
			nodesInClique++;
		}
	}

	// Copy solution to an array with only the necessary space 
	int* finalSolution = (int*)malloc((nodesInClique + 1) * sizeof(int)); // Allocate space so the first position of the array stores its own size
	finalSolution[0] = nodesInClique;
	memcpy(&finalSolution[1], solution, nodesInClique * sizeof(int));
	free(solution);
	return finalSolution;
}

bool validateCliqueCompleteness(int numNodes, int* adyMat, int nodesInClique, int* clique) {
	int baseIndex;
	for (int i = 0; i < nodesInClique; i++) { // For every node in the clique
		baseIndex = clique[i] * numNodes;
		for (int j = i + 1; j < nodesInClique; j++) { // Check if every other node is connected to it
			if (!adyMat[baseIndex + clique[j]]) { // If not, they do not form a clique
				return false;
			}
		}
	}
	return true; // When every node in a subgraph is connected to every other node, the subgraph is a clique
}

bool validateCliqueMaximality(int numNodes, std::vector<std::vector<int>> adyList, int nodesInClique, int* clique) {
	// Check if there are any nodes in the graph that could join the clique
	/*
	int baseIndex;
	for (int i = 0; i < numNodes; i++) {
		baseIndex = i*numNodes;
		bool enterClique = true;
		for (int j = 0; j < nodesInClique; j++) {
			if (i != clique[j] && adyMat[baseIndex + clique[j]]) {
				enterClique = enterClique && true;
			}
			else {
				enterClique = enterClique && false;
				break;
			}
		}
		if (enterClique) {
			return false; // If a node can enter the clique, then it is not maximal
		}
	}
	return true; // If no node can enter the clique, then it is maximal
	*/
	std::vector<int> last_intersection = adyList.at(clique[0]);
	std::vector<int> curr_intersection;
	for (int i = 1; i < nodesInClique; i++) {
		std::set_intersection(last_intersection.begin(), last_intersection.end(), adyList.at(clique[i]).begin(), adyList.at(clique[i]).end(), std::back_inserter(curr_intersection));
		std::swap(last_intersection, curr_intersection);
		curr_intersection.clear();
	}
	if (last_intersection.empty()) {
		
		return true;
	}
	else {
		return false;
	}
}

bool validateSolution(int numNodes, int* adyMat, std::vector<std::vector<int>> adyList, int solSize, int* solution) {
	// A solution is valid if it is a maximal clique of the graph
	bool isValidSolution = validateCliqueCompleteness(numNodes, adyMat, solSize, solution); // Validate if it is a clique
	isValidSolution = isValidSolution && validateCliqueMaximality(numNodes, adyList, solSize, solution); // Validate if it is maximal
	return isValidSolution;
}

float calculateObjectiveF(int solSize, int* solution, int* pWeights, int* qWeights) {
	// Calculate the relation between the sums of the weights.
	int pAcum = 0;
	int qAcum = 0;
	for (int i = 0; i < solSize; i++) {
		pAcum += pWeights[solution[i]];
		qAcum += qWeights[solution[i]];
	}
	return (float)pAcum / (float)qAcum;
}


int main(int argc, char *argv[]) {
	int numNodes = 100;
	int maxWeight = 500;
	float density = 0.8f;
	bool verifyResults = true;
	srand(1);

	size_t nodeListSize = numNodes * sizeof(int);
	int* pWeights = (int*)malloc(nodeListSize);
	int* qWeights = (int*)malloc(nodeListSize);
	int* adyMat = (int*)calloc(numNodes * numNodes, sizeof(int));
	std::vector<std::vector<int>> adyList(numNodes);

	initData(numNodes, maxWeight, density, pWeights, qWeights, adyMat);
	createAdyListFromAdyMat(numNodes, adyMat, &adyList);
	printf("\nGraph initialization finished\n");
	//printGraph(numNodes, maxWeight, density, pWeights, qWeights, adyMat);
	//printAdyList(numNodes, adyList);

	int* solution = generateSolution(numNodes, adyMat);
	printf("\nRandom Solution: %i nodes", solution[0]);
	printArray(solution[0], &solution[1]);
	printf("Objective Function: %g\n", calculateObjectiveF(solution[0], solution + 1, pWeights, qWeights));
	printf("Valid Solution: %i\n", validateSolution(numNodes, adyMat, adyList, solution[0], solution + 1));

	solution = greedyConstructor(numNodes, pWeights, qWeights, adyMat, adyList);
	printf("\nGreedy Solution: %i nodes", solution[0]);
	printArray(solution[0], &solution[1]);
	printf("Objective Function: %g\n", calculateObjectiveF(solution[0], solution + 1, pWeights, qWeights));
	printf("Valid Solution: %i\n", validateSolution(numNodes, adyMat, adyList, solution[0], solution + 1));

	solution = GRASPV2(numNodes, pWeights, qWeights, adyMat, adyList);
	printf("\nGRASP V2 Solution: %i nodes", solution[0]);
	printArray(solution[0], &solution[1]);
	printf("Objective Function: %g\n", calculateObjectiveF(solution[0], solution + 1, pWeights, qWeights));
	printf("Valid Solution: %i\n", validateSolution(numNodes, adyMat, adyList, solution[0], solution + 1));


	free(pWeights);
	free(qWeights);
	free(adyMat);
	free(solution);
	getchar();
	return(0);
}