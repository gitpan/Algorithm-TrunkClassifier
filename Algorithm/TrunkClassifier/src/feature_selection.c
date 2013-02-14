/*
 * MaxClassPol.cpp
 *
 *  Created on: 26 feb 2012
 *      Author: Wolfcastle
 *
 *  Description
 *		This library implements the independent t-test and
 *		Maximum Class Polarization (MCP) feature selection methods
 */

#include "feature_selection.h"

/*
 * Description: Support function that checks if a sample belongs to a given class
 */
int checkClass(char* sampleName, char* group){
	unsigned int sameCount;
	unsigned int samplePos;
	unsigned int groupPos;
	for(samplePos = 0; samplePos < strlen(sampleName) - (strlen(group) - 1); samplePos++){
		sameCount = 0;
		for(groupPos = 0; groupPos < strlen(group); groupPos++){
			if(sampleName[samplePos + groupPos] == group[groupPos])
				sameCount++;
		}
		if(sameCount == strlen(group)){
			return 1;
		}
	}
	return 0;
}

/*
Description: Feature selection method that uses the independent t-test to select a feature
Parameters: (1) 2D expression data matrix, (2) number of matrix rows, (3) number of matrix columns (4) list of column names,
            (5) first group symbol and (6) second group symbol
Return value: Index of the top t-value feature (row in data matrix)
*/
int indTTest(double** expData, int numFeatures, int numSamples, char** sampleNames, char* NORMAL, char* MALIGN){
	
	//Determine class sizes
	int sample;
	int numNormal = 0;
	int numMalign = 0;
	for(sample = 0; sample < numSamples; sample++){
		if(checkClass(sampleNames[sample], NORMAL))
			numNormal++;
		else if(checkClass(sampleNames[sample], MALIGN))
			numMalign++;
	}
	
	//Establish class index arrays
	int normalIndexes[numNormal];
	int malignIndexes[numMalign];
	int normalCounter = 0;
	int malignCounter = 0;
	for(sample = 0; sample < numSamples; sample++){
		if(checkClass(sampleNames[sample], NORMAL)){
			normalIndexes[normalCounter] = sample;
			normalCounter++;
		}
		else if(checkClass(sampleNames[sample], MALIGN)){
			malignIndexes[malignCounter] = sample;
			malignCounter++;
		}
	}
	
	//Perform t-test
	int feature;
	double currentT = 0;
	double topT = 0;
	int bestFeature;
	for(feature = 0; feature < numFeatures; feature++){

		//Calcualte sum of normal samples
		double sumNormal = 0;
		int normalBuffer;
		for(sample = 0; sample < numNormal; sample++){
			normalBuffer = normalIndexes[sample];
			sumNormal += expData[feature][normalBuffer];
		}
		//Calculate sum of malign samples
		double sumMalign = 0;
		int malignBuffer;
		for(sample = 0; sample < numMalign; sample++){
			malignBuffer = malignIndexes[sample];
			sumMalign += expData[feature][malignBuffer];
		}

		//Calculate mean
		double normalMean = sumNormal / numNormal;
		double malignMean = sumMalign / numMalign;

		//Calculate variance of normal samples
		double normalVar = 0;
		for(sample = 0; sample < numNormal; sample++){
			normalBuffer = normalIndexes[sample];
			normalVar += pow(expData[feature][normalBuffer] - normalMean, 2);
		}

		//Calculate variance of malign samples
		double malignVar = 0;
		for(sample = 0; sample < numMalign; sample++){
			malignBuffer = malignIndexes[sample];
			malignVar += pow(expData[feature][malignBuffer] - malignMean, 2);
		}
		
		//Calculate t and compare to previous best
		if(normalVar == 0 && malignVar == 0){
			normalVar = 0.000001;
		}
		currentT = fabs(normalMean - malignMean) / sqrt ((normalVar / numNormal) + (malignVar / numMalign));
		if(currentT > topT){
			topT = currentT;
			bestFeature = feature;
		}
	}
	
	//Return index of best feature
	return bestFeature;
}
