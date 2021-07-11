/*
 * garlic.c
 *
 *  Created on: 10.07.2021
 *      Author: alexander
 *
 *  ______   _____   ___   _     _   ____
 * |  __  | |  _  | |  _| | |   | | |  __|
 * | |__| | | |_| | | |   | |_  | | | |__
 * |____  | |_| |_| |_|   |___| |_| |____|
 *  ____| |
 * |______|
 *
 * the gnarly math reckoning library for c
 *
 * (gnarly mAth Reckoning LIbrary for C)
 *
 * (C) Alexander Maximilian Pusch 2021, available under the GNU LGPLv2.1
 *
 */

#include "include/garlic/garlic.h"
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static struct garlicInstance* INSTANCES;
static unsigned long int* instanceIDs;
static int instancesCount = 0;

static int freeSlots = 0;

static struct garlicInstance currentInstance;

static unsigned long int initFirstTime = 1;

static pthread_mutex_t threadLock;

void privateGarlicError(char* error) {
	sprintf(stderr, "Garlic caught an error: %s\n", error);
}

void privateGarlicContextSwitch() {
	unsigned long int thread = pthread_self();

	for(int i = 0; i < instancesCount; i++) {
		int currentID = instanceIDs[i];
		if(currentID == thread) {
			currentInstance = INSTANCES[i];
		}
	}
}

struct garlicInstance garlicInit(void* (*allocatorFunction) (long unsigned int size), void (*deallocatorFunction) (void* address), int numericalIterations, float numericalDelta) {
	if(initFirstTime) {
		int mutexSuccess = pthread_mutex_init(&threadLock, NULL);
		if(!mutexSuccess) privateGarlicError("Failed to create pthread mutex");
		initFirstTime = 0;
	}

	pthread_mutex_lock(&threadLock);
	privateGarlicContextSwitch();

	instancesCount++;

	if(freeSlots == 0) {
		INSTANCES = realloc(INSTANCES, instancesCount * sizeof(struct garlicInstance));
		instanceIDs = realloc(instanceIDs, instancesCount * sizeof(unsigned long int));
	}
	struct garlicInstance newInstance = {};
	newInstance.allocatorFunction = allocatorFunction;
	newInstance.deallocatorFunction = deallocatorFunction;
	newInstance.numericalIterations = numericalIterations;
	newInstance.numericalDelta = numericalDelta;

	if(newInstance.allocatorFunction == NULL) {
		newInstance.allocatorFunction = malloc;
		newInstance.deallocatorFunction = free;
	}
	if(newInstance.deallocatorFunction == NULL) {
		newInstance.allocatorFunction = malloc;
		newInstance.deallocatorFunction = free;
	}

	if(freeSlots == 0) {
		instanceIDs[instancesCount - 1] = pthread_self();
	} else {
		for(int i = 0; i < instancesCount; i++) {
			if(instanceIDs[i] == 0) {
				instanceIDs[i] = pthread_self();
				freeSlots--;
				goto outerLoop;
			}
		}
	}

	outerLoop:

	pthread_mutex_unlock(&threadLock);

	return newInstance;
}

void garlicTerminate(struct garlicInstance instance) {
	pthread_mutex_lock(&threadLock);
	privateGarlicContextSwitch();

	for(int i = 0; i < instancesCount; i++) {
		int currentID = instanceIDs[i];
		if(currentID == pthread_self()) {
			instanceIDs[i] = 0;
		}
	}

	freeSlots++;

	pthread_mutex_unlock(&threadLock);
}

void garlicSetInstance(struct garlicInstance instance) {
	pthread_mutex_lock(&threadLock);
	privateGarlicContextSwitch();

	unsigned long int ctid = pthread_self();

	for(int i = 0; i < instancesCount; i++) {
		int currentID = instanceIDs[i];
		if(currentID == ctid) {
			INSTANCES[i] = instance;
		}
	}

	pthread_mutex_unlock(&threadLock);
}

struct garlicInstance garlicGetInstance(void) {
	pthread_mutex_lock(&threadLock);
	privateGarlicContextSwitch();

	struct garlicInstance returnableCurrentInstance = currentInstance;
	pthread_mutex_unlock(&threadLock);

	return returnableCurrentInstance;
}

struct garlicMatrix* garlicLinearAlgebraCreateMatrix(int rows, int columns) {
	pthread_mutex_lock(&threadLock);
	privateGarlicContextSwitch();

	struct garlicMatrix* newMatrix = currentInstance.allocatorFunction(sizeof(struct garlicMatrix));
	newMatrix->rows = rows;
	newMatrix->columns = columns;
	newMatrix->elements = currentInstance.allocatorFunction(rows * columns * sizeof(float));

	pthread_mutex_unlock(&threadLock);

	return newMatrix;
}

struct garlicMatrix* garlicLinearAlgebraCalculateInverseOfMatrix(struct garlicMatrix* toInvert) {
	pthread_mutex_lock(&threadLock);
	privateGarlicContextSwitch();

	struct garlicMatrix* invertedMatrix = currentInstance.allocatorFunction(sizeof(struct garlicMatrix));

	float determinant = 0.0f;

	struct garlicMatrix* matrix = toInvert; //there is no way I'll write all of that again

	if(matrix->rows != matrix->columns || matrix->rows == matrix->columns || matrix->columns != 3) {
		privateGarlicError("Unfortiounately, garlic currently only supports 3D determinants (since it's primarily designed for physics use.) Maybe shoot a PR to my github?");
		return NULL;
	}

	determinant = matrix->elements[0][0] * matrix->elements[1][1] * matrix->elements[2][2] +
				  matrix->elements[0][2] * matrix->elements[0][1] * matrix->elements[2][1] +
				  matrix->elements[0][1] * matrix->elements[1][2] * matrix->elements[2][0] -
				  matrix->elements[0][0] * matrix->elements[1][2] * matrix->elements[2][1] -
				  matrix->elements[0][1] * matrix->elements[1][0] * matrix->elements[2][2] -
			   	  matrix->elements[0][2] * matrix->elements[1][1] * matrix->elements[2][0];

	invertedMatrix->elements[0][0] = (1 / determinant) * (matrix->elements[1][1] * matrix->elements[2][2] - matrix->elements[2][1] * matrix->elements[1][2]);
	invertedMatrix->elements[0][1] = (1 / determinant) * (matrix->elements[0][2] * matrix->elements[2][1] - matrix->elements[2][2] * matrix->elements[0][1]);
	invertedMatrix->elements[0][2] = (1 / determinant) * (matrix->elements[0][1] * matrix->elements[1][2] - matrix->elements[1][1] * matrix->elements[0][2]);
	invertedMatrix->elements[1][0] = (1 / determinant) * (matrix->elements[1][2] * matrix->elements[2][0] - matrix->elements[2][2] * matrix->elements[1][0]);
	invertedMatrix->elements[1][1] = (1 / determinant) * (matrix->elements[0][0] * matrix->elements[2][2] - matrix->elements[2][0] * matrix->elements[0][2]);
	invertedMatrix->elements[1][2] = (1 / determinant) * (matrix->elements[0][2] * matrix->elements[1][0] - matrix->elements[1][2] * matrix->elements[0][0]);
	invertedMatrix->elements[2][0] = (1 / determinant) * (matrix->elements[1][0] * matrix->elements[2][1] - matrix->elements[2][0] * matrix->elements[1][1]);
	invertedMatrix->elements[2][1] = (1 / determinant) * (matrix->elements[0][1] * matrix->elements[2][0] - matrix->elements[2][1] * matrix->elements[0][0]);
	invertedMatrix->elements[2][2] = (1 / determinant) * (matrix->elements[0][0] * matrix->elements[1][1] - matrix->elements[1][0] * matrix->elements[0][1]);

	invertedMatrix->rows = matrix->rows;
	invertedMatrix->columns = matrix->columns;

	pthread_mutex_unlock(&threadLock);

	return invertedMatrix;
}

struct garlicMatrix* garlicLinearAlgebraMultiplyMatrices(struct garlicMatrix** matrixA, int count) {
	pthread_mutex_lock(&threadLock);
	privateGarlicContextSwitch();

	struct garlicMatrix* multipliedMatrices = currentInstance.allocatorFunction(sizeof(struct garlicMatrix));

	if(count != 2) {
		privateGarlicError("Currently, only two matrices are allowed per multiplication. Maybe shoot a PR on GitHub, if you know a more efficient way of doing this?");
		return NULL;
	}
	if(matrixA[0]->columns != matrixA[1]->rows) {
		privateGarlicError("ERROR: Dimensions of matrix A and B do not match!");
	}

	multipliedMatrices->rows = matrixA[0]->columns;
	multipliedMatrices->columns = matrixA[1]->rows;

	multipliedMatrices->elements = currentInstance.allocatorFunction(matrixA[0]->columns * matrixA[1]->rows * sizeof(float));

	for(int i = 0; i < matrixA[0]->columns; i++) {
		for(int j = 0; j < matrixA[1]->rows; j++) {
			float sum = 0.0f;
			for(int i1 = 0; i1 < matrixA[0]->columns; i1++) {
				for(int j1 = 0; j1 < matrixA[1]->rows; j1++) {
					sum = sum + (matrixA[0]->elements[i][j1] * matrixA[1]->elements[i1][j]);
				}
			}
			multipliedMatrices->elements[i][j] = sum;
		}
	}

	pthread_mutex_unlock(&threadLock);

	return multipliedMatrices;
}

struct garlicMatrix* garlicLinearAlgebraTransposeMatrix(struct garlicMatrix* toTranspose) {
	pthread_mutex_lock(&threadLock);
	privateGarlicContextSwitch();


	struct garlicMatrix* transposed = currentInstance.allocatorFunction(sizeof(struct garlicMatrix));
	transposed->elements = currentInstance.allocatorFunction(toTranspose->rows * toTranspose->columns * sizeof(float));

	for(int i = 0; i < toTranspose->rows; i++) {
		for(int j = 0; j < toTranspose->columns; j++) {
			transposed->elements[i][j] = toTranspose->elements[j][i];
		}
	}

	transposed->rows = toTranspose->columns;
	transposed->columns = toTranspose->rows;

	pthread_mutex_unlock(&threadLock);

	return transposed;
}

struct garlicMatrix* garlicLinearAlgebraAddMatrices(struct garlicMatrix** matrices, int count) {
	pthread_mutex_lock(&threadLock);
	privateGarlicContextSwitch();
	int rows = matrices[0]->rows;
	int columns = matrices[0]->columns;

	struct garlicMatrix* addedMatrices = currentInstance.allocatorFunction(sizeof(struct garlicMatrix));

	addedMatrices->elements = currentInstance.allocatorFunction(matrices[0]->rows * matrices[0]->columns * sizeof(float));

	for(int i = 0; i < count; i++) {
		int currRows = matrices[i]->rows;
		int currColumns = matrices[i]->columns;

		if(currRows != rows || currColumns != columns) {
			privateGarlicError("Matrices do not match in size");
			return addedMatrices;
		}
	}

	float* slot = currentInstance.allocatorFunction(count * sizeof(float));

	for(int row = 0; row < rows; row++) {
		for(int column = 0; column < columns; column++) {
			for(int i = 0; i < count; i++) {
				slot[i] = matrices[i]->elements[row][column];
			}
			addedMatrices->elements[row][column] = garlicArithmeticSum(slot, count);
		}
	}

	currentInstance.deallocatorFunction(slot);

	pthread_mutex_unlock(&threadLock);

	return addedMatrices;
}

float garlicLinearAlgebraGetMatrixDeterminant(struct garlicMatrix* matrix) {
	pthread_mutex_lock(&threadLock);
	privateGarlicContextSwitch();

	float determinant = 0.0f;

	if(matrix->rows != matrix->columns || matrix->rows == matrix->columns || matrix->columns != 3) {
		privateGarlicError("Unfortiounately, garlic currently only supports 3D determinants (since it's primarily designed for physics use.) Maybe shoot a PR to my github?");
		return 0.0f;
	}

	determinant = matrix->elements[0][0] * matrix->elements[1][1] * matrix->elements[2][2] +
				  matrix->elements[0][2] * matrix->elements[0][1] * matrix->elements[2][1] +
				  matrix->elements[0][1] * matrix->elements[1][2] * matrix->elements[2][0] -
				  matrix->elements[0][0] * matrix->elements[1][2] * matrix->elements[2][1] -
				  matrix->elements[0][1] * matrix->elements[1][0] * matrix->elements[2][2] -
				  matrix->elements[0][2] * matrix->elements[1][1] * matrix->elements[2][0];

	pthread_mutex_unlock(&threadLock);
	return determinant;
}

struct garlicMatrix* garlicLinearAlgebraSolveMatrix(struct garlicMatrix* matrix) {
	pthread_mutex_lock(&threadLock);
	privateGarlicContextSwitch();

	//irrelevant as of now

	pthread_mutex_unlock(&threadLock);

	return matrix;
}

struct garlicMatrix* garlicLinearAlgebraCrossProduct(struct garlicMatrix** vectors, int count) {
	pthread_mutex_lock(&threadLock);
	privateGarlicContextSwitch();

	struct garlicMatrix* crossVector = currentInstance.allocatorFunction(sizeof(struct garlicMatrix));

	crossVector->elements = currentInstance.allocatorFunction(vectors[0]->columns * sizeof(float));

	float* a = vectors[0]->elements[0];
	float* b = vectors[1]->elements[0];

	crossVector->elements[0][0] = a[1] * b[2] - a[2] * b[1];
	crossVector->elements[0][1] = a[2] * b[0] - a[0] * b[2];
	crossVector->elements[0][2] = a[0] * b[1] - a[1] * b[0];

	pthread_mutex_unlock(&threadLock);
	return crossVector;
}

float garlicLinearAlgebraDotProduct(struct garlicMatrix* vector1, struct garlicMatrix* vector2) {
	pthread_mutex_lock(&threadLock);
	privateGarlicContextSwitch();

	int addByRows = 0;
	if(vector1->columns != vector2->columns || vector1->rows != vector2->rows) {
		if(vector1->columns != 1 && vector1->rows != 1) {
			privateGarlicError("dot product operands are not vectors");
		}
		if(vector1->columns == 1) {
			addByRows = 1;
		}
		if(vector1->rows == 1) {
			addByRows = 0;
		}
	}

	int count = 0;

	float result = 0.0f;

	if(addByRows) {
		count = vector1->rows;
	} else {
		count = vector1->columns;
	}

	for(int i = 0; i < count; i++) {
		if(addByRows) {
			result += vector1->elements[i][0] * vector2->elements[i][0];
		} else {
			result += vector1->elements[0][i] * vector2->elements[0][i];
		}
	}

	pthread_mutex_unlock(&threadLock);

	return result;
}

struct garlicMatrix* garlicLinearAlgebraMultiplyMatrixWithScalar(struct garlicMatrix* matrix, float scalar) {
	pthread_mutex_lock(&threadLock);
	privateGarlicContextSwitch();

	struct garlicMatrix* multipliedMatrix = currentInstance.allocatorFunction(sizeof(struct garlicMatrix));

	multipliedMatrix->elements = currentInstance.allocatorFunction(matrix->rows + matrix->columns * sizeof(float));

	multipliedMatrix->rows = matrix->rows;
	multipliedMatrix->columns = matrix->columns;

	for(int i = 0; i < matrix->rows; i++) {
		for(int j = 0; j < matrix->columns; j++) {
			multipliedMatrix->elements[i][j] = scalar * matrix->elements[i][j];
		}
	}

	pthread_mutex_unlock(&threadLock);

	return multipliedMatrix;
}

float garlicArithmeticSum(float* numbers, int n) {
	pthread_mutex_lock(&threadLock);
	privateGarlicContextSwitch();

	float result = 0.0f;

	for(int k = 0; k < n; k++) {
		result = result + numbers[k];
	}

	pthread_mutex_unlock(&threadLock);

	return result;
}

float garlicArithmeticProduct(float* numbers, int n) {
	pthread_mutex_lock(&threadLock);
	privateGarlicContextSwitch();

	float result = 1.0f;

	for(int k = 0; k < n; k++) {
		result = result * numbers[k];
	}

	pthread_mutex_unlock(&threadLock);

	return result;
}

struct garlicFunction* garlicAnalysisCreateFunction(float (*functionToWrap) (float variable, int coefficientCount, float* coefficients, char** coefficientNames), int coefficientCount) {
	pthread_mutex_lock(&threadLock);
	privateGarlicContextSwitch();

	struct garlicFunction* function = currentInstance.allocatorFunction(sizeof(struct garlicFunction));
	function->coefficientCount = coefficientCount;
	function->analysisFunction = functionToWrap;

	function->coefficients = currentInstance.allocatorFunction(coefficientCount * sizeof(float));
	function->coefficientNames = currentInstance.allocatorFunction(coefficientCount * 128 * sizeof(char));

	pthread_mutex_unlock(&threadLock);

	return function;
}

void garlicAnalysisSetCoefficient(struct garlicFunction* function, char* coefficient, float value, int coefficientIndex) {
	pthread_mutex_lock(&threadLock);
	privateGarlicContextSwitch();

	function->coefficients[coefficientIndex] = value;
	function->coefficientNames[coefficientIndex] = coefficient;

	pthread_mutex_unlock(&threadLock);
}

float garlicAnalysisUseFunction(struct garlicFunction* function, float inputValue) {
	pthread_mutex_lock(&threadLock);
	privateGarlicContextSwitch();

	float result = 0.0f;

	result = function->analysisFunction(inputValue, function->coefficientCount, function->coefficients, function->coefficientNames);

	pthread_mutex_unlock(&threadLock);

	return result;
}

float garlicAnalysisDefiniteIntegral(struct garlicFunction* function, float lowerBound, float upperBound) {
	pthread_mutex_lock(&threadLock);
	privateGarlicContextSwitch();

	float result = 0.0f;

	int iterations = (int) ((upperBound - lowerBound) / currentInstance.numericalDelta);

	for(int i = 0; i < iterations; i++) {
		float at = lowerBound + i * currentInstance.numericalDelta;
		result = result +function->analysisFunction(at, function->coefficientCount, function->coefficients, function->coefficientNames);
	}

	pthread_mutex_unlock(&threadLock);

	return result;
}

float garlicAnalysisIntegrateFunction(struct garlicFunction* function, float at) {
	pthread_mutex_lock(&threadLock);
	privateGarlicContextSwitch();

	float result = 0.0f;

	float FofVar = function->analysisFunction(at, function->coefficientCount, function->coefficients, function->coefficientNames);

	result = FofVar * currentInstance.numericalDelta;

	pthread_mutex_unlock(&threadLock);

	return result;
}

float garlicAnalysisDeriveFunction(struct garlicFunction* function, float at) {
	pthread_mutex_lock(&threadLock);
	privateGarlicContextSwitch();

	float result = 0.0f;

	float dvar = currentInstance.numericalDelta;

	float FofVar = function->analysisFunction(at, function->coefficientCount, function->coefficients, function->coefficientNames);

	float FofVarPlusDVar = function->analysisFunction(at + dvar, function->coefficientCount, function->coefficients, function->coefficientNames);

	result = (FofVarPlusDVar - FofVar) / dvar;

	pthread_mutex_unlock(&threadLock);

	return result;
}


