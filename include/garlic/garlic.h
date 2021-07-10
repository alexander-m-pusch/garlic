/*
 * cnarli.h
 *
 *  Created on: 09.07.2021
 *      Author: alexander
 */

#ifndef INCLUDE_GARLIC_GARLIC_H_
#define INCLUDE_GARLIC_GARLIC_H_

/**
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

struct garlicMatrix {
	int rows;
	int columns;
	float** elements;
};

struct garlicFunction {
	float (*analysisFunction) (float variable, int coefficientCount, float* coefficients, char** coefficientNames);
	int coefficientCount;
	float* coefficients;
	char** coefficientNames;
};

struct garlicInstance {
	void* (*allocatorFunction) (long unsigned int size);
	void (*deallocatorFunction) (void* address);
	int numericalIterations;
	float numericalDelta;
};

struct garlicInstance garlicInit(void* (*allocatorFunction) (long unsigned int size), void (*deallocatorFunction) (void* address), int numericalIterations, float numericalDelta);

void garlicTerminate(struct garlicInstance instance);

void garlicSetInstance(struct garlicInstance instance);

struct garlicInstance garlicGetInstance(void);

struct garlicMatrix* garlicLinearAlgebraCreateMatrix(int rows, int columns);

struct garlicMatrix* garlicLinearAlgebraTransposeMatrix(struct garlicMatrix* toTranspose);

struct garlicMatrix* garlicLinearAlgebraCalculateInverseOfMatrix(struct garlicMatrix* toInvert);

struct garlicMatrix* garlicLinearAlgebraMultiplyMatrices(struct garlicMatrix** matrixA, int count);

struct garlicMatrix* garlicLinearAlgebraAddMatrices(struct garlicMatrix** matrices, int count);

float garlicLinearAlgebraGetMatrixDeterminant(struct garlicMatrix* matrix);

struct garlicMatrix* garlicLinearAlgebraSolveMatrix(struct garlicMatrix* matrix);

struct garlicMatrix* garlicLinearAlgebraCrossProduct(struct garlicMatrix** vectors, int count);

float garlicLinearAlgebraDotProduct(struct garlicMatrix* vector1, struct garlicMatrix* vector2);

struct garlicMatrix* garlicLinearAlgebraMultiplyMatrixWithScalar(struct garlicMatrix* matrix, float scalar);

float garlicArithmeticSum(float* numbers, int n);

float garlicArithmeticProduct(float* numbers, int n);

struct garlicFunction* garlicAnalysisCreateFunction(float (*functionToWrap) (float variable, int coefficientCount, float* coefficients, char** coefficientNames), int coefficientCount);

void garlicAnalysisSetCoefficient(struct garlicFunction* function, char* coefficient, float value, int coefficientIndex);

float garlicAnalysisUseFunction(struct garlicFunction* function, float inputValue);

float garlicAnalysisDefiniteIntegral(struct garlicFunction* function, float lowerBound, float upperBound);

float garlicAnalysisIntegrateFunction(struct garlicFunction* function, float at);

float garlicAnalysisDeriveFunction(struct garlicFunction* function, float at);

#endif /* INCLUDE_GARLIC_GARLIC_H_ */
