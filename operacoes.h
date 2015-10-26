#ifndef OPERACOES_H_INCLUDED
#define OPERACOES_H_INCLUDED

/**
Retorna uma matriz quadrada de ordem m x n
*/
double** criaMatriz(int m, int n);

/**
Retorna uma matriz identidade de tamanho n
*/
double** criaMatrizI(int n);

/**
Retorna a matriz transposta de uma matriz de ordem m x n
*/
double** transposta(double** matriz, int m, int n);

/**
Retorna um vetor de tamanho n
*/
double* criaVetor(int n);

/**
Realiza operação de substituição para trás com a intenção de resolver um sistema linear.
Sistema Ax = b, no qual A = matriz e x é o retorno da função.
A é uma matriz triangular superior de tamanho n, e b um vetor de tamanho n.
*/
double* backSub(double** matriz, double* b, int n);

/**
Realiza operação de substituição para frente com a intenção de resolver um sistema linear.
Sistema Ax = b, no qual A = matriz e x é o retorno da função.
A é uma matriz triangular inferior de tamanho n, e b um vetor de tamanho n.
*/
double* forwardSub(double** matriz, double* b, int n);

/**
Retorna a multiplicação entre duas matrizes A e B, onde ab = retorno.
A matriz A possuindo ordem m x n, e a matriz B p x q
*/
double** multiplica(double** a, int m, int n, double** b, int p, int q);

/**
Retorna a multiplicação de um vetor v de dimensão s por uma matriz a de ordem m x n.
*/
double* multiplicaVetor(double** a, int m, int n, double* v, int s);

/**
Realiza e retorna o produto escalar entre a e b. <a, b> =  retorno.
Ambos de tamanho n.
*/
double produtoEscalar(double* a, double* b, int n);

/**
Retorna a subtração do vetor b ao vetor a. Ambos de tamanho n
*/
double* subtraiVetores(double* a, double* b, int n);

/**
Retorna a soma entre dois vetores de tamanho n
*/
double* somaVetores(double* a, double* b, int n);

/**
Retorna a soma entre duas matrizes de tamanho n.
*/
double** somaMatrizes(double** a, double** b, int n);

/**
Retorna da matriz b à matriz a. Ambas de tamanho n.
*/
double** subtraiMatrizes(double** a, double** b, int n);

/**
Calcula a norma de Frobenius de uma matriz de tamanho n
*/
double normaF(double** matriz, int n);

/**
Calcula a norma linha (infinita) de uma matriz de tamanho n
*/
double normaLinha(double** matriz, int n);

/**
Calcula a norma coluna (norma 1) de uma matriz de tamanho n
*/
double normaColuna(double** matriz, int n);

/**
Calcula a norma dois de um vetor de tamanho n
*/
double normaDois(double* vetor, int n);

/**
Calcula o angulo entre dois vetores de tamanho n
*/
double angulo(double* v1, double* v2, int n);

/**
Verifica se a matriz A é uma matriz de Vandermonde
*/
int isVandermonde(double** a, int n);

/**
Calcula o determinante de uma matriz de Vandermonde
*/
double detVandermonde(double** a, int n);

/**
Retorna um vetor com os elementos de uma coluna selecionada em uma matriz
*/
double* getColuna(double** matriz, int coluna, int n);

/**
Calcula a norma infinito de um vetor
*/
double normaInfinito(double* vetor, int n);

/**
Calcula o erro relativo entre dois números reais x1 e x0:
erro = |x1 - x0| / |x1|
*/
double erroRelativo(double x1, double x0);

#endif // OPERACOES_H_INCLUDED