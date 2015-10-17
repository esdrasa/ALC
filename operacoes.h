#ifndef OPERACOES_H_INCLUDED
#define OPERACOES_H_INCLUDED

/**
Retorna uma matriz quadrada de tamanho n
*/
double** criaMatriz(int n);

/**
Retorna uma matriz identidade de tamanho n
*/
double** criaMatrizI(int n);

/**
Realiza a operação de transposição em uma matriz m de tamanho n
*/
void transpor(double** m, int n);

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
Ambos de tamanho n.
*/
double** multiplica(double** a, double** b, int n);

/**
Retorna a multiplicação de um vetor v por uma matriz a.
Ambos de tamanho n.
*/
double* multiplicaVetor(double** a, double* v, int n);

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


#endif // OPERACOES_H_INCLUDED
