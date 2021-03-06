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
Calcula a matriz transposta () de uma matriz de ordem m x n
*/
void transposta(double** matriz, double** t, int m, int n);

/**
Retorna um vetor de tamanho n
*/
double* criaVetor(int n);

/**
Realiza operação de substituição para trás com a intenção de resolver um sistema linear.
Sistema Ax = b, no qual A = matriz e x é a solução da função.
A é uma matriz triangular superior de tamanho n, e b um vetor de tamanho n.
*/
void backSub(double** matriz, double* b, double* x, int n);

/**
Realiza operação de substituição para frente com a intenção de resolver um sistema linear.
Sistema Ax = b, no qual A = matriz e x é a solução da função.
A é uma matriz triangular inferior de tamanho n, e b um vetor de tamanho n.
*/
void forwardSub(double** matriz, double* b, double* x, int n);

/**
Multiplica duas matrizes A e B, onde result = ab.
A matriz A possuindo ordem m x n, e a matriz B p x q
*/
void multiplica(double** a, int m, int n, double** b, int p, int q, double** result);

/**
Multiplica um vetor v de dimensão s por uma matriz a de ordem m x n.
result = Av
*/
void multiplicaVetor(double** a, int m, int n, double* v, double* result, int s);

/**
Realiza e retorna o produto escalar entre a e b. <a, b> =  retorno.
Ambos de tamanho n.
*/
double produtoEscalar(double* a, double* b, int n);

/**
Subtrai o vetor b do vetor a. Ambos de tamanho n.
result = a - b
*/
void subtraiVetores(double* a, double* b, double* result, int n);

/**
Soma dois vetores de tamanho n.
result = a + b
*/
void somaVetores(double* a, double* b, double* result, int n);

/**
Soma duas matrizes de tamanho n.
result = A + B
*/
void somaMatrizes(double** a, double** b, double** result, int n);

/**
Subtrai b da matriz a. Ambas de tamanho n.
result = A - B
*/
void subtraiMatrizes(double** a, double** b, double** result, int n);

/**
Calcula a norma de Frobenius de uma matriz de tamanho n
*/
double normaF(double** matriz, int n);

/**
Calcula a norma linha (infinita) de uma matriz de tamanho n.
A norma linha é o maior número no conjunto das somas dos valores absolutos
dos elementos de cada linha.
*/
double normaLinha(double** matriz, int n);

/**
Calcula a norma coluna (norma 1) de uma matriz de tamanho n.
A norma coluna é o maior número no conjunto das somas dos valores absolutos
dos elementos de cada coluna.
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
Transforma vColuna em um vetor com os elementos de uma coluna selecionada em uma matriz
*/
void getColuna(double** matriz, double* vColuna, int coluna, int n);

/**
Calcula a norma infinito de um vetor (maior valor desse vetor)
*/
double normaInfinito(double* vetor, int n);

/**
Calcula o erro relativo entre dois números reais x1 e x0:
erro = |x1 - x0| / |x1|
*/
double erroRelativo(double x1, double x0);

/**
 * Libera o espaço alocado para um ponteiro de ponteiro na memória
 */
void liberaMatriz(double** a, int nLinhas);

/**
 * Libera o espaço alocado para um ponteiro na memória
 */
void liberaVetor(double* a);

/**
 * Calcula a derivada d do polinômio y.
 */
void derivada(double* y, double* d, int n);

/**
 * Retorna o valor do polinomio P(x) quando x = abscissa.
 */
double ordenada(double* polinomio, double abscissa, int n);

/** 
Preenche com zeros uma matriz qualquer.
*/
void **matrizNula(double **matriz, int n, int m);

/**
Encontra o grau de cada vértice em uma matriz de adjacência.
*/
void encontraGraus(double **matrizA, int n, int m, double **matrizD);

/**
Encontra o determinante de uma matriz identidade.
*/
double determinanteI(double **matrizI, int n);

/**
 * Multiplica o "vetor" por um "escalar". O vetor "resultado" receberá o novo vetor.
 */
void multiplicaPorEscalar(double* vetor, double* resultado, double escalar, int n);

/**
 * Realiza a projeção ortogonal do vetor w sobre o vetor v.
 */
void projecaoOrtogonal(double* w, double* v, double* proj, int n);

/**
 * Atribui os elementos do vetor x ao vetor y.
 */
void atribui(double* y, double* x, int n);

/**
 * Verifica se um conjunto de vetores são ortonormais entre si, dois a dois
 */
int isOrtonormal(double** vetores, int dimensao);

/**
 * Verifica se uma matriz está em banda.
 */
int isBanda(double** matriz, int n);

/**
  *Divide um vetor por um escalar
*/
void divideVetorPorEscalar(double *v, int n, double escalar);
  
/**
Preenche a coluna de uma matriz com um vetor desejado
*/
void preencheColuna(double **matriz, int coluna, int n, double *vetor);

/**
Atualiza a velocidade no PSO
*/
void atualizaVel(double A1, double A2, double A3, double c1, double c2, int n, double *melhorPos, double *vetVel, double *vetPos, double *vetMB, double *vetResul);

/**
Atualiza a posição no PSO
*/
void atualizaPos(double *Satual, double *velAtual, int n);

/**
Preenche com zeros um vetor qualquer.
*/
void vetorNulo(double *vetor, int n);

#endif // OPERACOES_H_INCLUDED