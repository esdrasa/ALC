#ifndef METODOS_H_INCLUDED
#define METODOS_H_INCLUDED

/**
Calcula o resíduo (r) da solução de um sistema linear.
Um sistema linear Ax = b, onde A é a matriz e x a solução encontrada.
O cálculo do resíduo será: r = b - Ax
*/
void residuo(double** matriz, double* x, double* b, double* r, int n);

/**
Realiza a fatoração LU na matriz A de tamanho n.
Os argumentos l e u são as matrizes L e U, respectivamente.
l e u são iniciados antes da chamada à função.
A função retornará 0 caso algum valor na diagonal principal de u seja 0,
e retornará 1 não ocorra problemas.
*/
int lu(double** a, double** l, double** u, int n);

/**
Aplica o método iterativo de Jacobi para a resolução do sistema Ax = b.
tolerancia = o erro relativo aceitável entre a resulução de uma iteração anterior e a atual.
max = número máximo de iterações.
x = vetor solução alocado antes da chamada à função.
A função retorna 0 caso não seja possível resolver o sistema, e 1 caso possível.
*/
int jacobi(double** A, double* b, double* x, double tolerancia, unsigned long int iMax, int n);

/**
Aplica o método iterativo de Gauss-Seidel para a resolução do sistema Ax = b.
tolerancia = o erro relativo aceitável entre a resulução de uma iteração anterior e a atual.
max = número máximo de iterações.
x = vetor solução alocado antes da chamada à função.
A função retorna 0 caso não seja possível resolver o sistema, e 1 caso possível.
*/
int gaussSeidel(double** A, double* b, double* x, double tolerancia, unsigned long int iMax, int n);

/**
Aplica o SOR (Sucessive Over Relaxation) para a resolução do sistema Ax = b.
tolerancia = o erro relativo aceitável entre a resulução de uma iteração anterior e a atual.
max = número máximo de iterações.
x = vetor solução alocado antes da chamada à função.
w = termo que acelera a convergência para a solução.
A função retorna 0 caso não seja possível resolver o sistema, e 1 caso possível.
*/
int SOR(double** A, double* b, double* x, double tolerancia, double w, unsigned long int iMax, int n);

/**
Retorna o valor encontrado se o critério das Linhas da matriz A é válido
*/
double criterioLinhas(double** A, int n);

/**
Retorna o valor encontrado se o critério das Colunas da matriz A é válido
*/
double criterioColunas(double** A, int n);

/**
Retorna o valor encontrado se o critério da Norma é válido para a matriz A
*/
double criterioNorma(double** A, int n);

/**
Retorna o valor encontrado se o critério de Sassenfeld é válido para a matriz A
*/
double criterioSassenfeld(double** A, int n);

/**
Calcula o fator de Cholesky R da matriz A. Este fator será utilizado para resolver
um sistema Ax = b, onde A=R^tR; R^t * y = b; Rx = y
*/
int cholesky(double **A, double **R, int n);

/**
 * Método do Gradiente Conjugado para calcular a solução do sistema Ax = b.
 * tol = erro relativo aceitável entre a resolução de uma iteração anterior e a atual.
 * iMax = número máximo de iterações.
 */
int gradienteConjugado(double** A, double* b, double* x, double tol, unsigned long int iMax, int n);

/**
 * Calcula uma raíz de um polinômio, se houver.
 */
int newton(double *polinomio, double* raiz, double tolerancia, unsigned long int iMax, int n);

/**
 * Método de Gram-Schimdt para determinar uma base ortonormal de um conjunto de vetores.
 */
void schmidt(double** vetores, double** base, int dimensao);

/**
 * Método de Gram-Schimdt modificado.
 */
void schmidtModificado(double** vetores, double** base, int dimensao);

/**
	Metodo da potencia para calcular o Autovalor dominante de uma matriz.
	n = Tamanho da matrix nxn
	tol = tolerancia do metodo
*/
double potencia(int n, double tol, double **A, double *x);

/**
 * Simulated Annealing.
 * alfa = constante para variação do limite superior
 * ls = limite superior
 * li = limite inferior
 * nr = número de repetições a cada iteração
 */
void simulatedAnnealing(double** A, double* b, double* x, double alfa, double ls, double li, int nr, int n);
 
 /**
 MB é o vetor solução
 */
void PSO(double **matrizA, int n, double *vetorB, double *MB, double tol, int particulas);

#endif // METODOS_H_INCLUDED