#ifndef METODOS_H_INCLUDED
#define METODOS_H_INCLUDED

/**
Calcula e retorna o resíduo da solução de um sistema linear.
Um sistema linear Ax = b, onde A é a matriz e x a solução encontrada.
O cálculo do resíduo será: r = b - Ax
*/
double* residuo(double** matriz, double* x, double* b, int n);

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

#endif // METODOS_H_INCLUDED