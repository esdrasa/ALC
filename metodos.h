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

#endif // METODOS_H_INCLUDED
