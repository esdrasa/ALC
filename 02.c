#include <stdio.h>
#include "operacoes.h"
#include "metodos.h"
#include "saida.h"
#include "entrada.h"

int main()
{
    double **A, *b, **R, **Rt, *y, *x1, *x2, *r1, *r2, tol, w;
    int n, teste, teste2;
    unsigned long int iMax;

    printf("Digite o tamanho da matriz: ");
    scanf("%d%*c", &n);

    A = lerMatriz(n, n);

    b = lerVetor(n);

    R = criaMatriz(n, n);
    
    //Se a matriz A for definida positiva e simétrica, a matriz R (fator de Cholesky) será calculada.
    teste = cholesky(A, R, n);

    if(!teste)
    {
        printf("Nao foi possivel determinar o fator de Cholesky\n");
    }
    else
    {
        printf("Fator de Cholesky encontrado:\n");
        imprimeMatriz(R, n, n);
	
	/*
	 * Rt = matriz transposta de R.
	 * temos A e b em um sistema A*x = b, e queremos determinar x
	 * A = Rt*R
	 * Rt*R*x = b
	 * R*x = y
	 * Rt*y = b (y será resolvido com substituição para frente, pois Rt é triangular inferior).
	 * Ao descobrir os valores de y, podemos achar os valores de x, pois
	 * R*x = y (x será resolvido com substituição para trás, pois R é triangular superior).
	*/
	
	Rt = criaMatriz(n, n);
        transposta(R, Rt, n, n);
	
	y = criaVetor(n);
        forwardSub(Rt, b, y, n);
	
	x1 = criaVetor(n);
        backSub(R, y, x1, n);

        printf("Solucao do sistema obtida com o fator de Cholesky:\n");

        imprimeVetor(x1, n);
    }

    printf("Resolucao do sistema com o metodo SOR:\n");

    printf("\nDigite a tolerancia para o erro relativo: ");
    scanf("%lf%*c", &tol);

    printf("Digite um valor para w entre 0 e 2: ");
    scanf("%lf%*c", &w);

    printf("Digite o numero maximo de iteracoes permitidas: ");
    scanf("%lu%*c", &iMax);

    x2 = criaVetor(n);

    if(teste2 = SOR(A, b, x2, tol, w, iMax, n))
    {
        printf("\nSolucao do sistema com o metodo SOR:\n");
        imprimeVetor(x2, n);
    }
    else
    {
        printf("\nNao foi possivel resolver o sistema pelo metodo SOR.\n");
    }
    
    //Compara a norma do resíduo das soluções obtidas com os métodos SOR e Cholesky, caso ambos sejam tenham convergido
    if(teste && teste2)
    {
	r1 = criaVetor(n);
	r2 = criaVetor(n);
	
        residuo(A, x1, b, r1, n);
        residuo(A, x2, b, r2, n);

        printf("\nDiferenca entre os resultados obtidos com os metodos SOR e Cholesky:\n");
        printf("Normas do residuo das solucoes obtidas em cada metodo:\n");
        printf("Cholesky: %lf\n", normaDois(r1, n));
        printf("SOR: %lf\n", normaDois(r2, n));
    }

    getchar();

    return 0;
}