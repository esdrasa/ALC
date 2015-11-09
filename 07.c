#include <stdio.h>
#include "operacoes.h"
#include "metodos.h"
#include "saida.h"
#include "entrada.h"

int main()
{
    double **A, *b, *x, cLinha, cColuna, cNorma, cSassenfeld, tol;
    int n, escolha = 0, ok = 0;
    unsigned long int iMax;

    printf("Digite o tamanho do sistema: ");
    scanf("%d%*c", &n);

    x = criaVetor(n);

    A = lerMatriz(n, n);

    b = lerVetor(n);

    cLinha = criterioLinhas(A, n);
    cColuna = criterioColunas(A, n);
    cNorma = criterioNorma(A, n);
    cSassenfeld = criterioSassenfeld(A, n);

    printf("Criterio das Linhas: ");
    if(cLinha)
        printf("[X]");
    else
        printf("[ ]");
    printf(" %lf\n", cLinha);


    printf("Criterio das Colunas: ");
    if(cColuna)
        printf("[X]");
    else
        printf("[ ]");
    printf(" %lf\n", cColuna);


    printf("Criterio das Normas: ");
    if(cNorma)
        printf("[X]");
    else
        printf("[ ]");
    printf(" %lf\n", cNorma);


    printf("Criterio de Sassenfeld: ");
    if(cSassenfeld)
        printf("[X]");
    else
        printf("[ ]");
    printf(" %lf\n\n", cSassenfeld);
    
    //O algoritmo ficará parado no loop abaixo caso o usuário não escolha entre um dos dois métodos.
    while(escolha != 1 && escolha != 2)
    {
        printf("Escolha um metodo para resolver o sistema:\n");
        printf("[1] - Metodo de Jacobi\n");
        printf("[2] - Metodo de Gauss-Seidel\n");

        scanf("%d%*c", &escolha);
    }

    printf("\nDigite a tolerancia para o erro relativo: ");
    scanf("%lf%*c", &tol);

    printf("Digite o numero maximo de iteracoes permitidas: ");
    scanf("%lu%*c", &iMax);

    if(escolha == 1 && (cNorma || cColuna || cLinha))
    {
	//Se há garantia de convergência, então a tolerância pode ser zero.
        ok = jacobi(A, b, x, 0, iMax, n);
    }
    else if(escolha == 1 && !(cNorma || cColuna || cLinha))
    {
	//Nesse caso a tolerância tem que ser a escolhida pelo usuário, pois não garantia de convergência.
        ok = jacobi(A, b, x, tol, iMax, n);
    }
    else if(cNorma || cColuna || cLinha || cSassenfeld)
    {
	//Se há garantia de convergência, então a tolerância pode ser zero.
        ok = gaussSeidel(A, b, x, 0, iMax, n);
    }
    else if(!(cNorma || cColuna || cLinha || cSassenfeld))
    {
	//Nesse caso a tolerância tem que ser a escolhida pelo usuário, pois não garantia de convergência.
        ok = gaussSeidel(A, b, x, tol, iMax, n);
    }

    if(!ok)
    {
        printf("Nao foi possivel resolver o sistema.\n");
    }
    else
    {
        printf("Solucao do sistema: \n");
        imprimeVetor(x, n);
    }

    getchar();

    return 0;
}