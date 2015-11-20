#include <stdio.h>
#include "operacoes.h"
#include "metodos.h"
#include "saida.h"
#include "entrada.h"

int main()
{
    double **vetores, **base1, **base2;
    int n;
    
    printf("Digite a dimensao dos vetores: ");
    scanf("%d%*c", &n);
    
    //Lê um conjunto de vetores do arquivo matriz.txt
    vetores = lerMatriz(n, n);
    
    base1 = criaMatriz(n, n);
    
    base2 = criaMatriz(n, n);
    
    if(isOrtonormal(vetores, n))
    {
	printf("Os vetores são ortonormais.\n");
    }
    else
    {
	schmidt(vetores, base1, n);
	
	schmidtModificado(vetores, base2, n);
	
	printf("Conjunto de vetores ortonormalizados pelo metodo de Gram-Schimdt:\n");
	imprimeMatriz(base1, n, n);
	
	printf("Conjunto de vetores ortonormalizados pelo metodo de Gram-Schimdt modificado:\n");
	imprimeMatriz(base2, n, n);
    }
    
    getchar();
    
    return 0;
}