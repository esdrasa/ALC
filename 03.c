#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "metodos.h"
#include "operacoes.h"
#include "saida.h"
#include "entrada.h"

/** O resultado não é encontrado logo na primeira execução. Sei que devo armazenar a melhor resposta, porém estou tendo um problema estranho. 
A ideia é usar a norma do residuo da primeira solução como base. Caso a norma do resíduo de algum outro vetor gerado seja menor, os valores desse
vetor irão ser colocados no vetorXmelhor e a normaMelhor ira ser atualizada com a norma do residuo desse melhor vetor para sempre ter uma comparação





int main(void){
    //Simulated Annealing

    double **matrizA, *vetorB, *vetorXatual, *vetorXtemp, **matrizAt, *vetorXmelhor;
    double *residuo1, *residuo2, *residuoMelhor, variacaoEnergia, nRandom, normaMelhor;

    int n, m, i, j, k, l, iteracoes;

    printf("Digite a quantidade de linhas:\n");
    scanf("%d",&n);

    printf("Digite a quantidade de colunas:\n");
    scanf("%d",&m);

    printf("\nDigite a quantidade de iteracoes desejada\n");
    scanf("%d",&iteracoes);

    matrizA     = lerMatriz(n, m);
    vetorB      = lerVetor(n);

    vetorXatual = criaVetor(n);
    vetorXtemp  = criaVetor(n);
    vetorXmelhor = criaVetor(n);
    residuo1    = criaVetor(n);
    residuo2    = criaVetor(n);
    matrizAt   = criaMatriz(m, n);

    /// transpor a matriz A

    transposta(matrizA, matrizAt, n, m); /// A partir de agora será usada a matriz AtA como sendo a A




    /**PEGAR SOLUCAO INICIAL(ALEATÓRIA)*/
    srand((unsigned)time(NULL));

    for(i=0; i<n; i++) {

        vetorXatual[i] = rand() % 100;
        vetorXmelhor[i] = vetorXatual[i]; ///Aqui, assumo que o primeiro vetor gerado é o melhor
    }
    /// calcular a norma do residuo do possivel melhor vetorX.  Não consegui ainda pegar o melhor vetor utilizando a menor norma do residuo
    residuo(matrizAt, vetorXatual, vetorB, residuoMelhor, n);
    normaMelhor = normaDois(residuoMelhor, n);

        /// Definir a temperatura máxima e mínima

        float temperaturaMAX = 100;
        float temperaturaMIN = 0.01;
        double normaDoResiduo1;

        while(temperaturaMAX > temperaturaMIN) {

            /// gerando um vetortemp solução "aleatório"
            for(j=0; j<iteracoes; j++)
            {
                for(k=0; k<n; k++)
                {
                    vetorXtemp[k] = rand() % 100;
                }
                ///calculando o residuo com as duas solucoes
                residuo(matrizAt, vetorXatual, vetorB, residuo1, n);
                residuo(matrizAt, vetorXtemp, vetorB, residuo2, n);
                
                ///checar substituicao da melhor solucao. Era para funcionar de boa, mas acontece um erro. Não consegui atualizar a normaMelhor e o 
                vetor ao mesmo tempo
                if(normaMelhor > normaDoResiduo1)
                {
                    normaMelhor = normaDoResiduo1;
                    for(l=0; l<n; l++)
                    {
                        vetorXmelhor[l] = vetorXatual[l];
                    }
                }

                /// calculando a variação da 2 solucoes
                normaDoResiduo1 = normaDois(residuo1, n);
                variacaoEnergia = normaDois(residuo2, n) - normaDoResiduo1;

                ///Agora as duas condições para trocar o vetor solução atual

                if(variacaoEnergia <= 0)
                {
                    for(i=0; i<n; i++)
                    {
                        vetorXatual[i] = vetorXtemp[i];
                    }
                }
                else{
                        nRandom = rand() % 100;
                        nRandom /= 100; /// gera um numero aleatorio entre 0 e 1

                        if(nRandom < exp(-variacaoEnergia / temperaturaMAX)) /// checa se esse numero é menor
                        {
                            for(i=0; i<n; i++)
                            {
                                vetorXatual[i] = vetorXtemp[i];
                            }
                        }
                }
            }
            temperaturaMAX = temperaturaMAX * 0.9;
        }

        printf("\n\nSolucao do sistema\n\n");
        imprimeVetor(vetorXatual,n);

       /* printf("\n\n melhor solucao encontrada\n\n");
        imprimeVetor(vetorXmelhor,n);*/

        system("pause");



    return 0;
}