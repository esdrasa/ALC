#include <stdio.h>
#include "operacoes.h"
#include "metodos.h"
#include "saida.h"
#include "entrada.h"

int main() {
	FILE *arquivo;
	int escolha;
	char leitura;

	do {
	escolha = 0;
	system("clear");
	printf("Ola, escolha um dos seguintes topicos para ampliar seu conhecimento\n");
	printf("1) Fatoracao Cholesky\n");
	printf("2) Metodo das Potencias\n");
	printf("3) Metodo iterativo SOR\n");
	printf("4) Fatoracao LU\n");
	printf("5) PSO/Simulated annelling\n");
	printf("6) Sair\n\n");

	printf("Digite o numero do topico escolhido\n");
	//escolhendo o topico    
	while(escolha < 1 || escolha > 6) {
		scanf("%d",&escolha);
		printf("\n");
		if(escolha < 1 || escolha > 6 )
			printf("Escolha um topico existente\n");
	}

	if(escolha == 1)
	{
		system("clear");
		if( (arquivo = fopen("cholesky.txt","r") ) == NULL )
			printf("Nao foi possivel abrir o arquivo\n");
		else 
		{
				while(!feof(arquivo))
				{	
					fscanf(arquivo, "%c", &leitura);
					printf("%c",leitura );					
				}

				int teste,n;
				double **A, *b, **R, **Rt, *y, *x1;
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

			        liberaVetor(b);
			        liberaVetor(y);
			        liberaVetor(x1);
			        liberaMatriz(A,n);
			        liberaMatriz(Rt,n);
			        liberaMatriz(R,n);
			    }
			    
		}//fecha o else
		printf("Pressione enter para voltar ao Menu.\n");
		getchar();
	}//fecha o if da escolha 1

	 
	 if(escolha == 2)
	 {
	 	system("clear");
		if( (arquivo = fopen("potencia.txt","r") ) == NULL )
			printf("Nao foi possivel abrir o arquivo\n");

		else 
		{
				while(!feof(arquivo))
				{	
					fscanf(arquivo, "%c", &leitura);
					printf("%c",leitura );					
				}

				int n, i;
				double **A,
				        *x,
				        tol = 0.000001,
				        autovalor;

				//Lendo os dados necessários

				printf("Digite a dimensão da matriz\n");
				scanf("%d%*c", &n);
				printf("\n");

				A = lerMatriz(n, n);
				x = lerVetor(n);

				//Chamando o metodo da potencia.
				autovalor = potencia(n, tol, A, x);

				printf("\nO Autovalor Dominante e: %lf\n",autovalor );

				   //liberando espaço alocado
				liberaVetor(x);
				liberaMatriz(A, n);
		}//chave do else
		printf("Pressione enter para voltar ao Menu.\n");
		getchar();
	}//chave do if da escolha 2

	if(escolha == 3)
	{
		system("clear");
		if( (arquivo = fopen("sor.txt","r") ) == NULL )
			printf("Nao foi possivel abrir o arquivo\n");

		else 
		{
				while(!feof(arquivo))
				{	
					fscanf(arquivo, "%c", &leitura);
					printf("%c",leitura );					
				}

			double **A, *b, *x2, tol, w;
		    int n, teste;
		    unsigned long int iMax;

		    printf("Digite o tamanho da matriz: \n");
		    scanf("%d%*c", &n);

			A = lerMatriz(n, n);

    		b = lerVetor(n);

			printf("Digite a tolerancia para o erro relativo: \n");
		    scanf("%lf%*c", &tol);

		    printf("Digite um valor para w entre 0 e 2: \n");
		    scanf("%lf%*c", &w);

		    printf("Digite o numero maximo de iteracoes permitidas: \n");
		    scanf("%lu%*c", &iMax);

		    x2 = criaVetor(n);

		    if(teste = SOR(A, b, x2, tol, w, iMax, n))
		    {
		        printf("\nSolucao do sistema com o metodo SOR:\n");
		        imprimeVetor(x2, n);
		    }
		    else
		    {
		        printf("\nNao foi possivel resolver o sistema pelo metodo SOR.\n");
		    }
		    liberaMatriz(A,n);
		    liberaVetor(b);
		    liberaVetor(x2);
		}
		printf("Pressione enter para voltar ao Menu.\n");
		getchar();
	}//chave do if da escolha 3


	if(escolha == 4)
	{
		system("clear");
		if( (arquivo = fopen("lu.txt","r") ) == NULL )
		printf("Nao foi possivel abrir o arquivo\n");

		else
		{
				while(!feof(arquivo))
				{	
					fscanf(arquivo, "%c", &leitura);
					printf("%c",leitura );					
				}

				double **A, **l, **u, *b, *b2,*x;
				int n;


			    printf("Digite o tamanho da matriz: \n");
			    scanf("%d%*c", &n);

				A = lerMatriz(n, n);
    			b = lerVetor(n);
    			u = criaMatriz(n, n);
    			l = criaMatrizI(n);
    			b2 = criaVetor(n);
			    x = criaVetor(n);

				lu(A, l, u, n);
    
			    forwardSub(l, b, b2, n);
			    
			    backSub(u, b2, x, n);
			    
			    printf("Matriz L:\n");
			    imprimeMatriz(l, n, n);

			    printf("Matriz U:\n");
			    imprimeMatriz(u, n, n);			    
			    
			    printf("Solucao obtida com o metodo LU:\n");
			    imprimeVetor(x, n);

			    liberaMatriz(A,n);
			    liberaMatriz(l,n);
			    liberaMatriz(u,n);
			    liberaVetor(b);
			    liberaVetor(b2);
			    liberaVetor(x);
		}
		printf("Pressione enter para voltar ao Menu.\n");
		getchar();
	}//fecha o if da escolha 4

	if(escolha == 5) {
		if( (arquivo = fopen(".txt","r") ) == NULL )
			printf("Nao foi possivel abrir o arquivo\n");
	}

	if(escolha == 6) {
		return 0;
	}
		} while(1);

	fclose(arquivo);
	return 0;
}