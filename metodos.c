#include "operacoes.h"
#include "saida.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <float.h>

/*
 * Em métodos iterativos, quando a tolerancia passada como parâmetro for igual 0, será utilizada
 * a constante DBL_EPSILON, da biblioteca float.h, no lugar de zero, pois DBL_EPSILON é a menor
 * diferença possível, na linguagem em uso, entre 1 e o menor valor maior que 1.
*/


void residuo(double** matriz, double* x, double* b, double* r, int n)
{
    double* Ax = criaVetor(n);
    
    multiplicaVetor(matriz, n, n, x, Ax, n);
    
    subtraiVetores(b, Ax , r, n);
    
    liberaVetor(Ax);
}

int lu(double** a, double** l, double** u, int n)
{
    int k, i, j;

    /* Os dois loops abaixo não são necessários caso a matriz "A" seja utilizada para calcular a matriz "U".
    A não utilização deles pode vir a ser útil no tempo de processamento, já que estes loops apenas atribuem os
    valores de A à U.*/
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            u[i][j] = a[i][j];

    for(k = 1; k < n; k++)
    {
        if(u[k-1][k-1] == 0)
            return 0;
	
	/*k será igual a quantidade de equações já escalonadas. Assim, a linha linha k deverá ser escalonada
	nos próximos loops, pois a linhas i = 0...k-1 já estarão escalonadas e o objetivo é escalonar todas as linhas.
	*/
        for(i = k; i < n; i++)
        {
	    /*l[i][k-1] será o coeficiente utilizado para zerar cada elemento abaixo do elemento u[k-1][k-1] da
	    diagonal principal. Sendo u[k-1][k-1] o primeiro coeficente diferente de zero da equação já escalonada;
	    e u[i][k-1] o coeficiente a ser zerado na equação (linha) i.
	    */
            l[i][k-1] = u[i][k-1] / u[k-1][k-1];

	    /*
	    Para zerar os coeficientes abaixo do elemento u[k-1][k-1], multiplica-se o coeficiente l
	    descoberto acima por cada coeficiente da última equação já escalonada, e subtrai-se cada
	    multiplicação obtida anteriormente ao correspondente coeficiente das equações abaixo da linha k-1.
	    */
            for(j = k-1; j < n; j++)
            {
		/*
		 Lembrando que quando j = k-1, u[i][k-1] = u[i][k-1] - l[i][k-1]*u[k-1][k-1]. Mas
		 l[i][k-1] = u[i][k-1] / u[k-1][k-1]. Então: u[i][k-1] = u[i][k-1] - (u[i][k-1] / u[k-1][k-1]) * u[k-1][k-1].
		 Mas, (u[i][k-1] / u[k-1][k-1]) * u[k-1][k-1] = u[i][k-1]. Logo, u[i][k-1] será igual a zero, que é o objetivo.
		*/
                u[i][j] = u[i][j] - l[i][k-1]*u[k-1][j];
            }
        }
    }

    if(u[n-1][n-1] == 0)
        return 0;

    return 1;
}

int jacobi(double** A, double* b, double* x, double tolerancia, unsigned long int iMax, int n)
{
    int i, j;
    double erro = tolerancia + 1, *xant, soma;
    unsigned long int k = 0;
    
    if(tolerancia == 0)
	tolerancia = DBL_EPSILON;

    xant = criaVetor(n);

    for(i = 0; i < n; i++)
    {
        /*
        Se algum elemento da diagonal principal for zero, não haverá solução,
        a menos que as linha da matriz sejam permutadas.
        */
        if(A[i][i] == 0)
	{
	    liberaVetor(xant);
            return 0;
	}

        //primeiro vetor para solução do método iterativo
        x[i] = b[i] / A[i][i];
    }

    while(k < iMax && erro >= tolerancia)
    {
        erro = 0;

        for(i = 0; i < n; i++)
            xant[i] = x[i];

        for(i = 0; i < n; i++)
        {
            soma = 0;

            for(j = 0; j < n; j++)
                if(j != i)
                    soma += A[i][j] * xant[j];

            x[i] = (b[i] - soma) / A[i][i];

            if(erroRelativo(x[i], xant[i]) > erro)
                erro = erroRelativo(x[i], xant[i]);
        }

        k++;
    }

    /*
    Se erro for maior que a tolerancia no final das iterações, significa que
    não ocorreu convergência para o número máximo de iterações definidas.
    */
    if(erro >= tolerancia)
        return 0;
    
    for(i = 0; i < n; i++)
	    if(isnan(x[i]) || isinf(x[i]))
	        return 0;

    return 1;
}

int gaussSeidel(double** A, double* b, double* x, double tolerancia, unsigned long int iMax, int n)
{
    int i, j;
    double erro = tolerancia + 1, *xant, soma;
    unsigned long int k = 0;
    
    if(tolerancia == 0)
	tolerancia = DBL_EPSILON;

    xant = criaVetor(n);

    for(i = 0; i < n; i++)
    {
        /*
        Se algum elemento da diagonal principal for zero, não haverá solução,
        a menos que as linhas da matriz sejam permutadas.
        */
        if(A[i][i] == 0)
	{
	    liberaVetor(xant);
            return 0;
	}

        //primeiro vetor para solução do método iterativo
        x[i] = b[i] / A[i][i];
    }

    while(k < iMax && erro >= tolerancia)
    {

        erro = 0;

        for(i = 0; i < n; i++)
            xant[i] = x[i];

        for(i = 0; i < n; i++)
        {
            soma = 0;

            for(j = 0; j < n; j++)
                if(j != i)
                    soma += A[i][j] * x[j];

            x[i] = (b[i] - soma) / A[i][i];

            if(erroRelativo(x[i], xant[i]) > erro)
                erro = erroRelativo(x[i], xant[i]);
        }

        k++;
    }
    
    liberaVetor(xant);

    /*
    Se erro for maior que a tolerancia no final das iterações, significa que
    não ocorreu convergência para o número máximo de iterações definidas.
    */
    if(erro >= tolerancia)
        return 0;
        
    for(i = 0; i < n; i++)
	    if(isnan(x[i]) || isinf(x[i]))
	        return 0;

    return 1;
}

int SOR(double** A, double* b, double* x, double tolerancia, double w, unsigned long int iMax, int n)
{
    int i, j;
    double erro = tolerancia + 1, *xant, soma;
    unsigned long int k = 0;
    
    if(tolerancia == 0)
	tolerancia = DBL_EPSILON;

    xant = criaVetor(n);

    for(i = 0; i < n; i++)
    {
        /*
        Se algum elemento da diagonal principal for zero, não haverá solução,
        a menos que as linha da matriz sejam permutadas.
        */
        if(A[i][i] == 0)
	{
	    liberaVetor(xant);
            return 0;
	}

        //primeiro vetor para solução do método iterativo
        x[i] = b[i] / A[i][i];
    }

    while(k < iMax && erro >= tolerancia)
    {

        erro = 0;

        for(i = 0; i < n; i++)
            xant[i] = x[i];

        for(i = 0; i < n; i++)
        {
            soma = 0;

            for(j = 0; j < n; j++)
                if(j != i)
                    soma += A[i][j] * x[j];

            x[i] = ((1 - w) * x[i]) + (w * (b[i] - soma) / A[i][i]);

            if(erroRelativo(x[i], xant[i]) > erro)
                erro = erroRelativo(x[i], xant[i]);
        }

        k++;
    }
    
    liberaVetor(xant);

    /*
    Se erro for maior que a tolerancia no final das iterações, significa que
    não ocorreu convergência para o número máximo de iterações definidas.
    */
    if(erro >= tolerancia)
        return 0;
        
    for(i = 0; i < n; i++)
	    if(isnan(x[i]) || isinf(x[i]))
	        return 0;

    return 1;
}

double criterioLinhas(double** A, int n)
{
    double soma;
    int i, j;
    
    /*
     * Para que o critério das linhas seja satisfeito, o valor absoluto do elemento de uma linha i
     * pertencente à diagonal principal deve ser maior do que a soma dos valores absolutos de cada
     * elemento dessa mesma linha.
    */

    for(i = 0; i < n; i++)
    {
        soma = 0;
        for(j = 0; j < n; j++)
            if(i != j)
                soma += fabs(A[i][j]);

        if(soma >= fabs(A[i][i]))
            return 0;
    }

    return soma;
}

double criterioColunas(double** A, int n)
{
    double soma;
    int i, j;
    
    /*
     * Para que o critério das colunas seja satisfeito, a soma da divisão entre cada elemento da linha i 
     * de uma coluna e o elemento da diagonal principal que pertence à linha i deverá ser menor que 1.
    */

    for(j = 0; j < n; j++)
    {
        soma = 0;
        for(i = 0; i < n; i++)
            if(i != j)
                soma += (fabs(A[i][j]) / fabs(A[i][i]));

        if(soma >= 1)
            return 0;
    }

    return soma;
}

double criterioNorma(double** A, int n)
{
    double norma;
    
    //Se a norma de Frobenius da matriz A for menor que 1, o critério será satisfeito.

    norma = normaF(A,n);

    if(norma >= 1)
        return 0;

    return norma;
}

double criterioSassenfeld(double** A, int n)
{
    double *b, maior;
    int i, j;

    b = criaVetor(n);

    for(i = 0; i < n; i++)
    {
        b[i] = 0;
        for(j = 0; j < n; j++)
        {
            if(j < i)
            {
                b[i] += fabs(A[i][j]) * b[j];
            }
            else if(j > i)
            {
                b[i] += fabs(A[i][j]);
            }
        }
        b[i] /= A[i][i];
    }
    
    /*
    Busca o maior elemento do vetor b calculado acima. Caso esse elemento seja 
    menor que 1, o critério será satisfeito.
    */
    maior = normaInfinito(b, n);
    
    liberaVetor(b);

    if(maior >= 1)
        return 0;

    return maior;
}

int cholesky(double **A, double **R, int n)
{
    int i, j, k;
    
    //R será uma matriz triangular superior caso A seja definida positiva e simétrica.

    for(i = 0; i < n; i++)
    {
        R[i][i] = 0;

	//Em cada iteração i será calculado o elemento da diagonal principal pertencente à linha i.
	/*
	O loop abaixo calcula o somatório do quadrado dos elementos que estão acima 
	do elemento da diagonal principal pertencente à linha i.
	*/
        for(k = 0; k < i; k++)
            R[i][i] += R[k][i] * R[k][i];

        R[i][i] = A[i][i] - R[i][i];

        if(R[i][i] <= 0)
            return 0;

        R[i][i] = sqrt(R[i][i]);

        for(j = i+1; j < n; j++)
        {
            R[i][j] = 0;
	    /*
	     * O loop abaixo calcula o somatório da multiplicação entre um elemento pertencente à coluna i
	     * e outro pertencente à coluna j, sendo que ambos são pertencentes à uma linha k menor que i.
	    */
            for(k = 0; k < i; k++)
                R[i][j] += R[k][i] * R[k][j];

            R[i][j] = A[i][j] - R[i][j];
            R[i][j] /= R[i][i];
        }
    }

    return 1;
}

int gradienteConjugado(double** A, double* b, double* x, double tol, unsigned long int iMax, int n)
{
    double *x0, *r, *d, *c, erro = tol + 1, q, p;
    int i;
    unsigned long int k = 0;
    
    x0 = criaVetor(n);
    r = criaVetor(n);
    d = criaVetor(n);
    c = criaVetor(n);
    
    if(tol == 0)
	tol = DBL_EPSILON;
    
    for(i = 0; i < n; i++)
	x[i] = b[i] / A[i][i];
    
    multiplicaVetor(A, n, n, x, r, n); // r =  Ax
    subtraiVetores(r, b, r, n); //r = Ax - b
    
    for(i = 0; i < n; i++)
	d[i] = r[i];
    
    while(erro >= tol && k < iMax)
    {
	
	multiplicaVetor(A, n, n, d, c, n); // c = Ad
	q = (produtoEscalar(r, d, n) * (-1)) / produtoEscalar(d, c, n); // q = -<r, d> / <d, Ad>
	
	for(i = 0; i < n; i++){ // x = x0 + q*d
	    x0[i] = x[i];
	    x[i] = x0[i] + q*d[i];
	}
	
	multiplicaVetor(A, n, n, x, r, n); // r =  Ax
	subtraiVetores(r, b, r, n); //r = Ax - b
	
	p = (produtoEscalar(r, c, n) * (-1)) / produtoEscalar(d, c, n); // p = -<r, Ad> / <d, Ad>
	
	for(i = 0; i < n; i++) // d = r + p*d
	    d[i] = r[i] + p*d[i];
	
	subtraiVetores(x, x0, x0, n); // x0 = x - x0
	erro = normaDois(x0, n) / normaDois(x, n); // erro = |x - x0| / |x| 
	
	k++;
    }
    
    if(erro >= tol)
	return 0;
    
    liberaVetor(x0);
    liberaVetor(r);
    liberaVetor(d);
    liberaVetor(c);
    
    return 1;
}

int newton(double *polinomio, double* raiz, double tolerancia, unsigned long int iMax, int n)
{
    int i = 0;
    double x0, x1, *deriv, erro = tolerancia + 1;
    unsigned long int k = 0;
    
    if(tolerancia == 0)
	tolerancia = DBL_EPSILON * 10;
    
    deriv = criaVetor(n-1);
    
    x0 = 1;
    
    derivada(polinomio, deriv, n);
    
    while(erro > tolerancia && k < iMax)
    {
	double p, d;
	
	//Calcula o valor de p no ponto (x0, p) do polinômio dado
	p = ordenada(polinomio, x0, n);
	
	//Calcula a derivada no ponto (x0, p) do polinômio dado
	d = ordenada(deriv, x0, n-1);
	
	//Calcula novo valor de x com base na iteração anterior, ou no x0 inicial
	x1 = x0 - (p / d);
	    
	erro = fabs(x1 - x0);
	
	//x0 recebe x1, pois na proxima iteração, x1 atual será o x da iteração anterior (x0)
	x0 = x1;
	
	k++;
    }
    
    liberaVetor(deriv);
    
    if(k >= iMax)
	return 0;
    
    *raiz = x1;
    
    return 1; 
}

void schmidt(double** vetores, double** base, int dimensao)
{
    int i, j;
    double *proj, norma;
    
    proj = criaVetor(dimensao);
    
    for(i = 0; i < dimensao; i++)
    {
	atribui(base[i], vetores[i], dimensao);
	
	for(j = 0; j < i; j++)
	{
	    projecaoOrtogonal(vetores[i], base[j], proj, dimensao);
	    subtraiVetores(base[i], proj, base[i], dimensao);
	}
	
	norma = normaDois(base[i], dimensao);
	multiplicaPorEscalar(base[i], base[i], 1/norma, dimensao);
    }
    
    liberaVetor(proj);
}

void schmidtModificado(double** vetores, double** base, int dimensao)
{
    int i, j, t;
    double *proj, norma;
    
    proj = criaVetor(dimensao);
    
    for(i = 0; i < dimensao; i++)
    {
	atribui(base[i], vetores[i], dimensao);
	
	for(j = 0; j < i; j++)
	{
	    projecaoOrtogonal(base[i], base[j], proj, dimensao);
	    
	    subtraiVetores(base[i], proj, base[i], dimensao);
	}
	
	norma = normaDois(base[i], dimensao);
	multiplicaPorEscalar(base[i], base[i], 1/norma, dimensao);
    }
    
    liberaVetor(proj);
}

double potencia(int n, double tol, double **A, double *x) {
  //Definindo as variaveis que vao ser usadas
  double *q,
          sigma = 0,
          sigmaAux,
          *vetorAux,
          limite;

  int i = 0,
      j = 0,
      k = 0,
      l = 0,
      FLAG = 0;

  // Alocando os vetores
  q = criaVetor(n);
  vetorAux = criaVetor(n);
  
  //Definindo os valores iniciais das variaveis
  sigmaAux = normaInfinito(x, n);
  divideVetorPorEscalar(x, n, sigmaAux);

  //LOOP DO METODO
  while(FLAG == 0){
    //A*x = Aqi
    multiplicaVetor(A, n, n, x, q, 1);

    //Achando o sigma usando a norma infinito
    sigmaAux = normaInfinito(q, n);

    //Usei um vetor auxiliar para verificar o limite para a tolerancia do metodo
    for( l = 0; l< n; l++) 
      vetorAux[l] = x[l] - q[l]/sigmaAux;

    limite = normaInfinito(vetorAux, n);    

    //Transformando o vetor x em qi
    for(i = 0; i < n; i++)
      x[i] = q[i]/sigmaAux;

    //Verificacao da tolerancia
    if(fabs(limite) < tol) 
      FLAG = 1;

    //define sigma como o sigma auxiliar
    sigma = sigmaAux;

    //chama a funcao para imprimir
    imprimePotencia(x, sigma, n, k);
    k++;
  }
  //liberando o espaço alocado
  liberaVetor(q);
  liberaVetor(vetorAux);
  //o ultimo sigma vai ser o autovalor dominante
  return sigma;
}

/**
 * Função que retorna um número entre 0 e 1.
 */
double aleatorio()
{
    return (double)rand() / (double)RAND_MAX ;
}

/**
 * Função para perturbar uma solução Simulated Annealing
 */
void perturba(double* x, double* xp, double intervalo, int n)
{
    double rnd, rnd0;
    int i;
    
    //Número no intervalo [-intervalo/2, intervalo/2]
    rnd0 = aleatorio() * intervalo - (intervalo / 2);
    rnd = rnd0;
    
    for(i = 0; i < n; i++)
    {
	while(rnd == rnd0)
	    rnd = aleatorio() * intervalo - (intervalo / 2);
	
	rnd0 = rnd;
	
	xp[i] = x[i] + rnd;
    }
}

void simulatedAnnealing(double** A, double* b, double* x, double alfa, double ls, double li, int nr, int n)
{
    double *xp, *r, *rp, rnd0, rnd = 0, normaRp, normaR, variacao, intervalo = 2;
    int i, j;
    
    rnd0 = aleatorio();
    rnd = rnd0;
    
    xp = criaVetor(n);
    r = criaVetor(n);
    rp = criaVetor(n);
    
    //Busca uma solução inicial para x
    perturba(x, x, 130, n);
    
    while(ls > li)
    {
	for(i = 0; i < nr; i++)
	{
	    perturba(x, xp, intervalo, n);
	    
	    residuo(A, x, b, r, n);
	    residuo(A, xp, b, rp, n);
	    
	    normaRp = normaDois(rp, n);
	    normaR = normaDois(r, n);
	    
	    if(normaRp < 1)
		    intervalo = 0.2;
	    
	    if(normaRp <= normaR)
	    {
		atribui(x, xp, n);
	    }
	    else
	    {
		while(rnd == rnd0)
		    rnd = aleatorio();
		
		rnd0 = rnd;
		
		variacao = normaRp - normaR;
		
		if(rnd < exp(-variacao / ls))
		{
		    atribui(x, xp, n);
		}
	    }
	}
	
	ls *= alfa;
    }
    
    liberaVetor(xp);
    liberaVetor(r);
    liberaVetor(rp);
}

void PSO(double **matrizA, int n, double *vetorB, double *MB, double tol, int particulas){ /// MB é a solução do sistema

    double **MP, **S, **V, *vResiduo, normaMB, *vetMB, *getVetV, *getVetMP, *getVetS, *velAtual, *vetResiduoMPCalcu, *vetResiduoMPAtual;
    int i, j;
    double c1 = 2, c2 = 2, A1, A2, A3;

    vetMB              = criaVetor(n);
    vResiduo           = criaVetor(n);
    getVetV            = criaVetor(n);
    getVetMP           = criaVetor(n);
    getVetS            = criaVetor(n);
    velAtual           = criaVetor(n);
    vetResiduoMPCalcu  = criaVetor(n);
    vetResiduoMPAtual  = criaVetor(n);

    S  = criaMatriz(n, particulas);
    V  = criaMatriz(n, particulas);
    MP = criaMatriz(n, particulas);


    /// Inicializa os vetores das particulas
    for(i = 0; i < n; i++) /// cada coluna dessa matriz vai ser uma particula
    {
        for(j = 0; j < particulas; j++)
        {
            S[i][j] = aleatorio() * 10;
            V[i][j] = aleatorio() * 10;
            MP[i][j] = S[i][j]; /// a melhor posição vai ser a primeira posição
        }
    }


    /// encontrar o primeiro MB usando a primeira particula como base
    getColuna(S, MB, 0, n);
    residuo(matrizA, MB, vetorB, vResiduo, n);
    normaMB = normaDois(vResiduo, n);


    for(i = 0; i < particulas; i++) /// Comparar com as outras particulas para saber qual é a melhor do bando
    {
	getColuna(S, vetMB, i, n);
        residuo(matrizA, vetMB, vetorB, vResiduo, n);
        if(normaMB > normaDois(vResiduo,n))
        {
	    for(j = 0; j < n; j++)
            {
                MB[j] = vetMB[j];
	    }
            normaMB = normaDois(vResiduo,n);
        }
    } /// Já tenho o vetor MB

    while(tol < normaMB)
    {
        for(j = 0; j < particulas; j++)
        {
            A1 = aleatorio();

            A2 = aleatorio();

            A3 = aleatorio();

            getColuna(V, getVetV, j, n);
            getColuna(S, getVetS, j, n);
            getColuna(MP, getVetMP, j, n);

            atualizaVel(A1, A2, A3, c1, c2, n, getVetMP, getVetV, getVetS, MB, velAtual);
            preencheColuna(V, j, n, velAtual); ///Atualiza os vetores de velocidade de cada particula. Já que cada coluna é uma particula, o preencheColuna vai colocar o vetor com a velocidade atualizada em cada particula

            atualizaPos(getVetS, velAtual, n); /// atualiza a posicao de cada particula
            preencheColuna(S, j, n, getVetS); /// preenche a particula atualizada

            ///Atualiza o MP
            residuo(matrizA, getVetS, vetorB, vetResiduoMPCalcu, n);
            residuo(matrizA, getVetMP, vetorB, vetResiduoMPAtual, n);
	    
            if(normaDois(vetResiduoMPCalcu, n) < normaDois(vetResiduoMPAtual, n)) /// se a norma do residuo da nova posicao for menor então substitui na melhor posicao
            {
                preencheColuna(MP, j, n, getVetS);
            }
        }
            /// pegar o MB
        for(i = 0; i < particulas; i++) /// Comparar com as outras particulas para saber qual é a melhor do bando
        {
            getColuna(MP, vetMB, i, n);
            residuo(matrizA, vetMB, vetorB, vResiduo, n);
	    
            if(normaMB > normaDois(vResiduo, n))
            {
                for(j = 0; j < n; j++)
                {
                    MB[j] = vetMB[j];
                }
                
                normaMB = normaDois(vResiduo,n);
            }
        }
    }

    liberaMatriz(S,n);
    liberaMatriz(V,n);
    liberaMatriz(MP,n);

    liberaVetor(vetMB);
    liberaVetor(vResiduo);
    liberaVetor(getVetV);
    liberaVetor(getVetMP);
    liberaVetor(getVetS);
    liberaVetor(velAtual);
    liberaVetor(vetResiduoMPCalcu);
    liberaVetor(vetResiduoMPAtual);
}