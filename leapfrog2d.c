#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#define nx_ 71
#define ny_ 71

double gauss( int x, int y, int tam);

double     H( int x, int y, int tam);

void atualizar_onda(double u_new[][ny_], double u_now[][ny_], double u_old[][ny_], double lambda[][ny_], double a,  double b,  double c, int nx, int ny, double rx, double ry);

double delta_u( double rx, double ry, double lambda[][ny_],  double u_now[][ny_], int i,  int j, int  im1, int ip1, int jm1, int jp1);



/*        META DADOS   */
//
/*    nx_ e ny_ , definidos acima por praticidade,    */
/*    setam o tamanho dos vetores.                    */
/*                                                    */
/*    tmax   , tempo total da simulacao               */
/*    rx e ry, dt/dx e dy/dt respectivamente          */
/*    lambda , mapa do solo (profundidade em (x,y))   */
//
//
/*     u representa a altura do mar em (x,y), sendo:  */
//
/*    u_old[x][y] em t = t-1    */
/*    u_now[x][y] em t = t      */
/*    u_new[x][y] em t = t+1    */
//



void main()
{

    FILE *arq;
    arq = fopen("onda.txt", "w+");


    int tmax, i, j, t, nx, ny;

    nx = nx_;
    ny = ny_;
    tmax = 300;

    double u_new[nx][ny], u_old[nx][ny], u_now[nx][ny], lambda[nx][ny], rx, ry;

    rx = 0.25;
    ry = 0.25;

    //   incio loop da condição inicial,
    //   laco duplo e usado para percorrer matrix nx*ny

    for(j = 0 ; j < ny  ; j++)
    {
        for(i = 0 ; i < nx ; i++)
        {
            /* u_now inicial em forma de sino */
            u_now[i][j] = gauss(i,j,nx);
            /* lambda inicial em forma de H() */
            lambda[i][j] = H(i,j,nx);

        }
    }

    //   fim loop da condição inicial


    /*  calculo "sintetico" de u_old, pois a atualizacao */
    /*  nao e auto-iniciavel   */
    atualizar_onda(u_old, u_now, u_old, lambda, 0.5, 0, 0.5, nx, ny, rx, ry);




    //inicio do laço temporal

    for(t = 0 ; t < tmax ; t++)
    {

        //imprimindo no arquivo com laco dupl


        /* calculam-se novos valores */

        atualizar_onda( u_new, u_now, u_old, lambda, 1, 1, 1, nx, ny, rx, ry);


        /* matrizes antes novas se tornam antigas */
        /*   u_old = u_now, u_now = u_new   */


        memcpy(u_old,u_now, sizeof(double)*nx*ny);
        memcpy(u_now,u_new, sizeof(double)*nx*ny);


        /* OBS: Double tem o balor de 8 bytes na memoria. Como temos uma matriz de nx*ny, pegamos o tamanho*/
        /* de um double e multiplicamos pela dimensao da matriz */
        /* sintaxe memcpy(matriz a ser atualizada, matriz que passa o valor, tamanho em bytes da matriz)  */

    }

            for(j = 0 ; j < ny ; j++)
        {
            for(i = 0 ; i < nx  ; i++)
            {
                if(i == j){
                fprintf(arq, "%d %d %lf\n", i, j, u_now[i][j]);
                }
            }
        }

        /* Linhas em branco ao final de cada tempo para index do gnuplot */
        fprintf(arq,"\n\n");

    fclose(arq);
}






double gauss(int x, int y,int tam)
{
    /* curva inicial em forma de sino*/

    double A, xc, yc, gauss, sigmax, sigmay;
    //xc = (tam-1) / 2.;
    //yc = (tam-1) / 2.;
    xc = 0;
    yc = 0;

    A = 1;
    sigmax = 1;
    sigmay = 1;

    gauss = A * exp(-0.5*pow(((x-xc)/sigmax),2) -0.5*pow(((y-yc)/sigmay),2));

    return gauss;
}




/* a funcao H corresponde ao formato do terreno, retorna a profundidade em relação a aguas calmas */

double H(int x, int y,int tam)
{
    /* curva para o formato do terreno em forma de sino virado*/
    double A,xc,yc,h,sigmax,sigmay;
    //xc = (tam-1) / 2.;
    //yc = (tam-1) / 2.;

    A = 1;
    sigmax = 1;
    sigmay = 1;
    xc = 0;
    yc = 0;

    h = 1 - A * exp(-0.5*pow(((x-xc)/sigmax),2) -0.5*pow(((y-yc)/sigmay),2));

    return h;
}

void atualizar_onda(double u_new[][ny_], double u_now[][ny_], double u_old[][ny_], double lambda[][ny_], double a,  double b,  double c, int nx, int ny, double rx, double ry)
{
// DOUBLE U_NEW[][NY_)] SINALIZA PARA O C QUE A FUNÇÃO RECEBERA UMA MATRIZ( TECNICAMENTE RECEBERA O ENDEREÇO NA MEMORIA DA MATRIZ) POR ISSO NÃO É NECESSARIO RETORNAR NENHUM VALOR
// ESTA "TECNICA" É POSSIVEL POIS O NOME DA MATRIZ, NO CASO U_NEW, É NA VERDADE O ENDEREÇO DESSA MATRIZ NA MEMORIA. COMO ESTAMOS PASSANDO O ENDEREÇO NA MEMORIA, A FUNÇÃO CONSEGUE ALTERAR OS VALORES SEM NECESSIDADE DE RETORNO
// PRECISAMOS DECLARAR O [NY_] NA FUNÇÃO POR QUESTÕES TECNICAS. MAS PENSEM NESSA SINTAXE COMO: PASSANDO O ENDEREÇO DA MATRIX PARA A FUNÇÃO, SEM NECESSIDADE DE RETORNO. A FUNÇÃO É LIVRE PARA EDITAR A PROPRIA MATRIZ

    int i,j;

    // LOOP DUPLO PARA ATUALIZAR AS PARTES INTERNAS DA MATRIZ
    for(j = 1 ; j < ny - 1 ; j++)
    {
        for(i = 1 ; i < nx - 1  ; i++)

            // SEPAREI PARTE DA EQUAÇAO EM OUTRA FUNÇÃO PARA SIMPLIFICAR A VIDA. DELTA_U É SIMPLESMENTE UMA PARTE DA EQUAÇÃO ORIGINAL

            u_new[i][j] = a * 2 * u_now[i][j] - b * u_old[i][j] + c * delta_u(rx,ry,lambda,u_now,i,j,i-1,i+1,j-1,j+1);
    }

    //PRECISAMOS AGORA ATUALIZAR AS LINHAS E COLUNAS EXTERNAS DA MATRIZ, POIS ESTAS NÃO FORAM INCLUIDAS NO LOOP ANTERIOR. ESTA AÇÃO NÃO ATUALIZA AS PONTAS, OU QUINAS, DA MATRIZ


    i = 0; // ATUALZANDO A PRIMEIRA COLUNA DA MATRIZ
    for( j = 1; j < ny - 1  ; j++)
    {
        u_new[i][j] = a * 2 * u_now[i][j] - b * u_old[i][j] + c * delta_u(rx,ry,lambda,u_now,i,j,i+1,i+1,j-1,j+1);
    }

    i = nx - 1;// ATUALZANDO A ULTIMA COLUNA DA MATRIZ
    for( j = 1; j < ny - 1  ; j++)
    {
        u_new[i][j] = a * 2 * u_now[i][j] - b * u_old[i][j] + c * delta_u(rx,ry,lambda,u_now,i,j,i-1,i-1,j-1,j+1);
    }

    // ---------------------- //


    j = 0; // ATUALZANDO A PRIMEIRA LINHA DA MATRIZ
    for( i = 1; i < nx - 1 ; i++)
    {
        u_new[i][j] = a * 2 * u_now[i][j] - b * u_old[i][j] + c * delta_u(rx,ry,lambda,u_now,i,j,i-1,i+1,j+1,j+1);

    }


    j = ny - 1; // ATUALZANDO A ULTIMA LINHA DA MATRIZ
    for( i = 1; i < nx - 1 ; i++)
    {
        u_new[i][j] = a * 2 * u_now[i][j] - b * u_old[i][j] + c * delta_u(rx,ry,lambda,u_now,i,j,i-1,i+1,j-1,j-1); //possivel erro erro achado

    }



    // AGORA VAMOS ATUALIZAR as PONTAS DA MATRIZ. OU CANTOS, SE PREFERIR CHAMAR ASSIM


    // PONTA [0][0]
    i = 0;
    j = 0;
    u_new[i][j] = a * 2 * u_now[i][j] - b * u_old[i][j] + c * delta_u(rx,ry,lambda,u_now,i,j,i+1,i+1,j+1,j+1);


    // PONTA [nx - 1][0]
    i= nx - 1;
    j= 0;
    u_new[i][j] = a * 2 * u_now[i][j] - b * u_old[i][j] + c * delta_u(rx,ry,lambda,u_now,i,j,i-1,i-1,j+1,j+1);


    // PONTA [0][ny -1]
    i= 0;
    j= ny - 1;
    u_new[i][j] = a * 2 * u_now[i][j] - b * u_old[i][j] + c * delta_u(rx,ry,lambda,u_now,i,j,i+1,i+1,j-1,j-1);


    // PONTA [NX-1][ny -1] (rx,ry,lambda,u_now,i,j,i-1,i+1,j-1,j-1)
    i= nx - 1;
    j= ny - 1;
    u_new[i][j] = a * 2 * u_now[i][j] - b * u_old[i][j] + c * delta_u(rx,ry,lambda,u_now,i,j,i-1,i-1,j-1,j-1); // possivel erro




}

double delta_u( double rx, double ry, double lambda[][ny_],  double u_now[][ny_], int i,  int j, int  im1, int ip1, int jm1, int jp1)
{
    //CALCULAMOS AQUI  SEPARADO UMA PARTE DA EQUAÇÃO, POIS ELA MUDA DEPENDENDO DO QUE ESTAMOS CALCULANDO. SEJA UMA COLUNA INICIAL OU FINAL OU UM CANTO DA MATRIZ
    //PARA ISSO COLOQUEI VALORES AUXILIARES DE IM1, IP1, JM1, JP1 . P DE PLUS E M DE MINUS.POIS ESTES VALORES SAO TROCADOS EM VARIAS PARTES

    return (
               pow(rx,2) * ( ((0.5 * (lambda[ip1][j] + lambda[i][j])) * (u_now[ip1][j] - u_now[i][j]))
                             -  ((0.5 * (lambda[i][j] + lambda[im1][j])) * (u_now[i][j] - u_now[im1][j])))


               + pow(ry,2) * ( ((0.5 * (lambda[i][jp1] + lambda[i][j])) * (u_now[i][jp1] - u_now[i][j]))
                               -  ((0.5 * (lambda[i][j] + lambda[i][jm1])) * (u_now[i][j] - u_now[i][jm1]))));

}
