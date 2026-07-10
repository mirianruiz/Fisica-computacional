
#include <iostream>
#include <vector>
#include <math.h> // Para usar la función sqrt
#include <algorithm> // Para usar la función min
#include <fstream> // Para usar la función ofstream
#include <string>
#include<chrono>
#include <random> 
using namespace std;


vector<vector<vector<vector<double>>>> itereaccion_particulas(int num_p, vector<vector<vector<int>>> epsilon, int N);
double algoritmo_metropolis(double T, mt19937_64 & generator, vector<vector<vector<vector<double>>>> & w, vector<vector<int>> & S, ofstream & fich_neuronas, int N);
vector<vector<int>> generar_patron_deformado(vector<vector<vector<int>>> epsilon, mt19937_64 & generator, int N);


int main(){

    //Numero de patrones iniciales    
     int num_p=3;
     int N=30;

    
    //Leo los patrones iniciales y los introduzco en matrices

    ifstream fich_cero("0_m.txt");
    ifstream fich_uno("1_m.txt");
    ifstream fich_dos("2_m.txt");

    // genero un vector de patrones iniciales, cada elemento es una matriz de ceros inicialemnte, y cada matriz va a ser un patrón inicial diferente

    vector<vector<vector<int>>> patrones_iniciales(num_p, vector<vector<int>>(N, vector<int>(N, 0)));

    //Relleno el array con los patrones iniciales

    for (int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            fich_cero>>patrones_iniciales[0][i][j];
            fich_uno>>patrones_iniciales[1][i][j];
            fich_dos>>patrones_iniciales[2][i][j];
            

        }

    }

    fich_cero.close();
    fich_dos.close();
    fich_uno.close();


    //****************************Genero un patrón deformado*****************

    //Genero una semilla para el generador de números aleatorios, que se puede elegir como un número grande o se puede obtener del reloj del sistema para que sea diferente cada vez que se ejecute el programa
    unsigned seed1 = chrono::system_clock::now().time_since_epoch().count(); // obtain a seed from the system clock ## decide a rule to start reading random numbers from a table
    //Inicializo el generador de números aleatorios con la semilla obtenida
    mt19937_64 generator(seed1);
    vector<vector<int>> patron_deformado = generar_patron_deformado(patrones_iniciales, generator, N);

    //************************Genero un patron aleatorio**********************/
    vector<vector<int>> S(N, vector<int>(N));
        // Inicializar la matriz S con valores aleatorios de 0 y 1
         // Inicializar la matriz S con valores aleatorios de 0 y 1
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                uniform_int_distribution<int> dist(0, 1);
                S[i][j] = dist(generator);
            }
        }
    

    //Calculo la interacción entre las partículas para el patrón inicial
    vector<vector<vector<vector<double>>>> w = itereaccion_particulas(num_p, patrones_iniciales, N);

    

    //Aplico el algoritmo de Metropolis para cada temperatura T
    double T = 1e-4;


    //Guardo los estados de las neuronas generados de forma aleatoria para cada paso motecarlo

    ofstream fich_aleatorio("neuronas_aleatorios_varios_patrones_m.txt");

    //Llamo al algoritmo de Metropolis para el S aleatorio

    algoritmo_metropolis(T, generator, w, S, fich_aleatorio, N);


    
    //Guardo los estados de las neuronas generados a partir del patrón deformado
    ofstream fich_deformado("deformado_varios_patrones_m.txt");

    //Llamo al algoritmo de Metropolis para el patrón deformado
    algoritmo_metropolis(T, generator, w, patron_deformado, fich_deformado, N);


    return 0;
}

//************Genero una función para calcular las interacciones entre partículas 

vector<vector<vector<vector<double>>>> itereaccion_particulas(int num_p, vector<vector<vector<int>>> patrones_iniciales, int N) {
    vector <double> a(num_p, 0.0);
    
    for (int mu=0; mu<num_p; mu++){
        // Calcular a_elem como la sumatoria de todos los elementos dividida por N²
        double suma = 0.0;
        for (int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                suma += patrones_iniciales[mu][i][j];
            }
        }
        a[mu] = (1.0/(N*N)) * suma;

    }
    
    
    // Inicializar la matriz w
    vector<vector<vector<vector<double>>>> w(N, vector<vector<vector<double>>>(N, vector<vector<double>>(N, vector<double>(N, 0.0))));
    
        for (int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                for(int k = 0; k < N; k++){
                    for(int l = 0; l < N; l++){
                        if (i != k && j != l) {  // Si (i,j) es distinto de (k,l)
                            for (int mu=0; mu<num_p; mu++){
                            w[i][j][k][l] += (patrones_iniciales[mu][i][j]-a[mu]) * (patrones_iniciales[mu][k][l]-a[mu])*1.0/(N*N);
                            }
                        }
                    }
                }
            }
    
       }
    return w;
    }


//**************Genero una función para que para cada T

//Tengo que meterle el generador de números aleatorios

double algoritmo_metropolis(double T, mt19937_64 & generator, vector<vector<vector<vector<double>>>> & w, vector<vector<int>> & S, ofstream & fich_neuronas, int N) {

    int Paso_montecarlo = N * N;

    for (int fila = 0; fila < N; fila++) {
        for (int col = 0; col < N; col++) {
            fich_neuronas << S[fila][col];
            if (col < N - 1) {
                fich_neuronas << ' ';
            }
        }
        fich_neuronas << '\n';
    }
    fich_neuronas << '\n';
    


    for (int k = 0; k < 20 * Paso_montecarlo-1; k++) { //Repetir los pasos 1 a 3 un número suficiente de veces

    
    //===============================Paso 1===========================================

    //Escoger un elemento aleatorio de la matriz S, es decir, escoger una partícula al azar en la red
        
        uniform_int_distribution<int> distribution(0, N-1);
        int i = distribution(generator);
        int j = distribution(generator);


    //===============================Paso 2===========================================

    //Diferencia de energía al cambiar el espín del elemento escogido, es decir, la energía que se gana o se pierde al cambiar el espín del elemento escogido
        
        double delta_E = 0.0;

        if (S[i][j] == 1) {
            for (int k = 0; k < N; k++) {
               for (int l = 0; l < N; l++) {
                        delta_E += w[i][j][k][l] * (S[k][l] - 0.5);
                
                }
            }
        } else {
            for (int k = 0; k < N; k++) {
                for (int l = 0; l < N; l++) {
                    
                        delta_E += w[i][j][k][l] * (0.5 - S[k][l]);
                }
            }
        }
    

            
        double p=min(1.0, exp(-delta_E / T));

    //===============================Paso 3===========================================

    //Generar un aleatorio entre 0 y 1 para decidir si se acepta el cambio de espín o no

        //r = np.random.rand() //Genera un número aleatorio entre 0 y 1
        uniform_real_distribution<double> r_distribution(0.0, 1.0);
        double r = r_distribution(generator);

        if (r < p) { //Si el número aleatorio es menor que la probabilidad de aceptar el cambio de espín, entonces se acepta el cambio de espín

            S[i][j]=1-S[i][j]; //Cambiar el espín del elemento escogido, es decir, cambiar el valor de S[i][j] de 0 a 1 o de 1 a 0
    
        }

    //===============================Paso 4===========================================
    //Repetir los pasos 1 a 3 un número suficiente de veces para que el sistema alcance el equilibrio térmico, es decir, para que el sistema alcance un estado estable donde las propiedades macroscópicas no cambien con el tiempo

    //Guardar los valores de la matriz S en cada iteración para crear una animación de la evolución del sistema

        // Guardar solo al terminar un paso Monte Carlo completo (N*N intentos)
        if ((k + 1) % Paso_montecarlo == 0) {
            for (int fila = 0; fila < N; fila++) {
                for (int col = 0; col < N; col++) {
                    fich_neuronas << S[fila][col];
                    if (col < N - 1) {
                        fich_neuronas << ' ';
                    }
                }
                fich_neuronas << '\n';
            }
            fich_neuronas << '\n'; // Línea en blanco para separar cada paso Monte Carlo
        }
    }


    return 0;
}



//***************Genero una función para deformar el patrón iniciaL, uno de los patrones

vector<vector<int>> generar_patron_deformado(vector<vector<vector<int>>> patron_inicial, mt19937_64 & generator, int N) {
    vector<vector<int>> patron_deformado(N, vector<int>(N));
    patron_deformado = patron_inicial[2]; // Inicializar el patrón deformado con el patrón inicial

    for (int k = 0; k < N * N / 10; k++) {

            uniform_int_distribution<int> distribution(0, N-1);
            int i = distribution(generator);
            int j = distribution(generator);


            if (patron_inicial[2][i][j] == 1) {
                patron_deformado[i][j] = 0;
            } else {
                patron_deformado[i][j] = 1;
            }
        }
    

    return patron_deformado;

    }

