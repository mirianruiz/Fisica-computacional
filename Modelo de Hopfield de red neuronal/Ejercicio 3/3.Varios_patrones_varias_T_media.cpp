#include <iostream>
#include <vector>
#include <math.h> 
#include <algorithm> 
#include <fstream> 
#include <string>
#include <sstream>
#include <chrono>
#include <random> 

using namespace std;


vector<vector<vector<vector<double>>>> itereaccion_particulas_un_patron(int num_p, vector<vector<vector<int>>> patrones_iniciales, int N);
double algoritmo_metropolis(double T, mt19937_64 & generator, vector<vector<vector<vector<double>>>> & w, vector<vector<int>> S, ofstream & fich_neuronas, int N);
vector<vector<int>> generar_patron_deformado(vector<vector<vector<int>>> patron_inicial, mt19937_64 & generator, int N);

int main() {
    // Numero de patrones iniciales y tamaño de la red    
    int num_p = 3;
    int N = 30;
    
    // Leo los patrones iniciales y los introduzco en matrices
    ifstream fich_cero("0_m.txt");
    ifstream fich_uno("1_m.txt");
    ifstream fich_dos("2_m.txt");

    // Genero un vector de patrones iniciales
    vector<vector<vector<int>>> patrones_iniciales(num_p, vector<vector<int>>(N, vector<int>(N, 0)));

    // Relleno el array con los patrones iniciales
    for (int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            fich_cero >> patrones_iniciales[0][i][j];
            fich_uno >> patrones_iniciales[1][i][j];
            fich_dos >> patrones_iniciales[2][i][j];
        }
    }

    fich_cero.close();
    fich_uno.close();
    fich_dos.close();

    // Genero una semilla para el generador de números aleatorios a partir del reloj
    unsigned seed1 = chrono::system_clock::now().time_since_epoch().count(); 
    mt19937_64 generator(seed1);

    // Calculo la interacción entre las partículas para los patrones iniciales
    vector<vector<vector<vector<double>>>> w = itereaccion_particulas_un_patron(num_p, patrones_iniciales, N);
    
    // Genero un array de temperaturas
    vector<double> Temperaturas = { 0.01,0.02,0.03,0.04,0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2 };
    
    // NÚMERO DE EXPERIMENTOS
    int num_exp = 10; 

    // Bucle principal por temperatura
    for (double T : Temperaturas) {
        
        // Preparar archivos para esta temperatura
        ostringstream nombre_aleatorio;
        nombre_aleatorio << "f_neuronas_aleatorios_varios_patrones_T_" << T << ".txt";
        ofstream fich_aleatorio(nombre_aleatorio.str());

        ostringstream nombre_deformado;
        nombre_deformado << "f_deformado_varios_patrones_T_" << T << ".txt";
        ofstream fich_deformado(nombre_deformado.str());

        // Repetir el experimento 10 veces
        for (int exp = 0; exp < num_exp; exp++) {
            
            // 1. Generar un patrón aleatorio NUEVO para este experimento específico
            vector<vector<int>> S_aleatorio(N, vector<int>(N));
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    uniform_int_distribution<int> dist(0, 1);
                    S_aleatorio[i][j] = dist(generator);
                }
            }

            // 2. Generar un patrón deformado NUEVO para este experimento específico
            vector<vector<int>> patron_deformado_exp = generar_patron_deformado(patrones_iniciales, generator, N);

            // 3. Ejecutar metrópolis y guardar SOLAMENTE la última matriz en el archivo
            algoritmo_metropolis(T, generator, w, S_aleatorio, fich_aleatorio, N);
            algoritmo_metropolis(T, generator, w, patron_deformado_exp, fich_deformado, N);
        }
        
        fich_aleatorio.close();
        fich_deformado.close();
    }

    return 0;
}


// Función para calcular las interacciones entre partículas 
vector<vector<vector<vector<double>>>> itereaccion_particulas_un_patron(int num_p, vector<vector<vector<int>>> patrones_iniciales, int N) {
    vector <double> a(num_p, 0.0);
    
    for (int mu = 0; mu < num_p; mu++){
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
                    if (i != k && j != l) {  
                        for (int mu = 0; mu < num_p; mu++){
                            w[i][j][k][l] += (patrones_iniciales[mu][i][j]-a[mu]) * (patrones_iniciales[mu][k][l]-a[mu])*1.0/(N*N);
                        }
                    }
                }
            }
        }
   }
   return w;
}


// Función del algoritmo de Metrópolis
double algoritmo_metropolis(double T, mt19937_64 & generator, vector<vector<vector<vector<double>>>> & w, vector<vector<int>> S, ofstream & fich_neuronas, int N) {

    int Paso_montecarlo = N * N;

    // Ejecutar los 20 pasos Monte Carlo completos
    for (int k = 0; k < 20 * Paso_montecarlo; k++) { 

        // Paso 1: Escoger una partícula al azar
        uniform_int_distribution<int> distribution(0, N-1);
        int i = distribution(generator);
        int j = distribution(generator);

        // Paso 2: Diferencia de energía
        double delta_E = 0.0;

        if (S[i][j] == 1) {
            for (int r = 0; r < N; r++) { // Cambiado k por r para no pisar la k del bucle principal
               for (int l = 0; l < N; l++) {
                        delta_E += w[i][j][r][l] * (S[r][l] - 0.5);
                }
            }
        } else {
            for (int r = 0; r < N; r++) {
                for (int l = 0; l < N; l++) {
                        delta_E += w[i][j][r][l] * (0.5 - S[r][l]);
                }
            }
        }
            
        double p = min(1.0, exp(-delta_E / T));

        // Paso 3: Decidir si se acepta el cambio
        uniform_real_distribution<double> r_distribution(0.0, 1.0);
        double rand_val = r_distribution(generator);

        if (rand_val < p) { 
            S[i][j] = 1 - S[i][j]; 
        }

        // Paso 4: Guardar la matriz SOLAMENTE al finalizar todo el proceso (en el último intento del último paso Monte Carlo)
        if (k == 20 * Paso_montecarlo - 1) {
            for (int fila = 0; fila < N; fila++) {
                for (int col = 0; col < N; col++) {
                    fich_neuronas << S[fila][col];
                    if (col < N - 1) {
                        fich_neuronas << ' ';
                    }
                }
                fich_neuronas << '\n';
            }
            fich_neuronas << '\n'; // Línea en blanco para separar los resultados de cada experimento
        }
    }

    return 0;
}


// Función para deformar el patrón inicial
vector<vector<int>> generar_patron_deformado(vector<vector<vector<int>>> patron_inicial, mt19937_64 & generator, int N) {
    vector<vector<int>> patron_deformado(N, vector<int>(N));
    patron_deformado = patron_inicial[2]; 

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