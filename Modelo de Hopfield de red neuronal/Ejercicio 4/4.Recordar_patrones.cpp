#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>
#include <chrono>
#include <random>

using namespace std;



//**************Función para calcular la iteracción de las partículas */

vector<vector<vector<vector<double>>>> itereaccion_particulas(int num_p, vector<vector<vector<int>>> patrones_iniciales, int N) {
    vector<double> a(num_p, 0.0);
    
    for (int mu = 0; mu < num_p; mu++){
        double suma = 0.0;
        for (int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                suma += patrones_iniciales[mu][i][j];
            }
        }
        a[mu] = (1.0 / (N * N)) * suma;
    }
    

    vector<vector<vector<vector<double>>>> w(N, vector<vector<vector<double>>>(N, vector<vector<double>>(N, vector<double>(N, 0.0))));
    
    for (int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            for(int k = 0; k < N; k++){
                for(int l = 0; l < N; l++){
                    if (i != k || j != l) {  
                        for (int mu = 0; mu < num_p; mu++){
                            w[i][j][k][l] += (patrones_iniciales[mu][i][j] - a[mu]) * (patrones_iniciales[mu][k][l] - a[mu]) * 1.0 / (N * N);
                        }
                    }
                }
            }
        }
    }
    return w;
}

//*********************Funcion para calcular el solapamiento*****************/

double calcular_solapamiento(const vector<vector<int>>& S, const vector<vector<int>>& patron, double a, int N) {
    double suma = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            suma += (patron[i][j] - a) * (S[i][j] - a);
        }
    }
    return suma / (N * N * a * (1.0 - a));
}

//***********************Función para el algoritmo de Metropolis*****************/

void metropolis_rapido(double T, mt19937_64 & generator, const vector<vector<vector<vector<double>>>> & w, vector<vector<int>> & S, int N) {
    int Paso_montecarlo = N * N;
    int pasos_totales = 20 * Paso_montecarlo; 

    uniform_int_distribution<int> distribution(0, N - 1);
    uniform_real_distribution<double> r_distribution(0.0, 1.0);

    for (int k = 0; k < pasos_totales; k++) {
        int i = distribution(generator);
        int j = distribution(generator);

        double delta_E = 0.0;
        if (S[i][j] == 1) {
            for (int x = 0; x < N; x++) {
                for (int y = 0; y < N; y++) {
                    delta_E += w[i][j][x][y] * (S[x][y] - 0.5);
                }
            }
        } else {
            for (int x = 0; x < N; x++) {
                for (int y = 0; y < N; y++) {
                    delta_E += w[i][j][x][y] * (0.5 - S[x][y]);
                }
            }
        }

        double p = min(1.0, exp(-delta_E / T));
        double r = r_distribution(generator);

        if (r < p) {
            S[i][j] = 1 - S[i][j]; 
        }
    }
}

// *********************Función para ver cuantos patrones se recuerdan en función de la cantidad de patrones almacenados*****************

void patrones_recordados() {
    int N = 20; 
    double T = 1e-4;
    int max_P = 100; 
    
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937_64 generator(seed);
    uniform_int_distribution<int> dist(0, 1);

    ofstream fich_capacidad("memoria_3.txt");
    fich_capacidad << "P_almacenados P_recordados\n";

    for (int P = 1; P <= max_P; P++) {
        // Generar P patrones aleatorios
        vector<vector<vector<int>>> patrones(P, vector<vector<int>>(N, vector<int>(N)));
        for (int mu = 0; mu < P; mu++) {
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    patrones[mu][i][j] = dist(generator);
                }
            }
        }

        // Calcular matriz de pesos 'w'
        vector<vector<vector<vector<double>>>> w = itereaccion_particulas(P, patrones, N);

        // Calcular a
        vector<double> a(P, 0.0);
        for(int mu = 0; mu < P; mu++){
            double suma_a = 0;
            for(int i=0; i<N; i++) {
                for(int j=0; j<N; j++) {
                    suma_a += patrones[mu][i][j];
                }
            }
            a[mu] = suma_a / (N * N);
        }

        int patrones_recordados = 0;

        // Probar si recuerda cada uno de los P patrones
        for (int mu = 0; mu < P; mu++) {
            vector<vector<int>> S = patrones[mu]; 
            
            metropolis_rapido(T, generator, w, S, N); 
            
            double solapamiento_final = calcular_solapamiento(S, patrones[mu], a[mu], N);
            
            if (solapamiento_final > 0.75) {
                patrones_recordados++;
            }
        }

        fich_capacidad << P << " " << patrones_recordados << "\n";
        cout << " P = " << P << " | Recordados = " << patrones_recordados << endl;
    }
    fich_capacidad.close();
}

int main() {
    patrones_recordados();
    return 0;
}