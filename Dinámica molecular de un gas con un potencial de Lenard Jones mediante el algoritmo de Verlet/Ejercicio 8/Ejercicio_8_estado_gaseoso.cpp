//SE HA DE DILUIR LAS PARTÍCULAS, POR ELLO, SE AUMENTA EL TAMAÑO DE LA CAJA PARA REDUCIR LA DENSIDAD
// SE HA DE AUMENTAR LA TEMPERATURA, POR ELLO, SE AUMENTA EL MÓDULO DE LA VELOCIDAD INICIAL DE CADA 
//PARTÍCULA

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <chrono>
#include <random>
#include <algorithm>
#include <string>

using namespace std;


void aceleracion(const vector < double > &x, const vector < double > &y, vector < double > &a_x, vector < double > &a_y, double L, int N);
void algoritmo_verlet_modificado(vector < double > &x, vector < double > &y, vector < double > &v_x, vector < double > &v_y, double L, int N, int num_bins, double r_max, double r_min, double ancho_bin);
void histograma(const vector < double > &x, const vector < double > &y,double L, int N ,double r_max, double r_min, double ancho_bin, vector < int > &histograma);


int main()
{


    // Semilla para el generador de números aleatorios basada en el tiempo actual
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();  
    mt19937_64 generator(seed);  // Inicializo el generador de números aleatorios

    double L = 22.0; //Tamaño del lado de la caja cuadrada
    int N = 144; //Número de átomos dentro de la caja

  // ***************************************Defino el histograma ****************************

    double r_max=3.0;
    double r_min = 1.0;
    int num_bins = 40;
    double ancho_bin = (r_max - r_min) / num_bins;
    
    //**********************Posiciones iniciales de los átomos********************************/
    
   // Red cuadrada equiespaciada
    vector < double > x(N), y(N); 


    int i = 0;
    for (int fila = 0; fila < sqrt(N); fila++)
    {
        for (int columna = 0; columna < sqrt(N); columna++)
        {
            x[i] = columna;
            y[i] = fila;
            
            i++;
        }
    }

    //**********************Velocidades iniciales de los átomos********************************/

   // Genero un array de ángulos aleatorios entre 0 y 2pi para cada átomo
    vector < double > angulos(N);

    double pi = acos(-1.0); // Defino el valor de pi

    uniform_real_distribution<double> distribucion_angular(0.0, 2.0 * pi);

    for (int j = 0; j < N; j++) 
    {
        angulos[j] = distribucion_angular(generator);
    }


    // El módulo de la velocidad de los átomos es uno

    vector < double > v_x(N), v_y(N);

    for (int k = 0; k < N; k++) 
    {
        v_x[k] = 4 * cos(angulos[k]);
        v_y[k] =4* sin(angulos[k]);
    }

    //Llamo al algoritmo de Verlet para evolucionar el sistema de átomos en el tiempo

    algoritmo_verlet_modificado(x, y, v_x, v_y, L, N, num_bins, r_max, r_min, ancho_bin);

    return 0;
}

// ********************Función para el Algoritmo de Verlet modificado *******************************

void algoritmo_verlet_modificado(vector < double > &x, vector < double > &y, vector < double > &v_x, vector < double > &v_y, double L, int N, int num_bins, double r_max, double r_min, double ancho_bin) 
{
    //**************************Paso 0. Calculo h y t, defino histograma*************************"""
    double t = 0.0;
    double h = 0.001;

    //Defino un vector para almacenar las particulas que hay en cada bin, uno para el estado líquido
    //otro para el estado sólido, y los inicializo a cero

    vector < int > h_gas(num_bins, 0); 

    //****************Paso 1. Calculo la aceleración inicial de cada átomo a partir de las posiciones iniciales**************

    vector < double > a_x(N, 0.0), a_y(N, 0.0);
    aceleracion(x, y, a_x, a_y, L, N);
   
    //Quiero evolucionar el sistema hasta t=2500, si el paso temporal es h=0.001,
    //el número de pasos temporales 

    ofstream archivo_velocidades("Velocidades_atomos_ej_8.txt");

   
    for (int i = 0; i < N; i++) 
    {
        archivo_velocidades << v_x[i] << ' ' << v_y[i] << ' ';
    }
    archivo_velocidades << '\n';


      int P= int(500 / h);

    // Inicio el ciclo temporal
    for (int k = 1; k <= P; k++) 
    {
       //************************paso 2. Calculo wi  y la posición en t+h********************** """

        vector < double > w_x(N), w_y(N);

        for (int i = 0; i < N; i++) 
        {
            w_x[i] = v_x[i] + (h / 2.0) * a_x[i];
            w_y[i] = v_y[i] + (h / 2.0) * a_y[i];
        }

        for (int i = 0; i < N; i++) 
        {
            x[i] = x[i] + h * w_x[i];
            y[i] = y[i] + h * w_y[i];
        }

        // Se ha de tener en cuenta las condiciones periódicas al calcular las nuevas posiciones de los átomos
        for (int i = 0; i < N; i++) {
            if (x[i] < 0.0) x[i] += L;
            if (x[i] >= L) x[i] -= L;
            if (y[i] < 0.0) y[i] += L;
            if (y[i] >= L) y[i] -= L;
        }
        
        //"""*********************Paso 3. Calculo la aceleración en t+h haciendo uso de las posiciones en t+h calculadas arriba justo*****"""

        for (int i = 0; i < N; i++) {
            a_x[i] = 0.0;
            a_y[i] = 0.0;
        }

        aceleracion(x, y, a_x, a_y, L, N);

        //****************************Paso 4. Calculo la velocidad en t+h**********************
    
        for (int i = 0; i < N; i++) 
        {
            v_x[i] = w_x[i] + (h / 2) * a_x[i];
            v_y[i] = w_y[i] + (h / 2) * a_y[i];
        }

        //********************Paso 5. Sumo al tiempo h******************************************

        t =t + h;

        // Guardo en el histograma de sólido desde t = 200 hasta t = 500 segundos, situación de equilibrio

        if ((k >= 100000) && (k % 100==0))
        {
            histograma(x, y, L, N, r_max, r_min, ancho_bin, h_gas);
        }

     // Guardo las posiciones y velocidades en  los ficheros

        for (int i = 0; i < N; i++) 
        {
        
            archivo_velocidades << v_x[i] << ' ' << v_y[i] << ' ';
        }
       
        archivo_velocidades << '\n';


    }
    archivo_velocidades.close();

    
    ofstream archivo_histograma_gas("Histograma_gas_ej_8.txt");

    
    for (int i = 0; i < num_bins; i++) 
    {
        archivo_histograma_gas << h_gas[i] << '\n';
    }


    archivo_histograma_gas.close();
 

}


//*************************************Función para almacener los átomos en los histogramas*******************/ 

void histograma(const vector < double > &x, const vector < double > &y,double L, int N ,double r_max, double r_min, double ancho_bin, vector < int > &histograma)
{

     // Hay que calcular la distancia entre un atomo i y el resto j,
     // se coje como referencia el átomo i=0, y se calcula la distancia 
     // entre el átomo i=0 y el resto de átomos j=1,2,...,N-1.

    int i = 0;
    for (int j = 1; j < N; j++)
    {
        // Calculo la distancia entre los átomos teniendo 
        // en cuenta las condiciones periódicas de la caja

        double x_ij = x[i] - x[j];
        if (x_ij > L / 2.0)   x_ij -= L;
        if (x_ij < -L / 2.0)  x_ij += L;

        double y_ij = y[i] - y[j];
        if (y_ij > L / 2.0)    y_ij -= L;
        if (y_ij < -L / 2.0)   y_ij += L;

        double r_ij = sqrt(x_ij * x_ij + y_ij * y_ij);

        // Ahora lo que hago es ver si la distancia entre el atomo i=0 y el resto 
        //esta entre 1 y 3, si es así, calculo en que bin del histograma se encuentra
        // y lo almaceno en el vector histograma

        if ((r_ij < r_max) && (r_ij >= r_min)) 
        {
            int bin = (int)((r_ij - r_min)/ancho_bin);

            histograma[bin]++;
        }
    }
}


// ********************Función para el calculo de la aceleración inicial de cada átomo********************


void aceleracion(const vector < double > &x, const vector < double > &y, vector < double > &a_x, vector < double > &a_y, double L, int N){
    for (int i = 0; i < N; i++) 
    {
        for (int j = i + 1; j < N; j++) // si  j = i + 1, no se calcula dos veces la interacción entre el átomo i y el j, por la ley de acción y reacción
        {
            // Condiciones periódicas de la caja: si la distancia entre los átomos i y j es mayor que L/2, 
            //entonces por condiciones periodicas, la distancia más corta entre los átomos i y j es L - r_ij, 
            //donde r_ij es la distancia entre los átomos i y j sin tener en cuenta las condiciones periódicas de la caja
            
            double x_ij = x[i] - x[j];
            if (x_ij > L / 2.0)     x_ij -= L;
            if (x_ij < -L / 2.0)    x_ij += L;

            double y_ij = y[i] - y[j];
            if (y_ij > L / 2.0)    y_ij -= L;
            if (y_ij < -L / 2.0)   y_ij += L;

            //Distancia entre los átomos i y j, teniendo en cuenta las condiciones periódicas de la caja.
            double r_ij = sqrt(x_ij * x_ij + y_ij * y_ij);

            // Si la distancia es mayor que 3.0, la iteracción entre los atmos es nula 
            if ((r_ij > 0.0) && (r_ij < 3.0)) 
            {
                double inv_r_ij = 1.0 / r_ij;

                double term = 24.0 * (2.0 * pow(inv_r_ij, 14) - pow(inv_r_ij, 8));

                a_x[i] += term * x_ij; 
                a_y[i] += term * y_ij;

                // Por la tercera ley de Newton, la fuerza que el átomo j ejerce sobre el átomo i es igual y opuesta a la fuerza que el átomo i ejerce sobre el átomo j
                a_x[j] -= term * x_ij;
                a_y[j] -= term * y_ij;
            }
        }
    }
}

