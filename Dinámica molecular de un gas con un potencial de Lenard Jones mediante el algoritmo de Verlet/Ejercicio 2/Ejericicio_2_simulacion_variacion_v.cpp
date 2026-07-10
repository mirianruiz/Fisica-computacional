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
void algoritmo_verlet(vector < double > &x, vector < double > &y, vector < double > &v_x, vector < double > &v_y, double L, int N);

int main()
{


    // Semilla para el generador de números aleatorios basada en el tiempo actual
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();  
    mt19937_64 generator(seed);  // Inicializo el generador de números aleatorios

    double L = 22.0; //Tamaño del lado de la caja cuadrada
    int N = 100; //Número de átomos dentro de la caja
    
    //**********************Posiciones iniciales de los átomos********************************/
    
    vector<double> x(N), y(N); 
    
    // Distribución uniforme para colocar los átomos dentro de la caja de tamaño L
    uniform_real_distribution<double> distribucion_pos(0.0, L);

    // Se ha de tener en cuenta que los átomos no pueden colocarse demasiado cerca unos de otros,
    //pues el potencial de Lennard-Jones es muy repulsivo a distancias cortas, y estos átomos
    //saldrían despedidos a velocidades enormes, lo cual provoca que el simulador ’explote’, por ello
    //se generan posiciones aleatorias de los átomos, pero se comprueba que no estén demasiado cerca
    // unos de otros, y si lo están, se genera otra posición aleatoria para ese átomo hasta que se encuentre una posición buena.
    //Se escoge una distancia mínima entre átomos de 1.1, a partir de dicha distancia, el potencial 
    // es altamente repulsivo.


    double distancia_minima = 1.1; 

    for (int i = 0; i < N; i++) 
    {
        while (true) // Bucle infinito que solo se rompe si encontramos una posición buena
        {
            // 1. Generamos una posición aleatoria candidata
            double pos_x_candidata = distribucion_pos(generator);
            double pos_y_candidata = distribucion_pos(generator);

            int j = 0;
            
            // 2. Comprobamos la distancia con todos los átomos ya colocados (de 0 a i-1)
            while (j < i) 
            {
                double dx = pos_x_candidata - x[j];
                double dy = pos_y_candidata - y[j];

                // Aplicamos condiciones de contorno periódicas a la distancia
                if (dx > L / 2.0) dx -= L;
                if (dx < -L / 2.0) dx += L;
                if (dy > L / 2.0) dy -= L;
                if (dy < -L / 2.0) dy += L;

                double r_ij = sqrt(dx * dx + dy * dy);

                // Si está demasiado cerca, rompemos este bucle de comprobación inmediatamente
                if (r_ij < distancia_minima) 
                {
                    break; 
                }
                
                j++; // Si no hay choque, pasamos al siguiente átomo
            }

            // 3. Si j es igual a i, significa que el bucle while(j < i) terminó de 
            // revisar todos los átomos sin interrupciones
            if (j == i) 
            {
                x[i] = pos_x_candidata;
                y[i] = pos_y_candidata;
                break; // Rompemos el bucle while(true) porque ya hemos colocado el átomo i
            }
        }
    }
        
    

    //**********************Velocidades iniciales de los átomos********************************/

  
    vector < double > v_x(N), v_y(N);
    uniform_real_distribution<double> distribucion_vx(0.0, 1.0);

    for (int k = 0; k < N; k++) 
    {
        v_x[k] = distribucion_vx(generator); // vx aleatoria entre 0 y 1
        v_y[k] = 0.0;                        // vy nula
    }
   
    //Llamo al algoritmo de Verlet para evolucionar el sistema de átomos en el tiempo

    algoritmo_verlet(x, y, v_x, v_y, L, N);



    return 0;
}


// ********************Función para el Algoritmo de Verlet********************

void algoritmo_verlet(vector < double > &x, vector < double > &y, vector < double > &v_x, vector < double > &v_y, double L, int N) 
{
    //**************************Paso 0. Calculo h y t*************************"""
    double t = 0.0;
    double h = 0.002;

    //****************Paso 1. Calculo la aceleración inicial de cada átomo a partir de las posiciones iniciales**************

    vector < double > a_x(N, 0.0), a_y(N, 0.0);
    aceleracion(x, y, a_x, a_y, L, N);



    //Guardo las posiciones y velocidades de los átomos en ficheros, de forma que cada fila 
    // corresponda a un tiempo t, y de dos en dos columnas corresponda a un átomo
    //de forma que se escriba con x1 y1 x2 y2... y lo mismo con las velocidades.

    ofstream archivo_posiciones("Posiciones_atomos_vy_0_2.txt");
    ofstream archivo_velocidades("Velocidades_atomos_vy_0_2.txt");

   
    for (int i = 0; i < N; i++) 
    {
        archivo_posiciones << x[i] << ' ' << y[i] << ' ';
        archivo_velocidades << v_x[i] << ' ' << v_y[i] << ' ';
    }
    archivo_posiciones << '\n';
    archivo_velocidades << '\n';

   
    //Quiero evolucionar el sistema hasta t=250, si el paso temporal es h=0.002,
    //el número de pasos temporales 

      int P= int(250 / h);

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

        // Guardo las posiciones y velocidades en  los ficheros

        for (int i = 0; i < N; i++) 
        {
            archivo_posiciones << x[i] << ' ' << y[i] << ' ';
            archivo_velocidades << v_x[i] << ' ' << v_y[i] << ' ';
        }
        archivo_posiciones << '\n';
        archivo_velocidades << '\n';
    }

    archivo_posiciones.close();
    archivo_velocidades.close();
}

// ********************Función para el calculo de la aceleración inicial de cada átomo********************


void aceleracion(const vector < double > &x, const vector < double > &y, vector < double > &a_x, vector < double > &a_y, double L, int N) {
    for (int i = 0; i < N; i++) 
    {
        for (int j = i + 1; j < N; j++) // si  j = i + 1, no se calcula dos veces la interacción entre el átomo i y el j, por la ley de acción y reacción
        {
            // Condiciones periódicas de la caja: si la distancia entre los átomos i y j es mayor que L/2, 
            //entonces por condiciones periodicas, la distancia más corta entre los átomos i y j es L - r_ij, 
            //donde r_ij es la distancia entre los átomos i y j sin tener en cuenta las condiciones periódicas de la caja
            
            double x_ij = x[i] - x[j];
            if (x_ij > L / 2.0)  x_ij -= L;
            if (x_ij < -L / 2.0)  x_ij += L;

            double y_ij = y[i] - y[j];
            if (y_ij > L / 2.0) y_ij -= L;
            if (y_ij < -L / 2.0) y_ij += L;

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