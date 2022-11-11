#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <set>
#include <math.h>
#define PHI 0.77351
//#include <bits/stdc++.h>

using namespace std;
hash<string> H1;
hash<long> H2;
int M = 16;                 // Error del 6.9%
vector<long> sketch1(M, 0); // Create Sketch & inicialized in 0
vector<long> sketch2(M, 0);
int r(long number)
{
    int count = 0;
    while ((number & 1) == 1)
    {
        count++;
        number = number >> 1;
    }
    return count;
}
long R(long x)
{
    return ~x & (x + 1);
}

void update(string Gen, vector<long> &sketch)
{
    unsigned long x = H1(Gen);
    unsigned long k = H2(x) % (M - 1);
    sketch[k] = sketch[k] | R(x);
}

double long estimacion(vector<long> sketch)
{
    int sum = 0;
    for (int i = 0; i < M; i++) // R = 2^r ----> Entonces r = log2(R);
        // sum += log2(R(sketch[i])); //]r(sketch[i]);    // r(sketch[i]);
        sum += r(sketch[i]);
    cout << "Suma" << sum << endl;
    double media = 1.0 * ((double)sum / (double)M);
    cout << "Media " << media << endl;
    double result = (M * pow(2.0, media)) / (double)PHI;
    return result;
}
// Hash retorna un numero de 64 o 32 bits dependiendo de la arquitectura de la maquina

// Si queremos que algoritmo de HyperLogLog tenga un error de 9,2% necesitamos de 128 buckets, es decir M = 128

int main(int argc, char const *argv[])
{

    ifstream genomaA("Genoma1.txt");
    ifstream genomaB("Genoma2.txt");
    set<string> Real;
    string line;
    string kmer;
    unsigned long hasheado;
    int pos;
    int countZeros;
    double estimado = 0;
    if (genomaA.is_open())
    {
        cout << "El archivo se ha abierto correctamente" << endl;

        int cantidad = 0;
        while (genomaA.peek() != EOF)
        {
            getline(genomaA, line);
            if (line[0] == '>')
            {
                continue;
            }
            int x = line.size();
            if (x >= 31)
            {
                for (int i = 0; i <= line.size() - 31; i++)
                {
                    string kmer = line.substr(i, 31);
                    update(kmer, sketch1);
                }
            }
            // cout << "Cantidad" << cantidad << endl;
        }
        cout << "Estimado Genoma A: " << estimacion(sketch1) << endl;
    }
    genomaA.close();
    if (genomaB.is_open())
    {
        cout << "El Genoma B se ha abierto correctamente" << endl;

        int cantidad = 0;
        while (genomaB.peek() != EOF)
        {
            getline(genomaB, line);
            if (line[0] == '>')
            {
                continue;
            }
            int x = line.size();
            if (x >= 31)
            {
                for (int i = 0; i <= line.size() - 31; i++)
                {
                    string kmer = line.substr(i, 31);
                    update(kmer, sketch2);
                }
            }
        }
        cout << "Estimado Genoma B: " << estimacion(sketch2) << endl;
    }
    genomaB.close();
    return 0;
}
/*#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <set>
#include <math.h>
#define PHI 0.77351
//#include <bits/stdc++.h>

using namespace std;
hash<string> H1;
hash<long> H2;
int M = 32;               // Error del 6.9%
vector<long> sketch(M, 0); // Create Sketch & inicialized in 0

int r(long number)
{
    int count = 0;
    while ((number & 1) == 1)
    {
        count++;
        number = number >> 1;
    }
    return count;
}
long R(long x)
{
    return ~x & (x + 1);
}

void update(string Gen)
{
    unsigned long x = H1(Gen);
    unsigned long k = H2(x) % (M - 1);
    sketch[k] = sketch[k] | R(x);
}

double long estimacion()
{
    int sum = 0;
    for (int i = 0; i < M; i++) // R = 2^r ----> Entonces r = log2(R);
        // sum += log2(R(sketch[i])); //]r(sketch[i]);    // r(sketch[i]);
        sum += r(sketch[i]);
    cout << "Suma" << sum << endl;
    double media = 1.0 * ((double)sum / (double)M);
    cout << "Media " << media << endl;
    double result = (M * pow(2.0, media)) / (double)PHI;
    return result;
}


int main(int argc, char const *argv[])
{

    ifstream genomaA("Genoma1.txt");
    ifstream genomaB("Genoma2.txt");
    set<string> Real;
    string line;
    string kmer;
    unsigned long hasheado;
    int pos;
    int countZeros;

    double estimado = 0;

    if (genomaA.is_open())
    {
        cout << "El archivo se ha abierto correctamente" << endl;

        int cantidad = 0;
        while (genomaA.peek() != EOF)
        {
            getline(genomaA, line);
            if (line[0] == '>')
            {
                continue;
            }
            int x = line.size();
            if (x >= 31)
            {
                for (int i = 0; i <= line.size() - 31; i++)
                {
                    string kmer = line.substr(i, 31);
                    Real.insert(kmer);
                    update(kmer);
                }
            }
            // cout << "Cantidad" << cantidad << endl;
        }
        cout << "Real :" << Real.size() << endl;
        cout << "Estimado :" << estimacion() << endl;
    }
    return 0;
}*/
