#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <set>
#include <math.h>

using namespace std;

// Hash retorna un numero de 64 o 32 bits dependiendo de la arquitectura de la maquina

int b = 14;
int M = pow(2,b);
hash<string> str_hash;
vector<int> sketch1(M, 0);
vector<int> sketch2(M, 0);

// Calculamos los parametros
double calculatedam(){
    double am = 0; 
    if(M == 16)
        am = 0.673;
    else if( M == 32)
        am = 0.697;
    else if( M == 64)
        am = 0.709;
    else
        am = 0.7213/(1.0+(1.079/(double)M));
    return am;
} 

// Union de sketchs manteniendo el maximo de ambos en el sketch1
void unionHyperLogLog()
{
    for (size_t i = 0; i < M; i++)
    {
        sketch1[i] = max(sketch1[i], sketch2[i]);
    }
}

//
unsigned long bitPosExtracted(unsigned long number)
{
    return (number >> (64 - b));
}

unsigned long bitExtracted(unsigned long number)
{
    return (number << b);
}

double estimateCardinality(vector<int> &sketch)
{
    double estimate = 0;
    double sumatoria1 = 0;
    double sumatoria2 = 0;
    for (int t = 0; t < M; t++)
    {
        double aux = (double)sketch[t];
        sumatoria2 = pow(2.0, -aux);

        sumatoria1 = sumatoria1 + sumatoria2;
    }
    estimate = 1.0 / sumatoria1;
    double am = calculatedam();
    estimate = estimate * am * pow((double)M, 2);
    double Erange = estimate;

    /************************************************************/

    double small_range = (5.0 / 2.0) * M;
    if (estimate <= small_range)
    {
        int numOfZeros = 0;
        for (int i = 0; i < M; i++)
        {
            if (sketch[i] == 0)
            {
                numOfZeros++;
            }
        }
        if (numOfZeros != 0)
        {
            Erange = M * (double)log2(M / numOfZeros);
        }
        else
        {
            Erange = estimate;
        }
    }

    return Erange;
}

void hyperLogLog(vector<int> &sketch, string kmer){
    unsigned long long hasheado;
    long pos = 0;
    int countZeros=0;
    hasheado = str_hash(kmer);
    pos = bitPosExtracted(hasheado);
    countZeros = __builtin_clz(bitExtracted(hasheado));
    if (sketch[pos] < countZeros + 1)
    {
        sketch[pos] = countZeros + 1;
    }
}

double jaccardHLL(double A, double B){
    unionHyperLogLog();
    double estimateUnion = estimateCardinality(sketch1);
    return (A + B - estimateUnion)/estimateUnion;
}

double cartesianHLL(double A, double B){
    return A*B;
}

double differenceHLL(double A, double B){
    unionHyperLogLog();
    double estimateUnion = estimateCardinality(sketch1);
    return A - (A + B - estimateUnion);
}

int main(int argc, char const *argv[])
{
    ifstream genomaA("Genoma1.txt");
    ifstream genomaB("Genoma2.txt");
    set<string> Real;
    string line;
    unsigned long long hasheado;
    long pos;
    int countZeros;
    double estimadoA = 0;
    double estimadoB = 0;
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
                    hyperLogLog(sketch1,kmer);
                }
            }
        }
        estimadoA = estimateCardinality(sketch1);
        cout << "Estimacion: " << estimadoA << endl;
    }
    else
    {
        cerr << "No se puso abrir el archivo '" << endl;
        return EXIT_FAILURE;
    }

    //Calculando cardinalidad para segundo genoma 
    if (genomaB.is_open())
    {
        cout << "El archivo se ha abierto correctamente" << endl;
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
                    hyperLogLog(sketch2,kmer);
                }
            }
        }
        estimadoB = estimateCardinality(sketch2);
        cout << "Estimacion: " << estimadoB << endl;
    }
    else
    {
        cerr << "No se puso abrir el archivo '" << endl;
        return EXIT_FAILURE;
    }

    //cout<<"El Jaccard de los genomas es :"<<jaccardHLL(estimadoA,estimadoB)<<endl;
    //cout<<"La cardinalidad de la diferencia del genoma A con el gneoma B es : " << differenceHLL( estimadoA, estimadoB)<<endl;
    //cout<<"La cardinalidad del producto cartesiano del genoma A y B es : " << cartesianHLL(estimadoA, estimadoB)<<endl;
}

