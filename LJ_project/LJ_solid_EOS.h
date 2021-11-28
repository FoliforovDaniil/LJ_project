#pragma once

#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
using namespace std;

double polysumm(string required_value, string type_of_lattice, double T, double Rho, double** BasisFunPar_FCC, double** BasisFunPar_HCP, int NTerms, int NParameters)
{
    double result = 0;
    if (required_value == "U")
    {
        if (type_of_lattice == "FCC")
        {
            for (int i = 0; i < NTerms; i++)
            {
                for (int j = 0; j < NParameters; j++)
                {
                    result += (j + 1.0) * BasisFunPar_FCC[i][j] * pow(T, j + 1) * pow(Rho, i);
                }
            }

        }

        if (type_of_lattice == "HCP")
        {
            for (int i = 0; i < NTerms; i++)
            {
                for (int j = 0; j < NParameters; j++)
                {
                    result += (j + 1.0) * BasisFunPar_HCP[i][j] * pow(T, j + 1) * pow(Rho, i);
                }
            }
        }
    }
    if (required_value == "Z")
    {
        if (type_of_lattice == "FCC")
        {
            for (int i = 0; i < NTerms; i++)
            {
                for (int j = 0; j < NParameters; j++)
                {
                    result += i * BasisFunPar_FCC[i][j] * pow(T, j + 1) * pow(Rho, i);
                }
            }
        }
        if (type_of_lattice == "HCP")
        {
            for (int i = 0; i < NTerms; i++)
            {
                for (int j = 0; j < NParameters; j++)
                {
                    result += i * BasisFunPar_HCP[i][j] * pow(T, j + 1) * pow(Rho, i);
                }
            }
        }
    }
    if (required_value == "A")
    {
        if (type_of_lattice == "FCC")
        {
            for (int i = 0; i < NTerms; i++)
            {
                for (int j = 0; j < NParameters; j++)
                {
                    result += BasisFunPar_FCC[i][j] * pow(T, j + 1) * pow(Rho, i);
                }
            }

        }

        if (type_of_lattice == "HCP")
        {
            for (int i = 0; i < NTerms; i++)
            {
                for (int j = 0; j < NParameters; j++)
                {
                    result += BasisFunPar_HCP[i][j] * pow(T, j + 1) * pow(Rho, i);
                }
            }
        }
    }
    return result;
}

double constant_of_integration(string required_value, string type_of_lattice, double Rho, double* b_FCC, double* b_HCP, double c0_FCC, double c0_HCP, int NTerms)
{
    double result = 0;

    if (required_value == "full")
    {
        if (type_of_lattice == "FCC")
        {
            result += c0_FCC;
            for (int i = 0; i < NTerms; i++)
            {
                result += b_FCC[i] * pow(Rho, i + 1);
            }
        }
        if (type_of_lattice == "HCP")
        {
            result += c0_HCP;
            for (int i = 0; i < NTerms; i++)
            {
                result += b_HCP[i] * pow(Rho, i + 1);
            }
        }
    }
    if (required_value == "derivative")
    {
        if (type_of_lattice == "FCC")
        {
            for (int i = 0; i < NTerms; i++)
            {
                result += (i + 1.0) * b_FCC[i] * pow(Rho, i);
            }
        }
        if (type_of_lattice == "HCP")
        {
            for (int i = 0; i < NTerms; i++)
            {
                result += (i + 1.0) * b_HCP[i] * pow(Rho, i);
            }
        }
    }
    return result;
}

vector<double> S_EOS(string type_of_lattice, double T, double Rho)
{
    //Matrices  that contains the parameters of the correlation
    //Matrices BasisFunPar
    double** BasisFunPar_FCC;
    double** BasisFunPar_HCP;
    int NTerms = 5, NParameters = 4;

    BasisFunPar_FCC = new double* [NTerms];
    for (int i = 0; i < NTerms; i++) BasisFunPar_FCC[i] = new double[NParameters];

    BasisFunPar_HCP = new double* [NTerms];
    for (int i = 0; i < NTerms; i++) BasisFunPar_HCP[i] = new double[NParameters];

    BasisFunPar_FCC[0][0] = 36.33808682;
    BasisFunPar_FCC[0][1] = -41.137808;
    BasisFunPar_FCC[0][2] = 25.74698103;
    BasisFunPar_FCC[0][3] = -7.807991646;
    BasisFunPar_FCC[1][0] = -107.1699988;
    BasisFunPar_FCC[1][1] = 128.050431;
    BasisFunPar_FCC[1][2] = -82.46948719;
    BasisFunPar_FCC[1][3] = 25.42569589;
    BasisFunPar_FCC[2][0] = 119.9321063;
    BasisFunPar_FCC[2][1] = -149.8142253;
    BasisFunPar_FCC[2][2] = 98.98941549;
    BasisFunPar_FCC[2][3] = -30.98224644;
    BasisFunPar_FCC[3][0] = -60.0831842;
    BasisFunPar_FCC[3][1] = 77.97346577;
    BasisFunPar_FCC[3][2] = -52.73174376;
    BasisFunPar_FCC[3][3] = 16.73649778;
    BasisFunPar_FCC[4][0] = 11.33640974;
    BasisFunPar_FCC[4][1] = -15.21781104;
    BasisFunPar_FCC[4][2] = 10.51355501;
    BasisFunPar_FCC[4][3] = -3.381006354;

    BasisFunPar_HCP[0][0] = 36.334506;
    BasisFunPar_HCP[0][1] = -41.129105;
    BasisFunPar_HCP[0][2] = 25.74887467;
    BasisFunPar_HCP[0][3] = -7.80991675;
    BasisFunPar_HCP[1][0] = -107.175742;
    BasisFunPar_HCP[1][1] = 128.051854;
    BasisFunPar_HCP[1][2] = -82.46970233;
    BasisFunPar_HCP[1][3] = 25.42605775;
    BasisFunPar_HCP[2][0] = 119.927892;
    BasisFunPar_HCP[2][1] = -149.8169985;
    BasisFunPar_HCP[2][2] = 98.98717633;
    BasisFunPar_HCP[2][3] = -30.98183575;
    BasisFunPar_HCP[3][0] = -60.082307;
    BasisFunPar_HCP[3][1] = 77.9725065;
    BasisFunPar_HCP[3][2] = -52.73275033;
    BasisFunPar_HCP[3][3] = 16.73701325;
    BasisFunPar_HCP[4][0] = 11.340422;
    BasisFunPar_HCP[4][1] = -15.217111;
    BasisFunPar_HCP[4][2] = 10.513569;
    BasisFunPar_HCP[4][3] = -3.38070975;

    //Matrix for parametr b
    double b_FCC[5];
    double b_HCP[5];

    b_FCC[0] = 107.235706;
    b_FCC[1] = -124.131009;
    b_FCC[2] = 77.69834533;
    b_FCC[3] = -24.920071;
    b_FCC[4] = 3.2250032;

    b_HCP[0] = 109.259855;
    b_HCP[1] = -127.542396;
    b_HCP[2] = 80.566189;
    b_HCP[3] = -26.11326875;
    b_HCP[4] = 3.4201126;

    //Parametr c

    double c0_FCC = -33.049893;
    double c2_FCC = -14.4539210435;
    double c4_FCC = 6.06594009827;

    double c0_HCP = -33.52540495;
    double c2_HCP = -14.4548972779;
    double c4_HCP = 6.06614688455;

    double U, Z, A;
    vector<double> CTP(20, 0);
    //double* CTP = new double[20];
    //for (int i = 0; i < 20; i++) CTP[i] = 0;
    

    if (type_of_lattice == "FCC")
    {
        U = 1.5 + (c2_FCC * pow(Rho, 2) + c4_FCC * pow(Rho, 4)) / T - polysumm("U", type_of_lattice, T, Rho, BasisFunPar_FCC, BasisFunPar_HCP, NTerms, NParameters);
        Z = T * Rho * (1 + (2 * c2_FCC * pow(Rho, 2) + 4 * c4_FCC * pow(Rho, 4)) / T + polysumm("Z", type_of_lattice, T, Rho, BasisFunPar_FCC, BasisFunPar_HCP, NTerms, NParameters) + Rho * constant_of_integration("derivative", type_of_lattice, Rho, b_FCC, b_HCP, c0_FCC, c0_HCP, NTerms));
        A = -1.5 * log(T) + (c2_FCC * pow(Rho, 2) + c4_FCC * pow(Rho, 4)) / T + polysumm("A", type_of_lattice, T, Rho, BasisFunPar_FCC, BasisFunPar_HCP, NTerms, NParameters) + constant_of_integration("full", type_of_lattice, Rho, b_FCC, b_HCP, c0_FCC, c0_HCP, NTerms);
        CTP[0] = Rho;
        CTP[1] = T;
        CTP[2] = 3;
        CTP[3] = Z;
        CTP[4] = U;
        CTP[9] = A;
    }
    if (type_of_lattice == "HCP")
    {
        U = 1.5 + (c2_HCP * pow(Rho, 2) + c4_HCP * pow(Rho, 4)) / T - polysumm("U", type_of_lattice, T, Rho, BasisFunPar_FCC, BasisFunPar_HCP, NTerms, NParameters);
        Z = T * Rho * (1 + (2 * c2_HCP * pow(Rho, 2) + 4 * c4_HCP * pow(Rho, 4)) / T + polysumm("Z", type_of_lattice, T, Rho, BasisFunPar_FCC, BasisFunPar_HCP, NTerms, NParameters) + Rho * constant_of_integration("derivative", type_of_lattice, Rho, b_FCC, b_HCP, c0_FCC, c0_HCP, NTerms));
        A = -1.5 * log(T) + (c2_HCP * pow(Rho, 2) + c4_HCP * pow(Rho, 4)) / T + polysumm("A", type_of_lattice, T, Rho, BasisFunPar_FCC, BasisFunPar_HCP, NTerms, NParameters) + constant_of_integration("Full", type_of_lattice, Rho, b_FCC, b_HCP, c0_FCC, c0_HCP, NTerms);
        CTP[0] = Rho;
        CTP[1] = T;
        CTP[2] = 4;
        CTP[3] = Z;
        CTP[4] = U;
        CTP[9] = A;
    }
    
    /*
    cout << setprecision(15) << endl << "Common Thermodynamic Properties:" << endl;

    cout << endl << "Pressure (total):                   " << Z << " (reduced, dimensionless)";
    cout << endl << "Internal Energy (residual):         " << U << " (reduced, dimensionless)";
    cout << endl << "Helmholtz energy(residual):         " << A << " (reduced, dimensionless)";
    */

    for (int i = 0; i < NTerms; i++) delete[] BasisFunPar_FCC[i];
    for (int i = 0; i < NTerms; i++) delete[] BasisFunPar_HCP[i];
    delete[] BasisFunPar_FCC;
    delete[] BasisFunPar_HCP;

    return CTP;
}

