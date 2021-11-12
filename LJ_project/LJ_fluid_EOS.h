#pragma once

#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

double PolyTerm(double tau, double delta, double n, double t, double d, int ordertau, int orderdelta)
{
    double product_tau = 1.0, product_delta = 1.0;

    for (int i = 0; i < ordertau; i++) product_tau = product_tau * (t - i);
    for (int i = 0; i < orderdelta; i++) product_delta = product_delta * (d - i);

    return n * product_tau * pow(tau, t - ordertau) * product_delta * pow(delta, d - orderdelta);
}

double ExpTerm(double tau, double delta, double n, double t, double d, double l, int ordertau, int orderdelta)
{
    double product_tau = 1.0, g = 1.0;

    for (int i = 0; i < ordertau; i++) product_tau = product_tau * (t - i);

    if (orderdelta == 0) return n * product_tau * pow(tau, t - ordertau) * exp(-g * pow(delta, l)) * pow(delta, d);
    if (orderdelta == 1) return n * product_tau * pow(tau, t - ordertau) * pow(delta, (d - 1.0)) * (d - g * l * pow(delta, l)) * exp(-g * pow(delta, l));
    if (orderdelta == 2) return n * product_tau * pow(tau, t - ordertau) * pow(delta, (d - 2.0)) * ((d - g * l * pow(delta, l)) * (d - 1.0 - g * l * pow(delta, l)) - pow(g, 2) * pow(l, 2) * pow(delta, l)) * exp(-g * pow(delta, l));
}

double GaussTerm(double tau, double delta, double n, double t, double d, double eta, double beta, double gamma, double epsilon, int ordertau, int orderdelta)
{
    double term00 = pow(tau, t) * pow(delta, d) * exp(-eta * pow((delta - epsilon), 2.0) - beta * pow((tau - gamma), 2.0));

    if (ordertau == 0 && orderdelta == 0) return n * term00;
    if (ordertau == 1 && orderdelta == 0) return n * term00 * (t / tau - 2.0 * beta * (tau - gamma));
    if (ordertau == 2 && orderdelta == 0) return n * term00 * (pow(t / tau - 2.0 * beta * (tau - gamma), 2.0) - t / tau / tau - 2.0 * beta);
    if (ordertau == 0 && orderdelta == 1) return n * term00 * (d / delta - 2.0 * eta * (delta - epsilon));
    if (ordertau == 0 && orderdelta == 2) return n * term00 * (pow(d / delta - 2.0 * eta * (delta - epsilon), 2.0) - d / delta / delta - 2.0 * eta);
    if (ordertau == 1 && orderdelta == 1) return n * term00 * (t / tau - 2.0 * beta * (tau - gamma)) * (d / delta - 2.0 * eta * (delta - epsilon));
}

double ReturnTermValue(int i, double tau, double delta, int ordertau, int orderdelta, double** BasisFunPar)
{
    if (fabs(BasisFunPar[i][0] - 1.0) < 0.01) return PolyTerm(tau, delta, BasisFunPar[i][1], BasisFunPar[i][2], BasisFunPar[i][3], ordertau, orderdelta);
    if (fabs(BasisFunPar[i][0] - 2.0) < 0.01) return ExpTerm(tau, delta, BasisFunPar[i][1], BasisFunPar[i][2], BasisFunPar[i][3], BasisFunPar[i][4], ordertau, orderdelta);
    if (fabs(BasisFunPar[i][0] - 3.0) < 0.01) return GaussTerm(tau, delta, BasisFunPar[i][1], BasisFunPar[i][2], BasisFunPar[i][3], BasisFunPar[i][5], BasisFunPar[i][6], BasisFunPar[i][7], BasisFunPar[i][8], ordertau, orderdelta);
}

int L_G_EOS(double Rho, double T)
{
    // Critical temperature and density:
    double Tc = 1.32, Rhoc = 0.31;

    //Matrix that contains the parameters of the correlation
    double** BasisFunPar;
    int NTerms = 23, NParameters = 9;

    BasisFunPar = new double* [NTerms];
    for (int i = 0; i < NTerms; i++) BasisFunPar[i] = new double[NParameters];

    // TermID
    BasisFunPar[0][0] = 1.0;
    BasisFunPar[1][0] = 1.0;
    BasisFunPar[2][0] = 1.0;
    BasisFunPar[3][0] = 1.0;
    BasisFunPar[4][0] = 1.0;
    BasisFunPar[5][0] = 1.0;
    BasisFunPar[6][0] = 2.0;
    BasisFunPar[7][0] = 2.0;
    BasisFunPar[8][0] = 2.0;
    BasisFunPar[9][0] = 2.0;
    BasisFunPar[10][0] = 2.0;
    BasisFunPar[11][0] = 2.0;
    BasisFunPar[12][0] = 3.0;
    BasisFunPar[13][0] = 3.0;
    BasisFunPar[14][0] = 3.0;
    BasisFunPar[15][0] = 3.0;
    BasisFunPar[16][0] = 3.0;
    BasisFunPar[17][0] = 3.0;
    BasisFunPar[18][0] = 3.0;
    BasisFunPar[19][0] = 3.0;
    BasisFunPar[20][0] = 3.0;
    BasisFunPar[21][0] = 3.0;
    BasisFunPar[22][0] = 3.0;
    // n
    BasisFunPar[0][1] = 0.005208073;
    BasisFunPar[1][1] = 2.186252000;
    BasisFunPar[2][1] = -2.161016000;
    BasisFunPar[3][1] = 1.452700000;
    BasisFunPar[4][1] = -2.041792000;
    BasisFunPar[5][1] = 0.186952860;
    BasisFunPar[6][1] = -0.090988445;
    BasisFunPar[7][1] = -0.497456100;
    BasisFunPar[8][1] = 0.109014310;
    BasisFunPar[9][1] = -0.800559220;
    BasisFunPar[10][1] = -0.568839000;
    BasisFunPar[11][1] = -0.620862500;
    BasisFunPar[12][1] = -1.466717700;
    BasisFunPar[13][1] = 1.891469000;
    BasisFunPar[14][1] = -0.138370100;
    BasisFunPar[15][1] = -0.386964500;
    BasisFunPar[16][1] = 0.126570200;
    BasisFunPar[17][1] = 0.605781000;
    BasisFunPar[18][1] = 1.179189000;
    BasisFunPar[19][1] = -0.477326790;
    BasisFunPar[20][1] = -9.921857500;
    BasisFunPar[21][1] = -0.574793200;
    BasisFunPar[22][1] = 0.003772923;
    // t
    BasisFunPar[0][2] = 1.000;
    BasisFunPar[1][2] = 0.320;
    BasisFunPar[2][2] = 0.505;
    BasisFunPar[3][2] = 0.672;
    BasisFunPar[4][2] = 0.843;
    BasisFunPar[5][2] = 0.898;
    BasisFunPar[6][2] = 1.294;
    BasisFunPar[7][2] = 2.590;
    BasisFunPar[8][2] = 1.786;
    BasisFunPar[9][2] = 2.770;
    BasisFunPar[10][2] = 1.786;
    BasisFunPar[11][2] = 1.205;
    BasisFunPar[12][2] = 2.830;
    BasisFunPar[13][2] = 2.548;
    BasisFunPar[14][2] = 4.650;
    BasisFunPar[15][2] = 1.385;
    BasisFunPar[16][2] = 1.460;
    BasisFunPar[17][2] = 1.351;
    BasisFunPar[18][2] = 0.660;
    BasisFunPar[19][2] = 1.496;
    BasisFunPar[20][2] = 1.830;
    BasisFunPar[21][2] = 1.616;
    BasisFunPar[22][2] = 4.970;
    // d
    BasisFunPar[0][3] = 4.0;
    BasisFunPar[1][3] = 1.0;
    BasisFunPar[2][3] = 1.0;
    BasisFunPar[3][3] = 2.0;
    BasisFunPar[4][3] = 2.0;
    BasisFunPar[5][3] = 3.0;
    BasisFunPar[6][3] = 5.0;
    BasisFunPar[7][3] = 2.0;
    BasisFunPar[8][3] = 2.0;
    BasisFunPar[9][3] = 3.0;
    BasisFunPar[10][3] = 1.0;
    BasisFunPar[11][3] = 1.0;
    BasisFunPar[12][3] = 1.0;
    BasisFunPar[13][3] = 1.0;
    BasisFunPar[14][3] = 2.0;
    BasisFunPar[15][3] = 3.0;
    BasisFunPar[16][3] = 3.0;
    BasisFunPar[17][3] = 2.0;
    BasisFunPar[18][3] = 1.0;
    BasisFunPar[19][3] = 2.0;
    BasisFunPar[20][3] = 3.0;
    BasisFunPar[21][3] = 1.0;
    BasisFunPar[22][3] = 1.0;
    // l
    BasisFunPar[0][4] = 0.0;
    BasisFunPar[1][4] = 0.0;
    BasisFunPar[2][4] = 0.0;
    BasisFunPar[3][4] = 0.0;
    BasisFunPar[4][4] = 0.0;
    BasisFunPar[5][4] = 0.0;
    BasisFunPar[6][4] = 1.0;
    BasisFunPar[7][4] = 2.0;
    BasisFunPar[8][4] = 1.0;
    BasisFunPar[9][4] = 2.0;
    BasisFunPar[10][4] = 2.0;
    BasisFunPar[11][4] = 1.0;
    BasisFunPar[12][4] = 0.0;
    BasisFunPar[13][4] = 0.0;
    BasisFunPar[14][4] = 0.0;
    BasisFunPar[15][4] = 0.0;
    BasisFunPar[16][4] = 0.0;
    BasisFunPar[17][4] = 0.0;
    BasisFunPar[18][4] = 0.0;
    BasisFunPar[19][4] = 0.0;
    BasisFunPar[20][4] = 0.0;
    BasisFunPar[21][4] = 0.0;
    BasisFunPar[22][4] = 0.0;
    // eta
    BasisFunPar[0][5] = 0;
    BasisFunPar[1][5] = 0;
    BasisFunPar[2][5] = 0;
    BasisFunPar[3][5] = 0;
    BasisFunPar[4][5] = 0;
    BasisFunPar[5][5] = 0;
    BasisFunPar[6][5] = 0;
    BasisFunPar[7][5] = 0;
    BasisFunPar[8][5] = 0;
    BasisFunPar[9][5] = 0;
    BasisFunPar[10][5] = 0;
    BasisFunPar[11][5] = 0;
    BasisFunPar[12][5] = 2.067;
    BasisFunPar[13][5] = 1.522;
    BasisFunPar[14][5] = 8.82;
    BasisFunPar[15][5] = 1.722;
    BasisFunPar[16][5] = 0.679;
    BasisFunPar[17][5] = 1.883;
    BasisFunPar[18][5] = 3.925;
    BasisFunPar[19][5] = 2.461;
    BasisFunPar[20][5] = 28.2;
    BasisFunPar[21][5] = 0.753;
    BasisFunPar[22][5] = 0.82;
    // beta
    BasisFunPar[0][6] = 0;
    BasisFunPar[1][6] = 0;
    BasisFunPar[2][6] = 0;
    BasisFunPar[3][6] = 0;
    BasisFunPar[4][6] = 0;
    BasisFunPar[5][6] = 0;
    BasisFunPar[6][6] = 0;
    BasisFunPar[7][6] = 0;
    BasisFunPar[8][6] = 0;
    BasisFunPar[9][6] = 0;
    BasisFunPar[10][6] = 0;
    BasisFunPar[11][6] = 0;
    BasisFunPar[12][6] = 0.625;
    BasisFunPar[13][6] = 0.638;
    BasisFunPar[14][6] = 3.91;
    BasisFunPar[15][6] = 0.156;
    BasisFunPar[16][6] = 0.157;
    BasisFunPar[17][6] = 0.153;
    BasisFunPar[18][6] = 1.16;
    BasisFunPar[19][6] = 1.73;
    BasisFunPar[20][6] = 383;
    BasisFunPar[21][6] = 0.112;
    BasisFunPar[22][6] = 0.119;
    // gamma
    BasisFunPar[0][7] = 0;
    BasisFunPar[1][7] = 0;
    BasisFunPar[2][7] = 0;
    BasisFunPar[3][7] = 0;
    BasisFunPar[4][7] = 0;
    BasisFunPar[5][7] = 0;
    BasisFunPar[6][7] = 0;
    BasisFunPar[7][7] = 0;
    BasisFunPar[8][7] = 0;
    BasisFunPar[9][7] = 0;
    BasisFunPar[10][7] = 0;
    BasisFunPar[11][7] = 0;
    BasisFunPar[12][7] = 0.71;
    BasisFunPar[13][7] = 0.86;
    BasisFunPar[14][7] = 1.94;
    BasisFunPar[15][7] = 1.48;
    BasisFunPar[16][7] = 1.49;
    BasisFunPar[17][7] = 1.945;
    BasisFunPar[18][7] = 3.02;
    BasisFunPar[19][7] = 1.11;
    BasisFunPar[20][7] = 1.17;
    BasisFunPar[21][7] = 1.33;
    BasisFunPar[22][7] = 0.24;
    // epsilon
    BasisFunPar[0][8] = 0;
    BasisFunPar[1][8] = 0;
    BasisFunPar[2][8] = 0;
    BasisFunPar[3][8] = 0;
    BasisFunPar[4][8] = 0;
    BasisFunPar[5][8] = 0;
    BasisFunPar[6][8] = 0;
    BasisFunPar[7][8] = 0;
    BasisFunPar[8][8] = 0;
    BasisFunPar[9][8] = 0;
    BasisFunPar[10][8] = 0;
    BasisFunPar[11][8] = 0;
    BasisFunPar[12][8] = 0.2053;
    BasisFunPar[13][8] = 0.409;
    BasisFunPar[14][8] = 0.6;
    BasisFunPar[15][8] = 1.203;
    BasisFunPar[16][8] = 1.829;
    BasisFunPar[17][8] = 1.397;
    BasisFunPar[18][8] = 1.39;
    BasisFunPar[19][8] = 0.539;
    BasisFunPar[20][8] = 0.934;
    BasisFunPar[21][8] = 2.369;
    BasisFunPar[22][8] = 2.43;

    double tau, delta;

    tau = Tc / T;
    delta = Rho / Rhoc;

    //residual and dimensionless Helmholtz free energy derivatives
    double A00res = 0.0, A10res = 0.0, A01res = 0.0, A20res = 0.0, A11res = 0.0, A02res = 0.0;

    for (int i = 0; i < NTerms; i++)
    {
        A00res += ReturnTermValue(i, tau, delta, 0, 0, BasisFunPar);
        A10res += ReturnTermValue(i, tau, delta, 1, 0, BasisFunPar) * tau;
        A01res += ReturnTermValue(i, tau, delta, 0, 1, BasisFunPar) * delta;
        A20res += ReturnTermValue(i, tau, delta, 2, 0, BasisFunPar) * tau * tau;
        A11res += ReturnTermValue(i, tau, delta, 1, 1, BasisFunPar) * tau * delta;
        A02res += ReturnTermValue(i, tau, delta, 0, 2, BasisFunPar) * delta * delta;
    }

    //common thermodynamic properties
    double p, ures, hres, cvres, cpres, dpdrho, dpdt;

    ures = A10res * T;
    hres = (A01res + A10res) * T;
    p = T * Rho * (1.0 + A01res);
    dpdrho = T * (1.0 + 2.0 * A01res + A02res);
    dpdt = Rho * (1.0 + A01res - A11res);
    cvres = -A20res;
    cpres = -A20res + (1.0 + A01res - A11res) * (1.0 + A01res - A11res) / (1.0 + 2.0 * A01res + A02res) - 1.0;

    cout << setprecision(15) << endl << "Common Thermodynamic Properties:" << endl;

    cout << endl << "Pressure (total):                   " << p << " (reduced, dimensionless)";
    cout << endl << "Internal Energy (residual):         " << ures << " (reduced, dimensionless)";
    cout << endl << "Enthalpy (residual):                " << hres << " (reduced, dimensionless)";
    cout << endl << "Isochoric Heat Capacity (residual): " << cvres << " (reduced, dimensionless)";
    cout << endl << "Isobaric  Heat Capacity (residual): " << cpres << " (reduced, dimensionless)";
    cout << endl << "dp/d(rho)  at T=const. (total):     " << dpdrho << " (reduced, dimensionless)";
    cout << endl << "dp/dT  at rho=const. (total):       " << dpdt << " (reduced, dimensionless)" << endl;

    cout << endl << "Helmholtz Free Energy Derivatives:" << endl;

    cout << endl << "A00 (residual): " << A00res << " (reduced, dimensionless)";
    cout << endl << "A10 (residual): " << A10res << " (reduced, dimensionless)";
    cout << endl << "A01 (residual): " << A01res << " (reduced, dimensionless)";
    cout << endl << "A20 (residual): " << A20res << " (reduced, dimensionless)";
    cout << endl << "A11 (residual): " << A11res << " (reduced, dimensionless)";
    cout << endl << "A02 (residual): " << A02res << " (reduced, dimensionless)" << endl;

    cout << endl << "The mixed partial derivative Axy=d^(x+y)(F/(NkBT))/d(1/T)^x/d(rho)^y * (1/T)^x * (rho)^y,";
    cout << " where F is the Helmholtz free energy [J] (extensive), N is the number of particles,";
    cout << " kB is the Boltzmann constant [J/K], T is the temperature [K], and rho is the density [mol/m3].";
    cout << endl << "Note that Axy is dimensionless, independent from the choice of units (Lennard-Jones or SI)." << endl;

    for (int i = 0; i < NTerms; i++) delete[] BasisFunPar[i];
    delete[] BasisFunPar;

    return 0;
}

