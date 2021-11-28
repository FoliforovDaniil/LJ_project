#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "LJ_fluid_EOS.h"
#include "LJ_solid_EOS.h"

// calculations is an array which contains Common Thermodynamic Properties (reduced, dimensionless) in following order: 
// Density[0], Temperature[1], Phase[2](0 - VL, 1 - Fluid, 2 - LS, 3 - S(FCC), 4 - S(HCP){L -liquid, S - solid, V - vapour}), Pressure[3], Internal Energy[4], 
// Enthalpy[5], Isochoric Heat Capacity[6], Isobaric  Heat Capacity[7], dp/d(rho)  at T=const[8], dp/dT  at rho=const[9] 
// Helmholtz Free Energy[10] or its Derivatives[11-16] (A00, A10, A01, A20, A11, A02)
// The mixed partial derivative Axy=d^(x+y)(F/(NkBT))/d(1/T)^x/d(rho)^y * (1/T)^x * (rho)^y,
// where F is the Helmholtz free energy [J] (extensive), N is the number of particles,
// kB is the Boltzmann constant [J/K], T is the temperature [K], and rho is the density [mol/m3].

using namespace std;

//S_EOS(string type_of_lattice, double T, double Rho);
//L_G_EOS(double Rho, double T);

int argmin(vector<float>& data, float rho)
{

    float min_d = abs(data[0] - rho);
    int   index = 0;

    for (int i = 0; i < data.size(); i++)
    {
        if (abs(data[i] - rho) < min_d)
        {
            min_d = abs(data[i] - rho);
            index = i;
        }
    }

    return index;
}

bool cross(float* a, float* b, float* p)
{
    float f_k = (a[1] - b[1]) / (a[0] - b[0]);
    float f_b = a[1] - f_k * a[0];
    //cout << f_k << endl;
    //cout << f_b << endl;

    if (p[1] > (f_k * p[0] + f_b)) return true;
    else                           return false;
}

string check_value
(
    float rho, float T,
    vector<float>& data_ro_gas, vector<float>& data_T, vector<float>& data_ro_liq,
    vector<float>& data1_ro_liq, vector<float>& data1_T, vector<float>& data1_ro_solid,
    vector<float>& data2_ro_HCP, vector<float>& data2_T, vector<float>& data2_ro_FCC,
    vector<float>& data3_ro_gas, vector<float>& data3_T, vector<float>& data3_ro_solid
)
{
    float p[2] = { rho, T };
    int  idx, idx1, idx2, idx3;

    if (rho < 0.001944)
    {
        idx = argmin(data3_ro_gas, rho);

        if (data3_ro_gas[idx] - rho <= 0) idx1 = idx + 1;
        else                              idx1 = idx - 1;


        if (idx >= data3_ro_gas.size() - 1)
        {
            idx = data3_ro_gas.size() - 1;
            idx1 = idx - 1;
        }
        if (idx == 0) idx1 = 1;


        float a[2] = { data3_ro_gas[idx], data3_T[idx] };
        float b[2] = { data3_ro_gas[idx1], data3_T[idx1] };


        if (cross(a, b, p))  return string("L");
        else                 return string("NOT USING");
    }
    if (rho <= 0.31)
    {
        idx = argmin(data_ro_gas, rho);

        if (data_ro_gas[idx] - rho <= 0) idx1 = idx + 1;
        else                             idx1 = idx - 1;

        if (idx >= data_ro_gas.size() - 1)
        {
            idx = data_ro_gas.size() - 1;
            idx1 = idx - 1;
        }
        if (idx == 0) idx1 = 1;

        float a[2] = { data_ro_gas[idx], data_T[idx] };
        float b[2] = { data_ro_gas[idx1], data_T[idx1] };

        if (cross(a, b, p))  return string("L");
        else                 return string("NOT USING");
    }
    if (rho > 0.31)
    {
        if (rho <= 0.845999)
        {
            idx = argmin(data_ro_liq, rho);

            if (data_ro_liq[idx] - rho <= 0) idx1 = idx + 1;
            else                             idx1 = idx - 1;

            if (idx >= data_ro_liq.size() - 1)
            {
                idx = data_ro_liq.size() - 1;
                idx1 = idx - 1;
            }
            if (idx == 0) idx1 = 1;

            float  a[2] = { data_ro_liq[idx], data_T[idx] };
            float  b[2] = { data_ro_liq[idx1], data_T[idx1] };

            if (cross(a, b, p))  return string("L");
            else                 return string("NOT USING");
        }

        if (rho >= 0.85)
        {

            idx = argmin(data1_ro_liq, rho);

            if (data1_ro_liq[idx] - rho <= 0) idx1 = idx + 1;
            else                             idx1 = idx - 1;

            if (idx >= data1_ro_liq.size() - 1)
            {
                idx = data1_ro_liq.size() - 1;
                idx1 = idx - 1;
            }
            if (idx == 0) idx1 = 1;

            float a[2] = { data1_ro_liq[idx], data1_T[idx] };
            float b[2] = { data1_ro_liq[idx1], data1_T[idx1] };

            if (cross(a, b, p))  return string("L");
            else
            {
                idx2 = argmin(data1_ro_solid, rho);

                if (data1_ro_solid[idx2] - rho <= 0) idx3 = idx2 + 1;
                else                             idx3 = idx2 - 1;

                if (idx2 >= data1_ro_solid.size() - 1)
                {
                    idx2 = data1_ro_solid.size() - 1;
                    idx3 = idx - 1;
                }
                if (idx2 == 0) idx3 = 1;

                float a1[2] = { data1_ro_solid[idx2], data1_T[idx2] };
                float b1[2] = { data1_ro_solid[idx3], data1_T[idx3] };

                if (cross(a1, b1, p))  return string("LS");
                else
                {
                    if (rho >= 0.959398 and rho <= 1.247067)
                    {
                        idx = argmin(data3_ro_solid, rho);


                        if (data3_ro_solid[idx] - rho <= 0) idx1 = idx + 1;
                        else                                idx1 = idx - 1;

                        if (idx >= data3_ro_solid.size() - 1)
                        {
                            idx = data3_ro_solid.size() - 1;
                            idx1 = idx - 1;
                        }
                        if (idx == 0) idx1 = 1;

                        float a[2] = { data3_ro_solid[idx], data3_T[idx] };
                        float b[2] = { data3_ro_solid[idx1], data3_T[idx1] };

                        if (cross(a, b, p))
                        {
                            if (T < 0.14) return string("HCP");
                            else          return string("FCC");
                        }
                        else return string("NOT USING");
                    }

                    if (rho >= 1.25)
                    {
                        idx = argmin(data2_ro_HCP, rho);

                        if (data2_ro_HCP[idx] - rho <= 0) idx1 = idx + 1;
                        else                              idx1 = idx - 1;

                        if (idx >= data2_ro_HCP.size() - 1)
                        {
                            idx = data2_ro_HCP.size() - 1;
                            idx1 = idx - 1;
                        }
                        if (idx == 0) idx1 = 1;

                        float a[2] = { data2_ro_HCP[idx], data2_T[idx] };
                        float b[2] = { data2_ro_HCP[idx1], data2_T[idx1] };

                        if (cross(a, b, p))  return string("FCC");
                        else                 return string("HCP");
                    }
                }
            }
        }
    }
}


vector<double> function(float T, float Rho)
{
    FILE* file = fopen("z:\\data_dan\\Binodal_LJ_EOS_Vrabec2016.dat", "r");
    FILE* file1 = fopen("z:\\data_dan\\LJSolidLiqSadus2010.dat", "r");
    FILE* file2 = fopen("z:\\data_dan\\LJ_HCP_FCC_AdidharmaTan2016.dat", "r");
    FILE* file3 = fopen("z:\\data_dan\\Sublimation_LJAdidharmaTan2016.dat", "r");

    vector<float> 	data_ro_gas, data_T, data_ro_liq;
    vector<float> 	data1_ro_liq, data1_T, data1_ro_solid;
    vector<float> 	data2_ro_HCP, data2_T, data2_ro_FCC;
    vector<float>   data3_ro_gas, data3_T, data3_ro_solid;

    float rho, t, P, muT, ro_gas, ro_solid, dmu, ro_liq, ro_HCP, ro_FCC;
    char  solid_type[10];

    while (fscanf(file, "%f%f%f%f", &t, &P, &ro_gas, &ro_liq))
    {
        data_T.push_back(t);
        data_ro_gas.push_back(ro_gas);
        data_ro_liq.push_back(ro_liq);
    }
    while (fscanf(file1, "%f%f%f%f", &t, &P, &ro_liq, &ro_solid))
    {
        data1_ro_liq.push_back(ro_liq);
        data1_T.push_back(t);
        data1_ro_solid.push_back(ro_solid);
    }
    while (fscanf(file2, "%f%f%f%f%f%f", &t, &P, &muT, &ro_HCP, &ro_FCC, &dmu))
    {
        data2_ro_HCP.push_back(ro_HCP);
        data2_T.push_back(t);
        data2_ro_FCC.push_back(ro_FCC);
    }
    while (fscanf(file3, "%f%f%f%f%f%f%s", &t, &P, &muT, &ro_gas, &ro_solid, &dmu, solid_type))
    {
        data3_ro_gas.push_back(ro_gas);
        data3_T.push_back(t);
        data3_ro_solid.push_back(ro_solid);
    }

    string tmp;
    vector<double> calculations(20, 0);
    
    tmp = check_value
    (
        Rho, T,
        data_ro_gas, data_T, data_ro_liq,
        data1_ro_liq, data1_T, data1_ro_solid,
        data2_ro_HCP, data2_T, data2_ro_FCC,
        data3_ro_gas, data3_T, data3_ro_solid
    );
    cout << tmp << endl;
    //cout << T << endl;
    if (tmp != "NOT USING")
    {
        if (tmp == "L")
        {
            //cout << endl << "Current point is located in the fluid phase." << endl;
            calculations = L_G_EOS(Rho, T);
        }
        if (tmp == "FCC")
        {
            //cout << endl << "Current point is located in the solid phase with FCC lattice type." << endl;
            calculations = S_EOS("FCC", T, Rho);
        }
        if (tmp == "HCP")
        {
            //cout << endl << "Current point is located in the solid phase with HCP lattice type." << endl;
            calculations = S_EOS("HCP", T, Rho);
        }
        if (tmp == "LS")
        {
            calculations[0] = Rho;
            calculations[1] = T;
            calculations[2] = 2;
        }
    }
    else
    {
        if (Rho <= 0.845999 and T >= 0.692)
        {
            calculations[0] = Rho;
            calculations[1] = T;
            calculations[2] = 0;
        }
        else
        {
            calculations[0] = Rho;
            calculations[1] = T;
            calculations[2] = 2;
        }
    }

    fclose(file);
    fclose(file1);
    fclose(file2);
    fclose(file3);

    return calculations;
}

int main()
{
    for (float rho = 0; rho < 1.5; rho+= 0.25) 
    {
        for (float T = 0; T < 3; T += 0.25)
        {
            vector<double> a = function(T, rho);
            cout << " " << a[2] << " (" << rho<< ", " << T << ")" << endl;
        }
        //cout << endl;
    }
    return 0;
}