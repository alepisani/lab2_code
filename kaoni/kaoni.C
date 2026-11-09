#include <iostream>
#include <cmath>

using namespace std;

#define BR_true 5.583e-2
#define BR_true_err 0.024e-2

void kaoni()
{
    int iN3 = 1;
    int iN1 = 1;

    double N3 = (double)iN3;
    double N1 = (double)iN1;
    double Ntot = N3 + N1;

    double BR = N3 / (N3 + N1);
    double err_BR = 1/pow(Ntot, 2) * sqrt(N3 * N1 * Ntot);

    double err_BR_rel = err_BR / BR;
    double Z = (BR - BR_true) / sqrt(pow(err_BR, 2) + pow(BR_true_err, 2));

    cout << "\nBR: " << BR << " +/- " << err_BR << " (" << err_BR_rel*100. << "%)" << endl;
    
    cout << "Confronto con valore PDG: Z = " << Z << endl;

    if(err_BR_rel < 0.15 && fabs(Z) < 2)
        cout << "\n================================\n\n     E COFFEE BREAK SIGNORI\n\n================================\n\n";

}