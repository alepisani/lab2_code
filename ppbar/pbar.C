#include <iostream>
#include <cmath>



using namespace std;

void pbar()
{
    int iNp = 176;
    int iNd = 123;
    int iN0 = 1063;


    double L1_mis = 50.27;              //cm
    double L1_err = 0.1;                //cm
    double L1_true = 28.1743;           //cm
    double L1_true_err = 0;             //cm

    double L2_mis = 38.73;              //cm
    double L2_err = 0.1;                //cm
    double L2_true = 28.1793;           //cm
    double L2_true_err = 0;             //cm

    double L_mis = 74.17;               //cm
    double L_mis_err = 0.1;             //cm

    double K = 2.37e-23;               //cm3

    double Np = (double)iNp;
    double Nd = (double)iNd;
    double N0 = (double)iN0;

    double err_Np = sqrt(Np);
    double err_Nd = sqrt(Nd);
    double err_N0 = sqrt(N0);

    double G1 = L1_mis / L1_true;
    double G1_err = G1 * sqrt(pow(L1_err / L1_mis, 2) + pow(L1_true_err / L1_true, 2));
    double G2 = L2_mis / L2_true;
    double G2_err = G2 * sqrt(pow(L2_err / L2_mis, 2) + pow(L2_true_err / L2_true, 2));

    double G = 0.5 * (G1 + G2);
    double G_err = 0.5 * sqrt(pow(G1_err, 2) + pow(G2_err, 2));

    double Ltrue = L_mis / G;
    double Ltrue_err = Ltrue * sqrt(pow(L_mis_err / L_mis, 2) + pow(G_err / G, 2));

    cout << "\n<G>: " << G << " +/- " << G_err << endl;
    cout << "L_true (cm): " << Ltrue << " +/- " << Ltrue_err << endl;


    double logtot = log(N0) - log(N0 - Np - Nd);
    double err_logtot = 1./(N0-Np-Nd) * sqrt(Np + Nd + pow(Np+Nd,2)/N0);

    double sigma_tot = (K/Ltrue) * logtot;
    double err_sigma_tot = sigma_tot * sqrt(pow(Ltrue_err / Ltrue, 2) + pow(err_logtot/logtot, 2));

    double BR_p = Np / (Np + Nd);
    double err_BR_p = 1/pow(Np + Nd, 2) * sqrt(Np * Nd * (Np + Nd));
    double sigma_p = sigma_tot * BR_p;
    double err_sigma_p = sigma_p * sqrt(pow(err_sigma_tot / sigma_tot, 2) + pow(err_BR_p / BR_p, 2));

    double BR_n = Nd / (Np + Nd);
    double err_BR_n = 1/pow(Np + Nd, 2) * sqrt(Np * Nd * (Np + Nd));
    double sigma_n = sigma_tot * BR_n;
    double err_sigma_n = sigma_n * sqrt(pow(err_sigma_tot / sigma_tot, 2) + pow(err_BR_n / BR_n, 2));

    cout << "\n==========================================================\n\n";
    cout << "Sigma_total (mb): " << sigma_tot*1e27 << " +/- " << err_sigma_tot*1e27 << "   (" << (err_sigma_tot / sigma_tot * 100.) << "%)" << endl;
    cout << "Sigma_p (mb): " << sigma_p*1e27 << " +/- " << err_sigma_p*1e27 << "   (" << (err_sigma_p / sigma_p * 100.) << "%)" << endl;
    cout << "Sigma_n (mb): " << sigma_n*1e27 << " +/- " << err_sigma_n*1e27 << "   (" << (err_sigma_n / sigma_n * 100.) << "%)" << endl;
    cout << "\n==========================================================\n\n";

    cout << "Sigma_total Z value:" << (sigma_tot*1e27 - 176.97) / (err_sigma_tot*1e27) << endl;
    cout << "Sigma_p Z value:" << (sigma_p*1e27 - 97.810) / (err_sigma_p*1e27) << endl;

    
}



