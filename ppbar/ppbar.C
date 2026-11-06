#include <iostream>
#include <cmath>



using namespace std;

void ppbar()
{
    int iNp = 10;
    int iNd = 10;
    int iN_free = 10;

    double L1_mis = 1.;                 //cm
    double L1_err = 0.1;                //cm
    double L1_true = 1.;                //cm
    double L1_true_err = 0.1;           //cm

    double L2_mis = 1.;                 //cm
    double L2_err = 0.1;                //cm
    double L2_true = 1.;                //cm
    double L2_true_err = 0.1;           //cm

    double L_mis = 57.5;                //cm
    double L_mis_err = 0.2;             //cm

    double K = 2.039e-23;               //cm3

    double Np = (double)iNp;
    double Nd = (double)iNd;
    double N_free = (double)iN_free;

    double err_Np = sqrt(Np);
    double err_Nd = sqrt(Nd);
    double err_Nfree = sqrt(N_free);

    double N0 = Np + Nd + N_free;
    double err_N0 = sqrt(pow(err_Np,2) + pow(err_Nd,2) + pow(err_Nfree,2));

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

    double logtot = log(Np+Nd+N_free) - log(N_free);
    double err_logtot = 1./N0 * sqrt(Np + Nd + pow(Np+Nd,2)/N_free);
    double sigma_tot = (K/Ltrue) * logtot;
    double err_sigma_tot = sigma_tot * sqrt(pow(Ltrue_err / Ltrue, 2) + pow(err_logtot/logtot, 2));

    double logp = log(N0) - log(N_free + Nd);
    double err_logp = 1./N0 * sqrt(Np + Np*Np/(Nd+N_free));
    double sigma_p = K/Ltrue * logp;
    double err_sigma_p = sigma_p * sqrt(pow(Ltrue_err / Ltrue, 2) + pow(err_logp/logp, 2));

    double logn = log(N0) - log(N0 - Nd);
    double err_logn = 1./N0 * sqrt(Nd + Nd*Nd/(Np+N_free));
    double sigma_n = K/Ltrue * logn;
    double err_sigma_n = sigma_n * sqrt(pow(Ltrue_err / Ltrue, 2) + pow(err_logn/logn, 2));

    cout << "\n==========================================================\n\n";
    cout << "Sigma_total (mb): " << sigma_tot*1e27 << " +/- " << err_sigma_tot*1e27 << "   (" << (err_sigma_tot / sigma_tot * 100.) << "%)" << endl;
    cout << "Sigma_p (mb): " << sigma_p*1e27 << " +/- " << err_sigma_p*1e27 << "   (" << (err_sigma_p / sigma_p * 100.) << "%)" << endl;
    cout << "Sigma_n (mb): " << sigma_n*1e27 << " +/- " << err_sigma_n*1e27 << "   (" << (err_sigma_n / sigma_n * 100.) << "%)" << endl;
    cout << "\n==========================================================\n\n";
}



