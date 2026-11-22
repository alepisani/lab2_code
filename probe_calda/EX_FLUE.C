#include <vector>
#include <cmath>
#include "TGraphErrors.h"
#include "TCanvas.h"

void EX_FLUE()
{
   // ================= DATI 1: C - f =================
   std::vector<double> f   = { 1, 1.5, 2.3, 3.5, 5.3, 8, 12, 18, 27, 40,
                               60, 90, 140, 210, 300 };   // [kHz]
   std::vector<double> C_f = { 151.5, 150.7, 150.1, 150.2, 150.1, 149.93, 148.5, 148.6, 147.9, 145.4,
                               138.7, 129.1, 106, 80.7, 55.2 };   // [pF]
   std::vector<double> ef(f.size(), 0.);   // nessun errore su f
   std::vector<double> eC_f= { 0.2, 0.1, 0.1, 0.1, 0.1, 0.05, 0.3, 0.1, 0.1, 0.1,
                               0.2, 0.2, 1, 0.2, 0.2 };   // [pF]

   // ================= DATI 2: C - V =================
   std::vector<double> V   = { 2, 4, 6, 8, 10, 12, 14, 16, 18,
                               20, 22, 24, 25, 26, 26.5, 27, 27.5, 28, 28.5,
                               29, 29.5, 30, 30.5, 31, 31.5, 32, 32.5, 33, 34,
                               35, 36, 37, 38, 39 };   // [-V]
   std::vector<double> C_v = { 19.71, 19.22, 18.86, 18.68, 18.53, 18.41, 18.3, 18.2, 18.1,
                               18.01, 17.93, 17.62, 17.46, 16.38, 13.95, 12.25, 11.0, 10.05, 9.3,
                               8.67, 8.16, 7.72, 7.36, 7.06, 6.9, 6.75, 6.66, 6.61, 6.59,
                               6.57, 6.5, 6.45, 6.40, 6.39 };   // [pF]
   std::vector<double> eV(V.size(), 0.);   // nessun errore su V
   std::vector<double> eC_v= { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                               0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                               0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                               0.05, 0.1, 0.1, 0.1, 0.1 };   // [pF]

   // ================= 1/C^2 - V =================
   std::vector<double> invC2, e_invC2;
   for (size_t i = 0; i < C_v.size(); ++i) {
      invC2.push_back(1. / (C_v[i] * C_v[i]));
      e_invC2.push_back(2. * eC_v[i] / (C_v[i] * C_v[i] * C_v[i]));   // propagazione errore
   }

   // ================= GRAFICI =================
   // 1) C - f
   TGraphErrors *g1 = new TGraphErrors(f.size(), f.data(), C_f.data(), ef.data(), eC_f.data());
   g1->SetName("g_Cf");
   g1->SetTitle("C vs f;f [kHz];C [pF]");
   g1->SetMarkerStyle(20);
   g1->SetMarkerSize(0.9);

   // 2) C - V
   TGraphErrors *g2 = new TGraphErrors(V.size(), V.data(), C_v.data(), eV.data(), eC_v.data());
   g2->SetName("g_CV");
   g2->SetTitle("C vs V;-V [V];C [pF]");
   g2->SetMarkerStyle(20);
   g2->SetMarkerSize(0.9);

   // 3) 1/C^2 - V
   TGraphErrors *g3 = new TGraphErrors(V.size(), V.data(), invC2.data(), eV.data(), e_invC2.data());
   g3->SetName("g_invC2V");
   g3->SetTitle("1/C^{2} vs V;-V [V];1/C^{2} [pF^{-2}]");
   g3->SetMarkerStyle(20);
   g3->SetMarkerSize(0.9);

   TCanvas *c1 = new TCanvas("c1", "C vs f", 800, 600);
   c1->SetGrid(); c1->SetLogx();   // f in scala log
   g1->Draw("AP");
   c1->SaveAs("./plots/EX-FLUE_C_f.png");

   TCanvas *c2 = new TCanvas("c2", "C vs V", 800, 600);
   c2->SetGrid();
   g2->Draw("AP");
   c2->SaveAs("./plots/EX-FLUE_C_V.png");

   TCanvas *c3 = new TCanvas("c3", "1/C^2 vs V", 800, 600);
   c3->SetGrid();
   g3->Draw("AP");
   c3->SaveAs("./plots/EX-FLUE_invC2_V.png");
}