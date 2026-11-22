#include <vector>
#include <cmath>
#include "TGraphErrors.h"
#include "TCanvas.h"

void PiN()
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
   std::vector<double> V   = { 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.2,
                               1.4, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5,
                               6, 7, 8, 10, 12, 15, 20, 30, 40 };   // [-V]
   std::vector<double> C_v = { 21.6, 15.20, 9.34, 6.84, 6.07, 5.5, 5.18, 4.58, 4.17, 3.94,
                               3.75, 3.62, 3.49, 3.38, 3.19, 3.07, 3.00, 2.95, 2.91, 2.89,
                               2.85, 2.83, 2.81, 2.78, 2.77, 2.75, 2.72, 2.69, 2.67 };   // [pF]
   std::vector<double> eV(V.size(), 0.);   // nessun errore su V
   std::vector<double> eC_v= { 0.1, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                               0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                               0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05 };   // [pF]

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
   c1->SaveAs("./plots/PiN_C_f.png");

   TCanvas *c2 = new TCanvas("c2", "C vs V", 800, 600);
   c2->SetGrid();
   g2->Draw("AP");
   c2->SaveAs("./plots/PiN_C_V.png");

   TCanvas *c3 = new TCanvas("c3", "1/C^2 vs V", 800, 600);
   c3->SetGrid();
   g3->Draw("AP");
   c3->SaveAs("./plots/PiN_invC2_V.png");
}