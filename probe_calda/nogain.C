#include <vector>
#include <cmath>
#include "TGraphErrors.h"
#include "TCanvas.h"

void nogain()
{
   // ================= DATI 1: C - f =================
   std::vector<double> f   = { 0.5, 0.6, 0.7, 0.85, 1, 1.2, 1.45, 1.7, 2, 2.5,
                               3, 4, 5.3, 7.1, 9.5, 12, 15, 18, 22, 27,
                               35, 42, 55, 72, 90, 100, 120, 150, 200, 500 };   // [kHz]
   std::vector<double> C_f = { 97, 96, 96, 95, 95, 94.5, 94, 93.3, 93.3, 92.7,
                               91.8, 90.5, 89.2, 85.9, 82.3, 77.4, 71.6, 65.6, 59.15, 52.0,
                               43.40, 38.18, 32.08, 27.76, 25.3, 24.3, 23.0, 21.9, 21.1, 19.90 };   // [pF]
   std::vector<double> ef(f.size(), 0.);   // da capire errore su f
   std::vector<double> eC_f= { 1, 1, 1, 1, 1, 0.5, 0.5, 0.3, 0.1, 0.2,
                               0.2, 0.3, 0.1, 0.2, 0.2, 0.1, 0.1, 0.1, 0.05, 0.1,
                               0.05, 0.05, 0.05, 0.05, 0.2, 0.1, 0.05, 0.1, 0.1, 0.1 };   // [pF]

   // ================= DATI 2: C - V =================
   std::vector<double> V   = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45,
                               50, 55, 60, 65, 70, 75, 80, 85, 90, 95,
                               100, 105, 110, 115, 120, 125, 130, 135, 140 };   // [V]
   std::vector<double> C_v = { 134, 104, 86.3, 81.5, 78.4, 76.2, 74.1, 72.2, 70.0, 64.6,
                               58.9, 56.9, 55.0, 53.6, 52.2, 51.2, 50.2, 49.4, 48.6, 48.0,
                               47.5, 47.2, 47.1, 47.0, 46.9, 46.9, 46.8, 46.8, 46.8 };   // [pF]
   std::vector<double> eV(V.size(), 0.);   // da capire errore su V
   std::vector<double> eC_v= { 1, 1, 0.5, 0.5, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                               0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                               0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2 };   // [pF]

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
   g2->SetTitle("C vs V;V [V];C [pF]");
   g2->SetMarkerStyle(20);
   g2->SetMarkerSize(0.9);

   // 3) 1/C^2 - V
   TGraphErrors *g3 = new TGraphErrors(V.size(), V.data(), invC2.data(), eV.data(), e_invC2.data());
   g3->SetName("g_invC2V");
   g3->SetTitle("1/C^{2} vs V;V [V];1/C^{2} [pF^{-2}]");
   g3->SetMarkerStyle(20);
   g3->SetMarkerSize(0.9);

   TCanvas *c1 = new TCanvas("c1", "C vs f", 800, 600);
   c1->SetGrid(); c1->SetLogx();   // f in scala log
   g1->Draw("AP");
   c1->SaveAs("./plots/C_f.pdf");
   c1->SaveAs("./plots/C_f.tex");

   TCanvas *c2 = new TCanvas("c2", "C vs V", 800, 600);
   c2->SetGrid();
   g2->Draw("AP");
   c2->SaveAs("./plots/C_V.pdf");
   c2->SaveAs("./plots/C_V.tex");

   TCanvas *c3 = new TCanvas("c3", "1/C^2 vs V", 800, 600);
   c3->SetGrid();
   g3->Draw("AP");
   c3->SaveAs("./plots/1_C2_V.pdf");
   c3->SaveAs("./plots/1_C2_V.tex");
}