// Uso: root -l LGAD.C
#include <vector>
#include <cmath>
#include "TGraphErrors.h"
#include "TCanvas.h"

void LGAD()
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
   std::vector<double> V   = { 0, 2, 4, 6, 8, 10, 12, 14, 16, 18,
                               19, 20, 20.5, 21, 21.5, 22, 22.5, 23, 23.5, 24,
                               24.5, 25, 25.5, 26, 27, 28, 29, 30, 31, 32,
                               34, 36, 38, 40 };   // [-V]
   std::vector<double> C_v = { 307.2, 220.35, 190.11, 172.22, 159.35, 149, 140.29, 132.41, 124.93, 117.36,
                               113.36, 108.91, 106.51, 103.89, 100.97, 97.56, 93.32, 87.3, 75.92, 45.5,
                               14.17, 4.81, 3.82, 3.32, 2.99, 2.84, 2.79, 2.76, 2.74, 2.7,
                               2.72, 2.72, 2.71, 2.7 };   // [pF]
   std::vector<double> eV(V.size(), 0.);   // nessun errore su V
   std::vector<double> eC_v= { 0.1, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                               0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.3,
                               0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                               0.05, 0.05, 0.05, 0.05 };   // [pF]

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
   c1->SaveAs("./plots/LGAD_C_f.png");

   TCanvas *c2 = new TCanvas("c2", "C vs V", 800, 600);
   c2->SetGrid();
   g2->Draw("AP");
   c2->SaveAs("./plots/LGAD_C_V.png");

   TCanvas *c3 = new TCanvas("c3", "1/C^2 vs V", 800, 600);
   c3->SetGrid();
   g3->Draw("AP");
   c3->SaveAs("./plots/LGAD_invC2_V.png");
}