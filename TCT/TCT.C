#include <vector>
#include <string>
#include <cmath>
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TSystem.h"

void TCT()
{
   // ================= IMPOSTAZIONI =================
   const bool mostraCanvas   = true;   // true = apre anche le finestre, false = salva soltanto i file
   const bool erroriDaNoise  = false;   // true = usa A_noise_rms come errore su A nei plot A(V)

   // salva una canvas sia in .pdf che in .tex
   auto salva = [](TCanvas *c, const std::string &base) {
      c->SaveAs((base + ".pdf").c_str());
      c->SaveAs((base + ".tex").c_str());
   };

   // ================= DATI 1: A(z), scan 1 =================
   std::vector<double> z1 = { -3000, -2900, -2800, -2700, -2600, -2500, -2400, -2300, -2200, -2100,
                              -2000, -1900, -1800, -1700, -1600, -1500, -1400, -1300, -1200, -1100,
                              -1000, -900, -800, -700, -600, -500, -400, -300, -200, -100,
                              0, 100, 200, 300, 400, 500 };   // [um]
   std::vector<double> A_z1 = { 740, 764, 785, 805, 822, 836, 845, 850, 854, 853,
                                853, 855, 852, 848, 844, 835, 818, 793, 759, 724,
                                690, 664, 638, 623, 611, 601, 589, 578, 571, 565,
                                556.04, 550.15, 544.15, 540.01, 536.72, 532.05 };   // [pWb]
   std::vector<double> ez1(z1.size(), 0.), eA_z1(z1.size(), 0.);   // nessun errore

   // ================= DATI 2: A(z), scan 2 =================
   std::vector<double> z2   = { -2300, -2250, -2200, -2150, -2100, -2050, -2000, -1950, -1900, -1850, -1800 };   // [um]
   std::vector<double> A_z2 = { 546, 561, 579, 602, 617, 628, 625, 621, 606, 584, 544 };   // [pWb]
   std::vector<double> ez2(z2.size(), 0.), eA_z2(z2.size(), 0.);   // nessun errore

   // ================= DATI 3: Q(attenuazione) =================
   const double A_noise = 7.87;   // [pWb]
   std::vector<double> att = { 45.0, 45.5, 46.0, 46.5, 47.0, 47.5, 48.0, 48.5, 49.0, 49.5,
                               49.6, 49.7, 50.0, 50.1, 50.2, 50.3, 50.4, 50.5, 51.0 };   // [%]
   std::vector<double> A_att = { 485, 430, 361, 301, 243, 180, 118, 55.5, 23.3, 16.7,
                                 15.3, 13.9, 11.6, 11.7, 11.5, 10.8, 10.1, 10.0, 7.3 };   // [pWb]
   std::vector<double> Q, eatt(att.size(), 0.), eQ(att.size(), 0.);   // nessun errore
   for (size_t i = 0; i < A_att.size(); ++i)
      Q.push_back((A_att[i] - A_noise) / 4.7);   // [fC]

   // ================= DATI 4: A(V), set 1 =================
   std::vector<double> V1   = { 140, 150, 160, 170, 180, 190, 200, 210, 220,
                                230, 240, 250, 260, 270, 280, 290, 300 };   // [V]
   std::vector<double> A_V1 = { 63.3, 60.7, 66.1, 70.9, 72.1, 79.2, 81.5, 92.1, 95.8,
                                106.8, 119.9, 133.4, 145, 180.8, 213.6, 243.2, 275.8 };   // [pWb]
   std::vector<double> An_V1= { 11.2, 12.1, 11.9, 11.2, 11.8, 11.6, 12.0, 11.7, 11.6,
                                11.5, 11.1, 11.9, 12.4, 11.5, 12.3, 12.5, 11.8 };   // A_noise_rms [pWb]
   std::vector<double> eV1(V1.size(), 0.);
   std::vector<double> eA_V1 = erroriDaNoise ? An_V1 : std::vector<double>(V1.size(), 0.);

   // ================= DATI 5: A(V), set 2 =================
   std::vector<double> V2   = { 140, 160, 180, 200, 220, 240, 260, 280, 300 };   // [V]
   std::vector<double> A_V2 = { 332.9, 378.2, 394.3, 463.5, 550, 693.4, 858.1, 1074.4, 1704.7 };   // [pWb]
   std::vector<double> An_V2= { 12.7, 12.6, 13.2, 19.8, 18.7, 20.1, 19.4, 13.6, 13.4 };   // A_noise_rms [pWb]
   std::vector<double> eV2(V2.size(), 0.);
   std::vector<double> eA_V2 = erroriDaNoise ? An_V2 : std::vector<double>(V2.size(), 0.);

   // ================= GRAFICI =================
   // 1) A(z), scan 1
   TGraphErrors *g1 = new TGraphErrors(z1.size(), z1.data(), A_z1.data(), ez1.data(), eA_z1.data());
   g1->SetName("g_Az_1");
   g1->SetTitle("A vs z (scan 1);z [#mum];A [pWb]");
   g1->SetMarkerStyle(20);
   g1->SetMarkerSize(0.9);
   g1->SetMarkerColor(kBlack);
   g1->SetLineColor(kBlack);

   // 2) A(z), scan 2
   TGraphErrors *g2 = new TGraphErrors(z2.size(), z2.data(), A_z2.data(), ez2.data(), eA_z2.data());
   g2->SetName("g_Az_2");
   g2->SetTitle("A vs z (scan 2);z [#mum];A [pWb]");
   g2->SetMarkerStyle(21);
   g2->SetMarkerSize(0.9);
   g2->SetMarkerColor(kRed + 1);
   g2->SetLineColor(kRed + 1);

   // 3) Q(attenuazione)
   TGraphErrors *g3 = new TGraphErrors(att.size(), att.data(), Q.data(), eatt.data(), eQ.data());
   g3->SetName("g_Q_att");
   g3->SetTitle("Q vs attenuazione;Attenuazione [%];Q [fC]");
   g3->SetMarkerStyle(20);
   g3->SetMarkerSize(0.9);
   g3->SetMarkerColor(kBlack);
   g3->SetLineColor(kBlack);

   // 4) A(V), set 1
   TGraphErrors *g4 = new TGraphErrors(V1.size(), V1.data(), A_V1.data(), eV1.data(), eA_V1.data());
   g4->SetName("g_AV_1");
   g4->SetTitle("A vs V (set 1);V [V];A [pWb]");
   g4->SetMarkerStyle(20);
   g4->SetMarkerSize(0.9);
   g4->SetMarkerColor(kBlack);
   g4->SetLineColor(kBlack);

   // 5) A(V), set 2
   TGraphErrors *g5 = new TGraphErrors(V2.size(), V2.data(), A_V2.data(), eV2.data(), eA_V2.data());
   g5->SetName("g_AV_2");
   g5->SetTitle("A vs V (set 2);V [V];A [pWb]");
   g5->SetMarkerStyle(21);
   g5->SetMarkerSize(0.9);
   g5->SetMarkerColor(kRed + 1);
   g5->SetLineColor(kRed + 1);

   // ================= CANVAS =================
   TCanvas *c1 = new TCanvas("c_Az_1", "A vs z (scan 1)", 800, 600);
   c1->SetGrid();
   g1->Draw("AP");
   salva(c1, "./plots/TCT_A_z_1");

   TCanvas *c2 = new TCanvas("c_Az_2", "A vs z (scan 2)", 800, 600);
   c2->SetGrid();
   g2->Draw("AP");
   salva(c2, "./plots/TCT_A_z_2");

   TCanvas *c3 = new TCanvas("c_Q_att", "Q vs attenuazione", 800, 600);
   c3->SetGrid();
   g3->Draw("AP");
   salva(c3, "./plots/TCT_Q_att");

   TCanvas *c4 = new TCanvas("c_AV_1", "A vs V (set 1)", 800, 600);
   c4->SetGrid();
   g4->Draw("AP");
   salva(c4, "./plots/TCT_A_V_1");

   TCanvas *c5 = new TCanvas("c_AV_2", "A vs V (set 2)", 800, 600);
   c5->SetGrid();
   g5->Draw("AP");
   salva(c5, "./plots/TCT_A_V_2");
}