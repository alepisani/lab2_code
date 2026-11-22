#include <vector>
#include <string>
#include <cmath>
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TAxis.h"
#include "TSystem.h"

struct Dataset {
   std::string name;
   int color, marker;
   std::vector<double> f, C_f, eC_f;   // C - f   (f in kHz, C in pF)
   std::vector<double> V, C_v, eC_v;   // C - V   (|V| in V, C in pF)
};

void chips()
{
   // ================= IMPOSTAZIONI =================
   const double Vmax = 40.;           // massimo |V| mostrato in tutti i plot con la tensione
   const bool mostraCanvas = true;   // true = apre anche le finestre, false = salva soltanto i file

   if (!mostraCanvas) gROOT->SetBatch(kTRUE);
   gSystem->mkdir("plots", kTRUE);   // crea la cartella plots se non esiste

   std::vector<Dataset> sets(4);

   // ================= nogain =================
   sets[0].name   = "nogain";
   sets[0].color  = kBlack;
   sets[0].marker = 20;
   sets[0].f   = { 0.5, 0.6, 0.7, 0.85, 1, 1.2, 1.45, 1.7, 2, 2.5,
                   3, 4, 5.3, 7.1, 9.5, 12, 15, 18, 22, 27,
                   35, 42, 55, 72, 90, 100, 120, 150, 200, 500 };
   sets[0].C_f = { 97, 96, 96, 95, 95, 94.5, 94, 93.3, 93.3, 92.7,
                   91.8, 90.5, 89.2, 85.9, 82.3, 77.4, 71.6, 65.6, 59.15, 52.0,
                   43.40, 38.18, 32.08, 27.76, 25.3, 24.3, 23.0, 21.9, 21.1, 19.90 };
   sets[0].eC_f= { 1, 1, 1, 1, 1, 0.5, 0.5, 0.3, 0.1, 0.2,
                   0.2, 0.3, 0.1, 0.2, 0.2, 0.1, 0.1, 0.1, 0.05, 0.1,
                   0.05, 0.05, 0.05, 0.05, 0.2, 0.1, 0.05, 0.1, 0.1, 0.1 };
   sets[0].V   = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45,
                   50, 55, 60, 65, 70, 75, 80, 85, 90, 95,
                   100, 105, 110, 115, 120, 125, 130, 135, 140 };
   sets[0].C_v = { 134, 104, 86.3, 81.5, 78.4, 76.2, 74.1, 72.2, 70.0, 64.6,
                   58.9, 56.9, 55.0, 53.6, 52.2, 51.2, 50.2, 49.4, 48.6, 48.0,
                   47.5, 47.2, 47.1, 47.0, 46.9, 46.9, 46.8, 46.8, 46.8 };
   sets[0].eC_v= { 1, 1, 0.5, 0.5, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                   0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                   0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2 };

   // ================= LGAD =================
   sets[1].name   = "LGAD";
   sets[1].color  = kRed + 1;
   sets[1].marker = 21;
   sets[1].f   = { 1, 1.5, 2.3, 3.5, 5.3, 8, 12, 18, 27, 40,
                   60, 90, 140, 210, 300 };
   sets[1].C_f = { 151.5, 150.7, 150.1, 150.2, 150.1, 149.93, 148.5, 148.6, 147.9, 145.4,
                   138.7, 129.1, 106, 80.7, 55.2 };
   sets[1].eC_f= { 0.2, 0.1, 0.1, 0.1, 0.1, 0.05, 0.3, 0.1, 0.1, 0.1,
                   0.2, 0.2, 1, 0.2, 0.2 };
   sets[1].V   = { 0, 2, 4, 6, 8, 10, 12, 14, 16, 18,
                   19, 20, 20.5, 21, 21.5, 22, 22.5, 23, 23.5, 24,
                   24.5, 25, 25.5, 26, 27, 28, 29, 30, 31, 32,
                   34, 36, 38, 40 };
   sets[1].C_v = { 307.2, 220.35, 190.11, 172.22, 159.35, 149, 140.29, 132.41, 124.93, 117.36,
                   113.36, 108.91, 106.51, 103.89, 100.97, 97.56, 93.32, 87.3, 75.92, 45.5,
                   14.17, 4.81, 3.82, 3.32, 2.99, 2.84, 2.79, 2.76, 2.74, 2.7,
                   2.72, 2.72, 2.71, 2.7 };
   sets[1].eC_v= { 0.1, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                   0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.3,
                   0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                   0.05, 0.05, 0.05, 0.05 };

   // ================= PiN =================
   sets[2].name   = "PiN";
   sets[2].color  = kBlue + 1;
   sets[2].marker = 22;
   sets[2].f   = { 1, 1.5, 2.3, 3.5, 5.3, 8, 12, 18, 27, 40,
                   60, 90, 140, 210, 300 };
   sets[2].C_f = { 151.5, 150.7, 150.1, 150.2, 150.1, 149.93, 148.5, 148.6, 147.9, 145.4,
                   138.7, 129.1, 106, 80.7, 55.2 };
   sets[2].eC_f= { 0.2, 0.1, 0.1, 0.1, 0.1, 0.05, 0.3, 0.1, 0.1, 0.1,
                   0.2, 0.2, 1, 0.2, 0.2 };
   sets[2].V   = { 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.2,
                   1.4, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5,
                   6, 7, 8, 10, 12, 15, 20, 30, 40 };
   sets[2].C_v = { 21.6, 15.20, 9.34, 6.84, 6.07, 5.5, 5.18, 4.58, 4.17, 3.94,
                   3.75, 3.62, 3.49, 3.38, 3.19, 3.07, 3.00, 2.95, 2.91, 2.89,
                   2.85, 2.83, 2.81, 2.78, 2.77, 2.75, 2.72, 2.69, 2.67 };
   sets[2].eC_v= { 0.1, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                   0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                   0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05 };

   // ================= EX-FLUE =================
   sets[3].name   = "EX-FLUE";
   sets[3].color  = kGreen + 2;
   sets[3].marker = 23;
   sets[3].f   = { 1, 1.5, 2.3, 3.5, 5.3, 8, 12, 18, 27, 40,
                   60, 90, 140, 210, 300 };
   sets[3].C_f = { 151.5, 150.7, 150.1, 150.2, 150.1, 149.93, 148.5, 148.6, 147.9, 145.4,
                   138.7, 129.1, 106, 80.7, 55.2 };
   sets[3].eC_f= { 0.2, 0.1, 0.1, 0.1, 0.1, 0.05, 0.3, 0.1, 0.1, 0.1,
                   0.2, 0.2, 1, 0.2, 0.2 };
   sets[3].V   = { 2, 4, 6, 8, 10, 12, 14, 16, 18,
                   20, 22, 24, 25, 26, 26.5, 27, 27.5, 28, 28.5,
                   29, 29.5, 30, 30.5, 31, 31.5, 32, 32.5, 33, 34,
                   35, 36, 37, 38, 39 };
   sets[3].C_v = { 19.71, 19.22, 18.86, 18.68, 18.53, 18.41, 18.3, 18.2, 18.1,
                   18.01, 17.93, 17.62, 17.46, 16.38, 13.95, 12.25, 11.0, 10.05, 9.3,
                   8.67, 8.16, 7.72, 7.36, 7.06, 6.9, 6.75, 6.66, 6.61, 6.59,
                   6.57, 6.5, 6.45, 6.40, 6.39 };
   sets[3].eC_v= { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                   0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                   0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                   0.05, 0.1, 0.1, 0.1, 0.1 };

   // ================= FUNZIONI DI SERVIZIO =================
   // salva una canvas sia in .pdf che in .tex
   auto salva = [](TCanvas *c, const std::string &base) {
      c->SaveAs((base + ".pdf").c_str());
      c->SaveAs((base + ".tex").c_str());
   };

   // ================= MULTIGRAPH E LEGENDE (confronto) =================
   TMultiGraph *mg1 = new TMultiGraph("mg_Cf",     "Confronto C vs f;f [kHz];C [pF]");
   TMultiGraph *mg2 = new TMultiGraph("mg_CV",     "Confronto C vs V;|V| [V];C [pF]");
   TMultiGraph *mg3 = new TMultiGraph("mg_invC2V", "Confronto 1/C^{2} vs V;|V| [V];1/C^{2} [pF^{-2}]");

   // coordinate NDC (x1, y1, x2, y2) della legenda; gli assi occupano da 0.1 a 0.9
   TLegend *leg1 = new TLegend(0.11, 0.11, 0.27, 0.29);   // C-f      : basso a sinistra
   TLegend *leg2 = new TLegend(0.73, 0.71, 0.89, 0.89);   // C-V      : alto a destra
   TLegend *leg3 = new TLegend(0.73, 0.11, 0.89, 0.29);   // 1/C^2-V  : basso a destra
   TLegend *legs[3] = { leg1, leg2, leg3 };
   for (int j = 0; j < 3; ++j) {
      legs[j]->SetFillColor(kWhite);
      legs[j]->SetFillStyle(1001);
      legs[j]->SetBorderSize(1);
   }

   // ================= CICLO SUI SENSORI =================
   for (size_t k = 0; k < sets.size(); ++k) {
      Dataset &s = sets[k];
      int nf = s.f.size();

      // tengo solo i punti con |V| <= Vmax
      std::vector<double> Vc, Cc, eCc, invC2, e_invC2;
      for (size_t i = 0; i < s.V.size(); ++i) {
         if (s.V[i] > Vmax) continue;
         Vc.push_back(s.V[i]);
         Cc.push_back(s.C_v[i]);
         eCc.push_back(s.eC_v[i]);
         invC2.push_back(1. / (s.C_v[i] * s.C_v[i]));
         e_invC2.push_back(2. * s.eC_v[i] / (s.C_v[i] * s.C_v[i] * s.C_v[i]));   // propagazione errore
      }
      int nv = Vc.size();
      std::vector<double> ef(nf, 0.), eV(nv, 0.);   // nessun errore su f e V

      // --- 1) C - f ---
      TGraphErrors *g1 = new TGraphErrors(nf, s.f.data(), s.C_f.data(), ef.data(), s.eC_f.data());
      g1->SetName(("g_Cf_" + s.name).c_str());
      g1->SetTitle((s.name + ": C vs f;f [kHz];C [pF]").c_str());

      // --- 2) C - V ---
      TGraphErrors *g2 = new TGraphErrors(nv, Vc.data(), Cc.data(), eV.data(), eCc.data());
      g2->SetName(("g_CV_" + s.name).c_str());
      g2->SetTitle((s.name + ": C vs V;|V| [V];C [pF]").c_str());

      // --- 3) 1/C^2 - V ---
      TGraphErrors *g3 = new TGraphErrors(nv, Vc.data(), invC2.data(), eV.data(), e_invC2.data());
      g3->SetName(("g_invC2V_" + s.name).c_str());
      g3->SetTitle((s.name + ": 1/C^{2} vs V;|V| [V];1/C^{2} [pF^{-2}]").c_str());

      // stesso colore e stesso marker usati poi nel confronto
      TGraphErrors *gs[3] = { g1, g2, g3 };
      for (int j = 0; j < 3; ++j) {
         gs[j]->SetMarkerStyle(s.marker);
         gs[j]->SetMarkerSize(0.9);
         gs[j]->SetMarkerColor(s.color);
         gs[j]->SetLineColor(s.color);
      }

      // --- canvas singole ---
      TCanvas *ci1 = new TCanvas(("c_" + s.name + "_Cf").c_str(), (s.name + ": C vs f").c_str(), 800, 600);
      ci1->SetGrid(); ci1->SetLogx();   // f in scala log
      g1->Draw("AP");
      salva(ci1, "./plots/" + s.name + "_C_f");

      TCanvas *ci2 = new TCanvas(("c_" + s.name + "_CV").c_str(), (s.name + ": C vs V").c_str(), 800, 600);
      ci2->SetGrid();
      g2->Draw("AP");
      ci2->Update();
      g2->GetXaxis()->SetLimits(0., Vmax);
      ci2->Modified(); ci2->Update();
      salva(ci2, "./plots/" + s.name + "_C_V");

      TCanvas *ci3 = new TCanvas(("c_" + s.name + "_invC2V").c_str(), (s.name + ": 1/C^2 vs V").c_str(), 800, 600);
      ci3->SetGrid();
      g3->Draw("AP");
      ci3->Update();
      g3->GetXaxis()->SetLimits(0., Vmax);
      ci3->Modified(); ci3->Update();
      salva(ci3, "./plots/" + s.name + "_invC2_V");

      // --- copie dei grafici per il confronto (stessi stili) ---
      TGraphErrors *m1 = (TGraphErrors*)g1->Clone();
      TGraphErrors *m2 = (TGraphErrors*)g2->Clone();
      TGraphErrors *m3 = (TGraphErrors*)g3->Clone();
      mg1->Add(m1, "P");
      mg2->Add(m2, "P");
      mg3->Add(m3, "P");
      leg1->AddEntry(m1, s.name.c_str(), "p");
      leg2->AddEntry(m2, s.name.c_str(), "p");
      leg3->AddEntry(m3, s.name.c_str(), "p");
   }

   // ================= CANVAS DI CONFRONTO =================
   TCanvas *c1 = new TCanvas("c_conf_Cf", "Confronto C vs f", 800, 600);
   c1->SetGrid(); c1->SetLogx();   // f in scala log
   mg1->Draw("A");
   leg1->Draw();
   salva(c1, "./plots/confronto_C_f");

   TCanvas *c2 = new TCanvas("c_conf_CV", "Confronto C vs V", 800, 600);
   c2->SetGrid(); c2->SetLogy();   // togli SetLogy() per la scala lineare
   mg2->Draw("A");
   c2->Update();
   mg2->GetXaxis()->SetLimits(0., Vmax);
   leg2->Draw();
   c2->Modified(); c2->Update();
   salva(c2, "./plots/confronto_C_V");

   TCanvas *c3 = new TCanvas("c_conf_invC2V", "Confronto 1/C^2 vs V", 800, 600);
   c3->SetGrid(); c3->SetLogy();   // togli SetLogy() per la scala lineare
   mg3->Draw("A");
   c3->Update();
   mg3->GetXaxis()->SetLimits(0., Vmax);
   leg3->Draw();
   c3->Modified(); c3->Update();
   salva(c3, "./plots/confronto_invC2_V");
}