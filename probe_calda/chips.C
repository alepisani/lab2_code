#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TAxis.h"
#include "TSystem.h"
#include "TF1.h"

void chips()
{
   // =============== PHYSICAL VALUES ================
   const double A_nogain = 2.4;    //cm2
   const double A_gain = 0.01;     //cm2
   const double errRelA  = 0.05;   // errore relativo sull'area
   const double eps_Si   = 11.7 * 8.8541878128e-14;   // F/cm  (eps_r = 11.7)
   const double q_el     = 1.602176634e-19;           // C



   // ================= IMPOSTAZIONI =================
   // scegli su quali chip fare l'analisi (true = analizza, false = salta)
   const bool chip_nogain  = false;
   const bool chip_LGAD    = false;
   const bool chip_PiN     = false;
   const bool chip_EXFLUE  = true;

   const bool fareConfronto = false;   // true = fa i plot di confronto (con i chip attivi sopra)
   const double Vmax = 40.;           // massimo |V| mostrato in tutti i plot con la tensione
   const bool mostraCanvas = true;    // true = apre anche le finestre, false = salva soltanto i file
   const bool salva_plots = false;

   // ---- fit lineari su 1/C^2 - V ----
   const bool fareFit = true;         // true = fa i 3 fit lineari su ogni plot 1/C^2 - V
   const bool fareProfilo = true;     // true = fa il plot del profilo di drogaggio N_B - W per ogni chip
   const int colFitA = kMagenta + 1;  // colori delle tre rette
   const int colFitB = kOrange + 7;
   const int colFitC = kCyan + 2;
   // I range dei fit (in V) sono dentro il blocco di ogni chip, vicino al plot 1/C^2 - V:
   //    fitA_min/fitA_max, fitB_min/fitB_max, fitC_min/fitC_max

   if (!mostraCanvas) gROOT->SetBatch(kTRUE);
   gSystem->mkdir("plots", kTRUE);

   // salva una canvas sia in .pdf che in .tex
   auto salva = [](TCanvas *c, const std::string &base) {
      if(salva_plots){
         c->SaveAs((base + ".pdf").c_str());
         c->SaveAs((base + ".tex").c_str());
      }
   };

   // fit lineare (y = p0 + p1*x) di un grafico nel range [xmin, xmax]:
   // stampa i risultati nel terminale e restituisce la retta (da disegnare con "same")
   auto fitLineare = [](TGraphErrors *g, double xmin, double xmax, int col, const char *tag) -> TF1* {
      TF1 *fit = new TF1(Form("fit_%s_%s", g->GetName(), tag), "pol1", xmin, xmax);
      fit->SetLineColor(col);
      fit->SetLineWidth(2);
      g->Fit(fit, "RNQ");   // R = usa il range della funzione, N = non disegna, Q = silenzioso
      double chi2ndf = (fit->GetNDF() > 0) ? fit->GetChisquare() / fit->GetNDF() : 0.;
      printf("[%s] fit %s  V in [%g, %g]:  p0 = %g +- %g,  p1 = %g +- %g,  chi2/ndf = %g\n",
             g->GetName(), tag, xmin, xmax,
             fit->GetParameter(0), fit->GetParError(0),
             fit->GetParameter(1), fit->GetParError(1), chi2ndf);
      return fit;
   };

   // concentrazione di droganti da una retta di fit su 1/C^2 - V:
   //   Nb = 2 / (eps_Si * A^2 * q * pendenza)
   // la pendenza p1 del fit e' in [pF^-2 / V] -> la converto in [F^-2 / V] (x 1e24); A in cm^2 -> Nb in cm^-3
   // errore: (dNb/Nb)^2 = (2 dA/A)^2 + (dpendenza/pendenza)^2
   auto stampaNb = [&](TF1 *fit, double A, const std::string &chip, const char *tag) {
      double slope  = fit->GetParameter(1) * 1e24;
      double eslope = fit->GetParError(1)  * 1e24;
      if (slope == 0.) {
         printf("[%s] Nb (fit %s): pendenza nulla, Nb non calcolabile\n", chip.c_str(), tag);
         return;
      }
      double Nb    = 2. / (eps_Si * A * A * q_el * std::fabs(slope));
      double relEr = std::sqrt(std::pow(2. * errRelA, 2) + std::pow(eslope / slope, 2));
      printf("[%s] Nb (fit %s) = %.3e +- %.1e cm^-3   (A = %g cm^2)\n",
             chip.c_str(), tag, Nb, Nb * relEr, A);
   };

   // profilo di drogaggio N_B(W) da C - V
   //   pendenza tra il punto i e il punto i+salto:  s = [1/C^2(j) - 1/C^2(i)] / [V(j) - V(i)]   (pF^-2 / V)
   //   N_B = 2 / (eps_Si * A^2 * q * s)                       [cm^-3]
   //   W   = eps_Si * A / C, calcolata come media dei due punti [um]
   // propagazione degli errori (V senza errore):
   //   e(1/C^2) = 2 eC / C^3,    e(s) = sqrt(e1^2 + e2^2) / dV
   //   (eN/N)^2 = (e(s)/s)^2 + (2 eA/A)^2
   //   eW^2     = [ 1/2 sqrt( (W1 eC1/C1)^2 + (W2 eC2/C2)^2 ) ]^2 + (W eA/A)^2
   // i punti con pendenza <= 0 (rumore) vengono scartati
   auto profiloDrogaggio = [&](const std::vector<double> &V, const std::vector<double> &C, const std::vector<double> &eC,
                               double A, int salto, const std::string &chip, int col, int mk) -> TGraphErrors* {
      std::vector<double> W, eW, Nb, eNb;
      for (size_t i = 0; i + salto < V.size(); ++i) {
         size_t j = i + salto;
         double dV = V[j] - V[i];
         if (dV <= 0.) continue;
         double y1  = 1. / (C[i] * C[i]),  y2  = 1. / (C[j] * C[j]);
         double ey1 = 2. * eC[i] / (C[i] * C[i] * C[i]);
         double ey2 = 2. * eC[j] / (C[j] * C[j] * C[j]);
         double s   = (y2 - y1) / dV;                          // [pF^-2 / V]
         double es  = std::sqrt(ey1 * ey1 + ey2 * ey2) / dV;
         if (s <= 0.) continue;
         double N   = 2. / (eps_Si * A * A * q_el * s * 1e24); // 1e24: pF^-2 -> F^-2
         double eN  = N * std::sqrt(std::pow(es / s, 2) + std::pow(2. * errRelA, 2));
         double W1  = eps_Si * A / (C[i] * 1e-12);             // [cm]
         double W2  = eps_Si * A / (C[j] * 1e-12);
         double Wm  = 0.5 * (W1 + W2);
         double eWs = 0.5 * std::sqrt(std::pow(W1 * eC[i] / C[i], 2) + std::pow(W2 * eC[j] / C[j], 2));
         double eWm = std::sqrt(eWs * eWs + std::pow(Wm * errRelA, 2));
         W.push_back(Wm * 1e4);  eW.push_back(eWm * 1e4);      // cm -> um
         Nb.push_back(N);        eNb.push_back(eN);
      }
      TGraphErrors *g = new TGraphErrors(W.size(), W.data(), Nb.data(), eW.data(), eNb.data());
      g->SetName(("g_profilo_" + chip).c_str());
      g->SetTitle((chip + ": profilo di drogaggio;W [#mum];N_{B} [cm^{-3}]").c_str());
      g->SetMarkerStyle(mk);
      g->SetMarkerSize(0.9);
      g->SetMarkerColor(col);
      g->SetLineColor(col);
      return g;
   };

   // ================= MULTIGRAPH E LEGENDE (confronto) =================
   TMultiGraph *mg1 = new TMultiGraph("mg_Cf",     "Confronto C vs f;f [kHz];C [pF]");
   TMultiGraph *mg2 = new TMultiGraph("mg_CV",     "Confronto C vs V;|V| [V];C [pF]");
   TMultiGraph *mg3 = new TMultiGraph("mg_invC2V", "Confronto 1/C^{2} vs V;|V| [V];1/C^{2} [pF^{-2}]");

   // coordinate NDC (x1, y1, x2, y2); gli assi occupano da 0.1 a 0.9
   TLegend *leg1 = new TLegend(0.11, 0.11, 0.27, 0.23);   // C-f     : basso a sinistra
   TLegend *leg2 = new TLegend(0.73, 0.74, 0.89, 0.89);   // C-V     : alto a destra
   TLegend *leg3 = new TLegend(0.73, 0.11, 0.89, 0.26);   // 1/C^2-V : basso a destra
   TLegend *legs[3] = { leg1, leg2, leg3 };
   for (int j = 0; j < 3; ++j) {
      legs[j]->SetFillColor(kWhite);
      legs[j]->SetFillStyle(1001);
      legs[j]->SetBorderSize(1);
   }

   // =====================================================================
   // ============================  nogain  ===============================
   // =====================================================================
   if (chip_nogain) {
      std::string name = "nogain";
      int color = kBlack, marker = 20;
      const double A = A_nogain;   // area usata per Nb [cm^2]
      const int nSalto = 3;   // derivata tra il punto i e il punto i+nSalto (aumentalo per ridurre il rumore)

      // range dei fit lineari su 1/C^2 - V [V]  (CAMBIALI TU)
      const double fitA_min = 19,   fitA_max = 39;
      const double fitB_min = 60,  fitB_max = 97;
      const double fitC_min = 101,  fitC_max = 140;

      std::vector<double> f   = { 0.5, 0.6, 0.7, 0.85, 1, 1.2, 1.45, 1.7, 2, 2.5,
                                  3, 4, 5.3, 7.1, 9.5, 12, 15, 18, 22, 27,
                                  35, 42, 55, 72, 90, 100, 120, 150, 200, 500 };
      std::vector<double> C_f = { 97, 96, 96, 95, 95, 94.5, 94, 93.3, 93.3, 92.7,
                                  91.8, 90.5, 89.2, 85.9, 82.3, 77.4, 71.6, 65.6, 59.15, 52.0,
                                  43.40, 38.18, 32.08, 27.76, 25.3, 24.3, 23.0, 21.9, 21.1, 19.90 };
      std::vector<double> eC_f= { 1, 1, 1, 1, 1, 0.5, 0.5, 0.3, 0.1, 0.2,
                                  0.2, 0.3, 0.1, 0.2, 0.2, 0.1, 0.1, 0.1, 0.05, 0.1,
                                  0.05, 0.05, 0.05, 0.05, 0.2, 0.1, 0.05, 0.1, 0.1, 0.1 };
      std::vector<double> V   = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45,
                                  50, 55, 60, 65, 70, 75, 80, 85, 90, 95,
                                  100, 105, 110, 115, 120, 125, 130, 135, 140 };
      std::vector<double> C_v = { 134, 104, 86.3, 81.5, 78.4, 76.2, 74.1, 72.2, 70.0, 64.6,
                                  58.9, 56.9, 55.0, 53.6, 52.2, 51.2, 50.2, 49.4, 48.6, 48.0,
                                  47.5, 47.2, 47.1, 47.0, 46.9, 46.9, 46.8, 46.8, 46.8 };
      std::vector<double> eC_v= { 1, 1, 0.5, 0.5, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                                  0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                                  0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2 };

      // (qui nogain usa tutti i punti, senza il taglio su Vmax)
      std::vector<double> Vc, Cc, eCc, invC2, e_invC2;
      for (size_t i = 0; i < V.size(); ++i) {
         Vc.push_back(V[i]);
         Cc.push_back(C_v[i]);
         eCc.push_back(eC_v[i]);
         invC2.push_back(1. / (C_v[i] * C_v[i]));
         e_invC2.push_back(2. * eC_v[i] / (C_v[i] * C_v[i] * C_v[i]));
      }
      int nf = f.size();
      int nv = Vc.size();
      std::vector<double> ef(nf, 0.), eV(nv, 0.);

      // grafici
      TGraphErrors *g1 = new TGraphErrors(nf, f.data(), C_f.data(), ef.data(), eC_f.data());
      g1->SetName(("g_Cf_" + name).c_str());
      g1->SetTitle((name + ": C vs f;f [kHz];C [pF]").c_str());

      TGraphErrors *g2 = new TGraphErrors(nv, Vc.data(), Cc.data(), eV.data(), eCc.data());
      g2->SetName(("g_CV_" + name).c_str());
      g2->SetTitle((name + ": C vs V;|V| [V];C [pF]").c_str());

      TGraphErrors *g3 = new TGraphErrors(nv, Vc.data(), invC2.data(), eV.data(), e_invC2.data());
      g3->SetName(("g_invC2V_" + name).c_str());
      g3->SetTitle((name + ": 1/C^{2} vs V;|V| [V];1/C^{2} [pF^{-2}]").c_str());

      TGraphErrors *gs[3] = { g1, g2, g3 };
      for (int j = 0; j < 3; ++j) {
         gs[j]->SetMarkerStyle(marker);
         gs[j]->SetMarkerSize(0.9);
         gs[j]->SetMarkerColor(color);
         gs[j]->SetLineColor(color);
      }

      // canvas singole
      TCanvas *ci1 = new TCanvas("c_nogain_Cf", "nogain: C vs f", 800, 600);
      ci1->SetGrid(); ci1->SetLogx();
      g1->Draw("AP");
      salva(ci1, "./plots/" + name + "_C_f");

      TCanvas *ci2 = new TCanvas("c_nogain_CV", "nogain: C vs V", 800, 600);
      ci2->SetGrid();
      g2->Draw("AP");
      ci2->Update();
      ci2->Modified(); ci2->Update();
      salva(ci2, "./plots/" + name + "_C_V");

      TCanvas *ci3 = new TCanvas("c_nogain_invC2V", "nogain: 1/C^2 vs V", 800, 600);
      ci3->SetGrid();
      g3->Draw("AP");
      ci3->Update();
      // --- 3 fit lineari su 1/C^2 - V ---
      if (fareFit) {
         TF1 *fA = fitLineare(g3, fitA_min, fitA_max, colFitA, "A");
         TF1 *fB = fitLineare(g3, fitB_min, fitB_max, colFitB, "B");
         TF1 *fC = fitLineare(g3, fitC_min, fitC_max, colFitC, "C");
         stampaNb(fA, A, name, "A");
         stampaNb(fB, A, name, "B");
         stampaNb(fC, A, name, "C");
         ci3->cd();
         fA->Draw("same");
         fB->Draw("same");
         fC->Draw("same");
      }
      ci3->Modified(); ci3->Update();
      salva(ci3, "./plots/" + name + "_invC2_V");

      // --- profilo di drogaggio N_B - W ---
      if (fareProfilo) {
         TGraphErrors *g4 = profiloDrogaggio(Vc, Cc, eCc, A, nSalto, name, color, marker);
         if (g4->GetN() > 0) {
            TCanvas *ci4 = new TCanvas(("c_" + name + "_profilo").c_str(), (name + ": profilo di drogaggio").c_str(), 800, 600);
            ci4->SetGrid(); ci4->SetLogy();   // togli SetLogy() per la scala lineare
            g4->Draw("AP");
            salva(ci4, "./plots/" + name + "_profilo");
         }
      }

      // confronto: nogain compare solo in C - f
      TGraphErrors *m1 = (TGraphErrors*)g1->Clone();
      mg1->Add(m1, "P");
      leg1->AddEntry(m1, name.c_str(), "p");
   }

   // =====================================================================
   // =============================  LGAD  ================================
   // =====================================================================
   if (chip_LGAD) {
      std::string name = "LGAD";
      int color = kRed + 1, marker = 21;
      const double A = A_gain;   // area usata per Nb [cm^2]
      const int nSalto = 7;   // derivata tra il punto i e il punto i+nSalto (aumentalo per ridurre il rumore)

      // range dei fit lineari su 1/C^2 - V [V]  (CAMBIALI TU)
      const double fitA_min = 5,   fitA_max = 20;
      const double fitB_min = 25,  fitB_max = 28;
      const double fitC_min = 30,  fitC_max = 40;

      std::vector<double> f   = { 1, 1.5, 2.3, 3.5, 5.3, 8, 12, 18, 27, 40,
                                  60, 90, 140, 210, 300 };
      std::vector<double> C_f = { 151.5, 150.7, 150.1, 150.2, 150.1, 149.93, 148.5, 148.6, 147.9, 145.4,
                                  138.7, 129.1, 106, 80.7, 55.2 };
      std::vector<double> eC_f= { 0.2, 0.1, 0.1, 0.1, 0.1, 0.05, 0.3, 0.1, 0.1, 0.1,
                                  0.2, 0.2, 1, 0.2, 0.2 };
      std::vector<double> V   = { 0, 2, 4, 6, 8, 10, 12, 14, 16, 18,
                                  19, 20, 20.5, 21, 21.5, 22, 22.5, 23, 23.5, 24,
                                  24.5, 25, 25.5, 26, 27, 28, 29, 30, 31, 32,
                                  34, 36, 38, 40 };
      std::vector<double> C_v = { 307.2, 220.35, 190.11, 172.22, 159.35, 149, 140.29, 132.41, 124.93, 117.36,
                                  113.36, 108.91, 106.51, 103.89, 100.97, 97.56, 93.32, 87.3, 75.92, 45.5,
                                  14.17, 4.81, 3.82, 3.32, 2.99, 2.84, 2.79, 2.76, 2.74, 2.7,
                                  2.72, 2.72, 2.71, 2.7 };
      std::vector<double> eC_v= { 0.1, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                                  0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.3,
                                  0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                                  0.05, 0.05, 0.05, 0.05 };

      // tengo solo i punti con |V| <= Vmax
      std::vector<double> Vc, Cc, eCc, invC2, e_invC2;
      for (size_t i = 0; i < V.size(); ++i) {
         if (V[i] > Vmax) continue;
         Vc.push_back(V[i]);
         Cc.push_back(C_v[i]);
         eCc.push_back(eC_v[i]);
         invC2.push_back(1. / (C_v[i] * C_v[i]));
         e_invC2.push_back(2. * eC_v[i] / (C_v[i] * C_v[i] * C_v[i]));
      }
      int nf = f.size();
      int nv = Vc.size();
      std::vector<double> ef(nf, 0.), eV(nv, 0.);

      // grafici
      TGraphErrors *g1 = new TGraphErrors(nf, f.data(), C_f.data(), ef.data(), eC_f.data());
      g1->SetName(("g_Cf_" + name).c_str());
      g1->SetTitle((name + ": C vs f;f [kHz];C [pF]").c_str());

      TGraphErrors *g2 = new TGraphErrors(nv, Vc.data(), Cc.data(), eV.data(), eCc.data());
      g2->SetName(("g_CV_" + name).c_str());
      g2->SetTitle((name + ": C vs V;|V| [V];C [pF]").c_str());

      TGraphErrors *g3 = new TGraphErrors(nv, Vc.data(), invC2.data(), eV.data(), e_invC2.data());
      g3->SetName(("g_invC2V_" + name).c_str());
      g3->SetTitle((name + ": 1/C^{2} vs V;|V| [V];1/C^{2} [pF^{-2}]").c_str());

      TGraphErrors *gs[3] = { g1, g2, g3 };
      for (int j = 0; j < 3; ++j) {
         gs[j]->SetMarkerStyle(marker);
         gs[j]->SetMarkerSize(0.9);
         gs[j]->SetMarkerColor(color);
         gs[j]->SetLineColor(color);
      }

      // canvas singole
      TCanvas *ci1 = new TCanvas("c_LGAD_Cf", "LGAD: C vs f", 800, 600);
      ci1->SetGrid(); ci1->SetLogx();
      g1->Draw("AP");
      salva(ci1, "./plots/" + name + "_C_f");

      TCanvas *ci2 = new TCanvas("c_LGAD_CV", "LGAD: C vs V", 800, 600);
      ci2->SetGrid();
      g2->Draw("AP");
      ci2->Update();
      g2->GetXaxis()->SetLimits(0., Vmax);
      ci2->Modified(); ci2->Update();
      salva(ci2, "./plots/" + name + "_C_V");

      TCanvas *ci3 = new TCanvas("c_LGAD_invC2V", "LGAD: 1/C^2 vs V", 800, 600);
      ci3->SetGrid();
      g3->Draw("AP");
      ci3->Update();
      g3->GetXaxis()->SetLimits(0., Vmax);
      // --- 3 fit lineari su 1/C^2 - V ---
      if (fareFit) {
         TF1 *fA = fitLineare(g3, fitA_min, fitA_max, colFitA, "A");
         TF1 *fB = fitLineare(g3, fitB_min, fitB_max, colFitB, "B");
         TF1 *fC = fitLineare(g3, fitC_min, fitC_max, colFitC, "C");
         stampaNb(fA, A, name, "A");
         stampaNb(fB, A, name, "B");
         stampaNb(fC, A, name, "C");
         ci3->cd();
         fA->Draw("same");
         fB->Draw("same");
         fC->Draw("same");
      }
      ci3->Modified(); ci3->Update();
      salva(ci3, "./plots/" + name + "_invC2_V");

      // --- profilo di drogaggio N_B - W ---
      if (fareProfilo) {
         TGraphErrors *g4 = profiloDrogaggio(Vc, Cc, eCc, A, nSalto, name, color, marker);
         if (g4->GetN() > 0) {
            TCanvas *ci4 = new TCanvas(("c_" + name + "_profilo").c_str(), (name + ": profilo di drogaggio").c_str(), 800, 600);
            ci4->SetGrid(); ci4->SetLogy();   // togli SetLogy() per la scala lineare
            g4->Draw("AP");
            salva(ci4, "./plots/" + name + "_profilo");
         }
      }

      // confronto: LGAD compare in C - f, C - V e 1/C^2 - V
      TGraphErrors *m1 = (TGraphErrors*)g1->Clone();
      TGraphErrors *m2 = (TGraphErrors*)g2->Clone();
      TGraphErrors *m3 = (TGraphErrors*)g3->Clone();
      mg1->Add(m1, "P");  leg1->AddEntry(m1, name.c_str(), "p");
      mg2->Add(m2, "P");  leg2->AddEntry(m2, name.c_str(), "p");
      mg3->Add(m3, "P");  leg3->AddEntry(m3, name.c_str(), "p");
   }

   // =====================================================================
   // =============================  PiN  =================================
   // =====================================================================
   if (chip_PiN) {
      std::string name = "PiN";
      int color = kBlue + 1, marker = 22;
      const double A = A_gain;   // area usata per Nb [cm^2]
      const int nSalto = 5;   // derivata tra il punto i e il punto i+nSalto (aumentalo per ridurre il rumore)

      // range dei fit lineari su 1/C^2 - V [V]  (CAMBIALI TU)
      const double fitA_min = 0.5,   fitA_max = 2.4;
      const double fitB_min = 5.2,   fitB_max = 40;
      std::vector<double> V   = { 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.2,
                                  1.4, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5,
                                  6, 7, 8, 10, 12, 15, 20, 30, 40 };
      std::vector<double> C_v = { 21.6, 15.20, 9.34, 6.84, 6.07, 5.5, 5.18, 4.58, 4.17, 3.94,
                                  3.75, 3.62, 3.49, 3.38, 3.19, 3.07, 3.00, 2.95, 2.91, 2.89,
                                  2.85, 2.83, 2.81, 2.78, 2.77, 2.75, 2.72, 2.69, 2.67 };
      std::vector<double> eC_v= { 0.1, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                                  0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                                  0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05 };

      // tengo solo i punti con |V| <= Vmax
      std::vector<double> Vc, Cc, eCc, invC2, e_invC2;
      for (size_t i = 0; i < V.size(); ++i) {
         if (V[i] > Vmax) continue;
         Vc.push_back(V[i]);
         Cc.push_back(C_v[i]);
         eCc.push_back(eC_v[i]);
         invC2.push_back(1. / (C_v[i] * C_v[i]));
         e_invC2.push_back(2. * eC_v[i] / (C_v[i] * C_v[i] * C_v[i]));
      }
      int nv = Vc.size();
      std::vector<double> eV(nv, 0.);

      // grafici
      TGraphErrors *g2 = new TGraphErrors(nv, Vc.data(), Cc.data(), eV.data(), eCc.data());
      g2->SetName(("g_CV_" + name).c_str());
      g2->SetTitle((name + ": C vs V;|V| [V];C [pF]").c_str());

      TGraphErrors *g3 = new TGraphErrors(nv, Vc.data(), invC2.data(), eV.data(), e_invC2.data());
      g3->SetName(("g_invC2V_" + name).c_str());
      g3->SetTitle((name + ": 1/C^{2} vs V;|V| [V];1/C^{2} [pF^{-2}]").c_str());

      TGraphErrors *gs[2] = { g2, g3 };
      for (int j = 0; j < 2; ++j) {
         gs[j]->SetMarkerStyle(marker);
         gs[j]->SetMarkerSize(0.9);
         gs[j]->SetMarkerColor(color);
         gs[j]->SetLineColor(color);
      }

      // canvas singole
      TCanvas *ci2 = new TCanvas("c_PiN_CV", "PiN: C vs V", 800, 600);
      ci2->SetGrid();
      g2->Draw("AP");
      ci2->Update();
      g2->GetXaxis()->SetLimits(0., Vmax);
      ci2->Modified(); ci2->Update();
      salva(ci2, "./plots/" + name + "_C_V");

      TCanvas *ci3 = new TCanvas("c_PiN_invC2V", "PiN: 1/C^2 vs V", 800, 600);
      ci3->SetGrid();
      g3->Draw("AP");
      ci3->Update();
      g3->GetXaxis()->SetLimits(0., Vmax);
      // --- 3 fit lineari su 1/C^2 - V ---
      if (fareFit) {
         TF1 *fA = fitLineare(g3, fitA_min, fitA_max, colFitA, "A");
         TF1 *fB = fitLineare(g3, fitB_min, fitB_max, colFitB, "B");
         stampaNb(fA, A, name, "A");
         stampaNb(fB, A, name, "B");
         ci3->cd();
         fA->Draw("same");
         fB->Draw("same");
      }
      ci3->Modified(); ci3->Update();
      salva(ci3, "./plots/" + name + "_invC2_V");

      // --- profilo di drogaggio N_B - W ---
      if (fareProfilo) {
         TGraphErrors *g4 = profiloDrogaggio(Vc, Cc, eCc, A, nSalto, name, color, marker);
         if (g4->GetN() > 0) {
            TCanvas *ci4 = new TCanvas(("c_" + name + "_profilo").c_str(), (name + ": profilo di drogaggio").c_str(), 800, 600);
            ci4->SetGrid(); ci4->SetLogy();   // togli SetLogy() per la scala lineare
            g4->Draw("AP");
            salva(ci4, "./plots/" + name + "_profilo");
         }
      }

      // confronto: PiN compare in C - V e 1/C^2 - V
      TGraphErrors *m2 = (TGraphErrors*)g2->Clone();
      TGraphErrors *m3 = (TGraphErrors*)g3->Clone();
      mg2->Add(m2, "P");  leg2->AddEntry(m2, name.c_str(), "p");
      mg3->Add(m3, "P");  leg3->AddEntry(m3, name.c_str(), "p");
   }

   // =====================================================================
   // ===========================  EX-FLUE  ===============================
   // =====================================================================
   if (chip_EXFLUE) {
      std::string name = "EX-FLUE";
      int color = kGreen + 2, marker = 23;
      const double A = A_gain;   // area usata per Nb [cm^2]
      const int nSalto = 6;   // derivata tra il punto i e il punto i+nSalto (aumentalo per ridurre il rumore)

      // range dei fit lineari su 1/C^2 - V [V]  (CAMBIALI TU)
      const double fitA_min = 5,   fitA_max = 24.5;
      const double fitB_min = 26,  fitB_max = 31;
      const double fitC_min = 32.8,  fitC_max = 39;

      std::vector<double> V   = { 2, 4, 6, 8, 10, 12, 14, 16, 18,
                                  20, 22, 24, 25, 26, 26.5, 27, 27.5, 28, 28.5,
                                  29, 29.5, 30, 30.5, 31, 31.5, 32, 32.5, 33, 34,
                                  35, 36, 37, 38, 39 };
      std::vector<double> C_v = { 19.71, 19.22, 18.86, 18.68, 18.53, 18.41, 18.3, 18.2, 18.1,
                                  18.01, 17.93, 17.62, 17.46, 16.38, 13.95, 12.25, 11.0, 10.05, 9.3,
                                  8.67, 8.16, 7.72, 7.36, 7.06, 6.9, 6.75, 6.66, 6.61, 6.59,
                                  6.57, 6.5, 6.45, 6.40, 6.39 };
      std::vector<double> eC_v= { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                                  0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                                  0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                                  0.05, 0.1, 0.1, 0.1, 0.1 };

      // tengo solo i punti con |V| <= Vmax
      std::vector<double> Vc, Cc, eCc, invC2, e_invC2;
      for (size_t i = 0; i < V.size(); ++i) {
         if (V[i] > Vmax) continue;
         Vc.push_back(V[i]);
         Cc.push_back(C_v[i]);
         eCc.push_back(eC_v[i]);
         invC2.push_back(1. / (C_v[i] * C_v[i]));
         e_invC2.push_back(2. * eC_v[i] / (C_v[i] * C_v[i] * C_v[i]));
      }
      int nv = Vc.size();
      std::vector<double> eV(nv, 0.);

      // grafici
      TGraphErrors *g2 = new TGraphErrors(nv, Vc.data(), Cc.data(), eV.data(), eCc.data());
      g2->SetName(("g_CV_" + name).c_str());
      g2->SetTitle((name + ": C vs V;|V| [V];C [pF]").c_str());

      TGraphErrors *g3 = new TGraphErrors(nv, Vc.data(), invC2.data(), eV.data(), e_invC2.data());
      g3->SetName(("g_invC2V_" + name).c_str());
      g3->SetTitle((name + ": 1/C^{2} vs V;|V| [V];1/C^{2} [pF^{-2}]").c_str());

      TGraphErrors *gs[2] = { g2, g3 };
      for (int j = 0; j < 2; ++j) {
         gs[j]->SetMarkerStyle(marker);
         gs[j]->SetMarkerSize(0.9);
         gs[j]->SetMarkerColor(color);
         gs[j]->SetLineColor(color);
      }

      // canvas singole
      TCanvas *ci2 = new TCanvas("c_EXFLUE_CV", "EX-FLUE: C vs V", 800, 600);
      ci2->SetGrid();
      g2->Draw("AP");
      ci2->Update();
      g2->GetXaxis()->SetLimits(0., Vmax);
      ci2->Modified(); ci2->Update();
      salva(ci2, "./plots/" + name + "_C_V");

      TCanvas *ci3 = new TCanvas("c_EXFLUE_invC2V", "EX-FLUE: 1/C^2 vs V", 800, 600);
      ci3->SetGrid();
      g3->Draw("AP");
      ci3->Update();
      g3->GetXaxis()->SetLimits(0., Vmax);
      // --- 3 fit lineari su 1/C^2 - V ---
      if (fareFit) {
         TF1 *fA = fitLineare(g3, fitA_min, fitA_max, colFitA, "A");
         TF1 *fB = fitLineare(g3, fitB_min, fitB_max, colFitB, "B");
         TF1 *fC = fitLineare(g3, fitC_min, fitC_max, colFitC, "C");
         stampaNb(fA, A, name, "A");
         stampaNb(fB, A, name, "B");
         stampaNb(fC, A, name, "C");
         ci3->cd();
         fA->Draw("same");
         fB->Draw("same");
         fC->Draw("same");
      }
      ci3->Modified(); ci3->Update();
      salva(ci3, "./plots/" + name + "_invC2_V");

      // --- profilo di drogaggio N_B - W ---
      if (fareProfilo) {
         TGraphErrors *g4 = profiloDrogaggio(Vc, Cc, eCc, A, nSalto, name, color, marker);
         if (g4->GetN() > 0) {
            TCanvas *ci4 = new TCanvas(("c_" + name + "_profilo").c_str(), (name + ": profilo di drogaggio").c_str(), 800, 600);
            ci4->SetGrid(); ci4->SetLogy();   // togli SetLogy() per la scala lineare
            g4->Draw("AP");
            salva(ci4, "./plots/" + name + "_profilo");
         }
      }

      // confronto: EX-FLUE compare in C - V e 1/C^2 - V
      TGraphErrors *m2 = (TGraphErrors*)g2->Clone();
      TGraphErrors *m3 = (TGraphErrors*)g3->Clone();
      mg2->Add(m2, "P");  leg2->AddEntry(m2, name.c_str(), "p");
      mg3->Add(m3, "P");  leg3->AddEntry(m3, name.c_str(), "p");
   }

   // =====================================================================
   // =========================  CANVAS DI CONFRONTO  =====================
   // =====================================================================
   if (fareConfronto) {

      // C - f (si disegna solo se c'e' almeno un chip che la contiene)
      if (mg1->GetListOfGraphs() && mg1->GetListOfGraphs()->GetSize() > 0) {
         TCanvas *c1 = new TCanvas("c_conf_Cf", "Confronto C vs f", 800, 600);
         c1->SetGrid(); c1->SetLogx();
         mg1->Draw("A");
         leg1->Draw();
         salva(c1, "./plots/confronto_C_f");
      }

      // C - V
      if (mg2->GetListOfGraphs() && mg2->GetListOfGraphs()->GetSize() > 0) {
         TCanvas *c2 = new TCanvas("c_conf_CV", "Confronto C vs V", 800, 600);
         c2->SetGrid(); c2->SetLogy();   // togli SetLogy() per la scala lineare
         mg2->Draw("A");
         c2->Update();
         mg2->GetXaxis()->SetLimits(0., Vmax);
         leg2->Draw();
         c2->Modified(); c2->Update();
         salva(c2, "./plots/confronto_C_V");
      }

      // 1/C^2 - V
      if (mg3->GetListOfGraphs() && mg3->GetListOfGraphs()->GetSize() > 0) {
         TCanvas *c3 = new TCanvas("c_conf_invC2V", "Confronto 1/C^2 vs V", 800, 600);
         c3->SetGrid(); c3->SetLogy();   // togli SetLogy() per la scala lineare
         mg3->Draw("A");
         c3->Update();
         mg3->GetXaxis()->SetLimits(0., Vmax);
         leg3->Draw();
         c3->Modified(); c3->Update();
         salva(c3, "./plots/confronto_invC2_V");
      }
   }
}