#include <TGraphErrors.h>
#include <iostream>
#include <TCanvas.h>
#include <TLine.h>

void baseline_vbias()
{
    double Vbias[]    = {0.001, 1.196, 2.399, 3.599, 4.797};
    double err_bias[] = {0.01,  0.01,  0.01,  0.01,  0.01};

    double baseline[] = {346.6, 281.6, 246.2, 227.4, 214.0};
    double err_base[] = {5.0,   5.0,   5.0,   5.0,   5.0};
    double signal[]   = {26.2,  41.8,  55.8,  61.2,  67.0};
    double err_sig[]  = {1.0,   1.0,   1.0,   1.0,   1.0};

    // Nota: errori su bias e signal discussi in laboratorio, dipendono da utilizzo dei cursori e oscillazioni di Vbias rispetto al valore nominale
    // Errore su baseline messo 5, vediamo poi se sistemare con RMS


    TCanvas *c1 = new TCanvas("c1", "Baseline vs Vbias", 800, 600);
    TGraphErrors *gr = new TGraphErrors(5, Vbias, baseline, err_bias, err_base);
    gr->SetTitle("Baseline vs Vbias;Vbias [V];Baseline [mV]");
    gr->SetMarkerStyle(21);
    gr->SetMarkerColor(kBlue);
    gr->Draw("AP");
    c1->Update();
    TLine *line = new TLine(4.8, gPad->GetUymin(), 4.8, gPad->GetUymax());
    line->SetLineColorAlpha(kGray+1, 0.8);
    line->SetLineWidth(2);
    line->SetLineStyle(2);
    line->Draw("same");

    TCanvas *c2 = new TCanvas("c2", "Signal vs Vbias", 800, 600);
    TGraphErrors *gr2 = new TGraphErrors(5, Vbias, signal, err_bias, err_sig);
    gr2->SetTitle("Signal vs Vbias #left(V_{H} = 1.2 V#right);Vbias [V];Signal [mV]");
    gr2->SetMarkerStyle(21);
    gr2->SetMarkerColor(kRed);
    gr2->Draw("AP");
    c2->Update();
    TLine *line2 = new TLine(4.8, gPad->GetUymin(), 4.8, gPad->GetUymax());
    line2->SetLineColorAlpha(kGray+1, 0.8);
    line2->SetLineWidth(2);
    line2->SetLineStyle(2);
    line2->Draw("same");
}