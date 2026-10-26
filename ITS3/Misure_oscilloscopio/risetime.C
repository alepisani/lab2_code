#include <iostream>
#include <TCanvas.h>
#include <TGraph.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TLine.h>
#include <TF1.h>
#include <fstream>

using namespace std;

#define OUT false

void fit(TGraph *g, double tmin, double tmax, double Vrms, bool print, double& slope, double& intercept);
void normalize(TGraph *g, double t50)
{
    int n = g->GetN();
    double *x = g->GetX();
    double *y = g->GetY();
    for(int i=0;i<n;i++)
    {
        x[i] -= t50;
    }
}

void risetime()
{
    double tmin = -25.;
    double tmax = 25.;
    double t_baseline = 5.;

    const int n_signals = 4;

    TFile *inputfile = new TFile("segnali_clean.root", "READ");

    TGraph *g1 = (TGraph*)inputfile->Get("Segnale 1");
    TGraph *g2 = (TGraph*)inputfile->Get("Segnale 2");
    TGraph *g3 = (TGraph*)inputfile->Get("Segnale 3");
    TGraph *g4 = (TGraph*)inputfile->Get("Segnale 4");
    //TGraph *g5 = (TGraph*)inputfile->Get("Segnale 5");

    double baseline[n_signals], signal[n_signals], Vrms[n_signals];
    
    for (int i = 0; i < n_signals; i++)
    {
        TGraph *g = (TGraph*)inputfile->Get(Form("Segnale %d", i+1));
        int n = g->GetN();
        double *x = g->GetX();
        double *y = g->GetY();
        double sum = 0;
        int count = 0;
        do
        {
            sum += y[count];
            count++;
        }while(x[count] < tmin + t_baseline && count < n);
        baseline[i] = sum / count;

        sum = 0;
        for(int j=0;j<count;j++)
        {
            sum += (y[j] - baseline[i]) * (y[j] - baseline[i]);
        }
        Vrms[i] = sqrt(sum / count);

        count = 0;
        sum = 0;
        do
        {
            sum += y[n-count-1];
            count++;
        }while(x[n-count-1] > tmax - t_baseline && count < n);
        signal[i] = sum / count;

        cout << "Ampiezza segnale " << i+1 << " (mV): " << signal[i] - baseline[i] << endl;
    }
    inputfile->Close();

    double dVdt[n_signals], intercept[n_signals];
    
    //Fidatevi, fittano
    fit(g1, -0.07, 0.14, Vrms[0], false, dVdt[0], intercept[0]);
    fit(g2, -0.5, 0, Vrms[1], false, dVdt[1], intercept[1]);
    fit(g3, -0.8, 0, Vrms[2], false, dVdt[2], intercept[2]);
    fit(g4, -0.5, 0.4, Vrms[3], false, dVdt[3], intercept[3]);

    double amplitudes[n_signals], amplitudes_err[n_signals], signal10[n_signals], signal50[n_signals];
    for(int i=0;i<n_signals;i++)
    {
        amplitudes[i] = signal[i] - baseline[i];
        signal10[i] = baseline[i] + 0.1 * amplitudes[i];
        signal50[i] = baseline[i] + 0.5 * amplitudes[i];
        amplitudes_err[i] = Vrms[i]*sqrt(2);
        amplitudes[i] = fabs(amplitudes[i]);
    }

    double time10[n_signals], time50[n_signals], risetime[n_signals], risetime_err[n_signals];
    //calcolo tempi a partire dalla retta
    for(int i=0;i<n_signals;i++)
    {
        time10[i] = (signal10[i] - intercept[i]) / dVdt[i];
        time50[i] = (signal50[i] - intercept[i]) / dVdt[i];

        risetime[i] = time50[i] - time10[i];
        risetime_err[i] = sqrt(2) * Vrms[i] / fabs(dVdt[i]);

        cout << "Risetime segnale " << i+1 << " (ns): " << risetime[i] << " +/- " << risetime_err[i] << endl;
    }

    TCanvas *c_signals = new TCanvas("c_signals", "Segnali", 800, 600);
    TMultiGraph *mg = new TMultiGraph();
    
    normalize(g1,time50[0]);
    normalize(g2,time50[1]);
    normalize(g3,time50[2]);
    normalize(g4,time50[3]);

    mg->Add(g1, "L");
    mg->Add(g2, "L");
    mg->Add(g3, "L");
    mg->Add(g4, "L");
    mg->SetTitle("Segnali all'oscilloscopio;Tempo (ns);Ampiezza (mV)");
    mg->Draw("A");
    
    TLegend *leg = new TLegend(0.12, 0.12, 0.32, 0.32);
    leg->AddEntry(g1, "Segnale 1", "l");
    leg->AddEntry(g2, "Segnale 2", "l");
    leg->AddEntry(g3, "Segnale 3", "l");
    leg->AddEntry(g4, "Segnale 4", "l");
    leg->Draw();

    TCanvas *c_risetime = new TCanvas("c_risetime", "Risetime vs Ampiezza", 800, 600);
    TGraphErrors *g_risetime = new TGraphErrors(n_signals, amplitudes, risetime, amplitudes_err, risetime_err);
    g_risetime->SetTitle("Risetime vs Ampiezza;Ampiezza (mV);Risetime (ns)");
    g_risetime->SetMarkerStyle(20);
    g_risetime->SetMarkerColor(kBlue);
    g_risetime->Draw("AP");

    c_signals->SaveAs("risultati/segnali.tex");
    c_signals->SaveAs("risultati/segnali.png");

    c_risetime->SaveAs("risultati/risetime.tex");
    c_risetime->SaveAs("risultati/risetime.png");

    if(OUT)
    {
        ofstream outputfile("risultati/risetime.txt");
        outputfile << "Segnale \tBaseline (mV) \tSegnale (mV) \tVrms (mV) \tSlope (mV/ns) \tT10 (ns) \tT50 (ns)" << endl;
        for(int i=0;i<n_signals;i++)
            outputfile << i+1 << "\t" << baseline[i] << "\t" << signal[i] << "\t" << Vrms[i] << "\t" << dVdt[i] << "\t" << time10[i] << "\t" << time50[i] << endl;
        outputfile.close();
    }

}

void fit(TGraph *g, double tmin, double tmax, double Vrms, bool print, double& slope, double& intercept)
{
    int n = g->GetN();
    double *x = g->GetX();
    double *y = g->GetY();
    double erry[n];
    for(int i=0;i<n;i++)
    {
        erry[i] = Vrms;
    }
    TGraphErrors *g_err = new TGraphErrors(n, x, y, nullptr, erry);
    TF1 *f = new TF1("f", "[0]*x + [1]", tmin, tmax);
    g_err->Fit(f, "RQ");
    slope = f->GetParameter(0);
    intercept = f->GetParameter(1);
    if(print)
    {
        TCanvas *c = new TCanvas("c", "c", 800, 600);
        g_err->Draw("AP");
        f->Draw("same");
        cout << "Slope: " << slope << endl;
        cout << "p-value: " << f->GetProb() << endl;
        cout << "DOF: " << f->GetNDF() << endl;
    }
}