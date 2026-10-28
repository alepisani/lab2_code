#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TCanvas.h>
#include <TH1D.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TF1.h>
#include <TLegend.h>
#include <TStyle.h>

using namespace std;

void energy_calib(int RebinFactor=46)
{
    TFile *inputfile = TFile::Open("data.root");
    TCanvas *c1 = new TCanvas("c1", "Scope Signals", 1600, 1200);
    c1->Divide(2, 2);
    TCanvas *c2 = new TCanvas("c2", "Energy Calibration", 1600, 1200);
    c2->Divide(2, 2);

    TH1D *histo_signal[4];
    TTree *inputTree[4];

    double minBin = 1.0;
    double maxBin = 91.0;

    const int p[4] = {5,6,9,10};
    double signal[4];
    

    double scope_res = 0.01;
    int Nbins = (maxBin - minBin) / scope_res;
    gStyle->SetOptStat(0);
    for(int i=0; i<4; i++)
    {
        histo_signal[i] = new TH1D(Form("Pixel %d", i), Form("Pixel %d", p[i]), Nbins, minBin-0.005, maxBin+0.005);
        histo_signal[i]->Rebin(RebinFactor);
        histo_signal[i]->GetXaxis()->SetTitle("Signal [mV]");
        histo_signal[i]->GetYaxis()->SetTitle("Counts");
        histo_signal[i]->SetLineColor(kBlack);
        histo_signal[i]->SetLineWidth(2);
        histo_signal[i]->SetDirectory(0);

        inputTree[i] = (TTree*)inputfile->Get(Form("pixel%d", p[i]));
        inputTree[i]->SetBranchAddress("signal", &signal[i]);
    }
    
    int Nevents = inputTree[0]->GetEntries();

    for(int i=0; i<Nevents; i++)
    {
        for(int j=0; j<4; j++)
        {
            inputTree[j]->GetEntry(i);
            histo_signal[j]->Fill(signal[j]);
        }
    }

    TH1D *histo_signal_clone[4];
    for(int i=0; i<4; i++)
    {
        c1->cd(i+1);
        histo_signal_clone[i] = (TH1D*)histo_signal[i]->Clone(Form("Pixel %d clone", i));
        histo_signal_clone[i]->SetDirectory(0);
        histo_signal_clone[i]->Draw();
    }

    double Norm_ka[4] = {106,114,112,101};
    double Norm_kb[4] = {20,9,13,11};

    double mu_ka[4] = {74,74,76,77};   //muuuu
    double mu_kb[4] = {80.5,81,83,84};

    double sigma_ka[4] = {0.9,0.9,0.9,0.9};
    double sigma_kb[4] = {1.9,0.9,0.9,0.9};

    double xmin_ka[4] = {72.25,72,75,75.5};
    double xmin_kb[4] = {77,78,81,81};

    double xmax_ka[4] = {78,77,79,80};
    double xmax_kb[4] = {83,83,86,88};


    double x_calib[2] = {5.895, 6.490};
    double y_calib[4][2];
    double y_calib_err[4][2];

    for(int i=0; i<4; i++)
    {
        c2->cd(i+1);
        histo_signal[i]->GetXaxis()->SetRangeUser(60,90);
        histo_signal[i]->Draw();
        TF1 *fit_ka = new TF1(Form("fit_ka_%d", i), "gaus", xmin_ka[i], xmax_ka[i]);
        fit_ka->SetParameters(Norm_ka[i], mu_ka[i], sigma_ka[i]);
        fit_ka->SetLineColor(kRed);
        fit_ka->SetLineWidth(2);
        histo_signal[i]->Fit(fit_ka, "RQ+");

        TF1 *fit_kb = new TF1(Form("fit_kb_%d", i), "gaus", xmin_kb[i], xmax_kb[i]);
        fit_kb->SetParameters(Norm_kb[i], mu_kb[i], sigma_kb[i]);
        fit_kb->SetLineColor(kBlue);
        fit_kb->SetLineWidth(2);
        histo_signal[i]->Fit(fit_kb, "RQ+");

        y_calib[i][0] = fit_ka->GetParameter(1);
        y_calib[i][1] = fit_kb->GetParameter(1);
        y_calib_err[i][0] = fit_ka->GetParameter(2);
        y_calib_err[i][1] = fit_kb->GetParameter(2);

        cout << "Pixel " << p[i] << ": " << endl;
        cout << "  - Ka: " << fit_ka->GetParameter(1) << " +/- " << fit_ka->GetParameter(2) << endl;
        cout << "        norm: " << fit_ka->GetParameter(0) << endl;
        cout << "        p-value: " << fit_ka->GetProb() << endl;
        cout << "  - Kb: " << fit_kb->GetParameter(1) << " +/- " << fit_kb->GetParameter(2) << endl;
        cout << "        norm: " << fit_kb->GetParameter(0) << endl;
        cout << "        p-value: " << fit_kb->GetProb() << endl;

        TLegend *leg = new TLegend(0.12, 0.73, 0.37, 0.88);
        leg->AddEntry(fit_ka, "K_{#alpha}", "l");
        leg->AddEntry(fit_kb, "K_{#beta}", "l");
        leg->SetBorderSize(0);
        leg->SetTextSize(0.05);
        leg->Draw();
        
    }

    cout << "\n\n================================ Calibration ================================" << endl;
    TCanvas *c3 = new TCanvas("c3", "Energy Calibration", 1600, 1200);
    TMultiGraph *mg = new TMultiGraph();
    TGraphErrors *gr[4];
    TF1 *fit_calib[4];

    TLegend *leg_calib = new TLegend(0.12, 0.73, 0.37, 0.88);
    leg_calib->SetBorderSize(0);

    int color[4] = {kRed, kBlue, kGreen+2, kOrange};
    for(int i=0; i<4; i++)
    {
        gr[i] = new TGraphErrors(2, x_calib, y_calib[i], nullptr, y_calib_err[i]);
        gr[i]->SetMarkerColor(color[i]);
        gr[i]->SetLineColor(color[i]);
        gr[i]->SetMarkerStyle(20);

        fit_calib[i] = new TF1(Form("fit_calib_%d", i), "[0]*x", -0.2, 6.6);
        fit_calib[i]->SetLineColor(color[i]);
        gr[i]->Fit(fit_calib[i], "RQ0");

        cout << "Pixel " << p[i] << ": " << endl;
        cout << "  - Calibration: " << fit_calib[i]->GetParameter(0) << " +/- " << fit_calib[i]->GetParError(0) << endl;
        cout << "        p-value: " << fit_calib[i]->GetProb() << endl;

        leg_calib->AddEntry(gr[i], Form("Pixel %d", p[i]), "ple");

        mg->Add(gr[i]);
    }

    TGraph *gr_origin = new TGraph(1);
    gr_origin->SetPoint(0, 0, 0);
    gr_origin->SetMarkerStyle(20);
    gr_origin->SetMarkerSize(1.2);
    gr_origin->SetMarkerColor(kBlack);
    
    mg->SetTitle("Energy Calibration;Energy [keV];Signal [mV]");
    mg->Draw("AP");

    leg_calib->Draw();

    mg->GetXaxis()->SetLimits(-0.5, 6.7);
    mg->SetMinimum(-5);
    mg->SetMaximum(90);

    for(int i=0; i<4; i++)
        fit_calib[i]->Draw("same");

    gr_origin->Draw("P same");

    c3->Modified();
    c3->Update();

    TCanvas *c4 = new TCanvas("c4", "Energy Calibration Zoom", 1600, 1200);
    TMultiGraph *mg_zoom = (TMultiGraph*)mg->Clone("mg_zoom");
    mg_zoom->SetTitle("Energy Calibration Zoom;Energy [keV];Signal [mV]");
    mg_zoom->Draw("AP");
    for(int i=0; i<4; i++)
        fit_calib[i]->Draw("same");   
    mg_zoom->GetXaxis()->SetLimits(5.8, 6.6);
    mg_zoom->SetMinimum(71);
    mg_zoom->SetMaximum(87);
    leg_calib->Draw();

// ---------- Istogrammi in keV ----------
    const int    NbinsE = Nbins;
    const double Emin = (minBin-0.005)/fit_calib[0]->GetParameter(0);
    const double Emax = (maxBin+0.005)/fit_calib[0]->GetParameter(0);

    double slope[4];
    TH1D *histo_energy[4];
    TH1D *histo_energy_all = new TH1D("h_energy_all", "Energy Spectru;Energy [keV];Counts", NbinsE, Emin, Emax);
    histo_energy_all->SetDirectory(0);
    histo_energy_all->SetLineColor(kBlack);
    histo_energy_all->SetLineWidth(2);
    histo_energy_all->Rebin(RebinFactor);

    for(int i=0; i<4; i++)
    {
        slope[i] = fit_calib[i]->GetParameter(0);

        histo_energy[i] = new TH1D(Form("h_energy_%d", i), Form("Pixel %d;Energy [keV];Counts", p[i]), NbinsE, Emin, Emax);
        histo_energy[i]->SetDirectory(0);
        histo_energy[i]->SetLineColor(kBlack);
        histo_energy[i]->SetLineWidth(2);
        histo_energy[i]->Rebin(RebinFactor);
    }

    for(int k=0; k<Nevents; k++)
    {
        for(int j=0; j<4; j++)
        {
            inputTree[j]->GetEntry(k);
            double E = signal[j] / slope[j];
            histo_energy[j]->Fill(E);
            histo_energy_all->Fill(E);
        }
    }

    // Canvas con i 4 pixel separati
    TCanvas *c5 = new TCanvas("c5", "Energy Spectra", 1600, 1200);
    c5->Divide(2, 2);
    for(int i=0; i<4; i++)
    {
        c5->cd(i+1);
        histo_energy[i]->Draw();
    }

    // Canvas con tutti i pixel sommati
    TCanvas *c6 = new TCanvas("c6", "Energy Spectrum - All Pixels", 1600, 1200);
    histo_energy_all->Draw();

    inputfile->Close();
}