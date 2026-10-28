#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TCanvas.h>
#include <TH1D.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TF2.h>
#include <TLegend.h>
#include <TStyle.h>

void p2p_time_correlation(){
    gStyle->SetOptStat(0);

    TFile *inputfile = TFile::Open("data.root");

    TCanvas *c1 = new TCanvas("c1", "p2p_corr", 1600, 1200);
    c1->Divide(2, 2);

    TH1D *histo_signal[4];
    TTree *inputTree[4];
    
    const int p[4] = {5, 6, 9, 10};

    double signal[4];
    double t_10[4];
    double t_50[4];
    
    double minBin = 1.0;
    double maxBin = 91.0;

    int RebinFactor = 50;
    double scope_res = 0.01;
    int Nbins = (maxBin - minBin) / scope_res;

    for(int i=0; i<4; i++)
    {
        histo_signal[i] = new TH1D(Form("Pixel %d", i), Form("Pixel %d", p[i]), Nbins, minBin-0.005, maxBin+0.005);
        histo_signal[i]->Rebin(RebinFactor);
        histo_signal[i]->GetXaxis()->SetTitle("V [mV]");
        histo_signal[i]->GetYaxis()->SetTitle("Counts");
        histo_signal[i]->SetLineColor(kBlack);
        histo_signal[i]->SetLineWidth(2);
        histo_signal[i]->SetDirectory(0);

        inputTree[i] = (TTree*)inputfile->Get(Form("pixel%d", p[i]));
        inputTree[i]->SetBranchAddress("signal", &signal[i]);
    }
}