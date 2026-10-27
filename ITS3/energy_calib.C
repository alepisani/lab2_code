#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TCanvas.h>
#include <TH1D.h>
#include <TStyle.h>

using namespace std;

void energy_calib()
{
    TFile *inputfile = TFile::Open("GRUPPO1_55Fe/data.root");
    TTree *input_scopeTree  = (TTree*)inputfile->Get("Scope data");
    TTree *input_ADCTree = (TTree*)inputfile->Get("ADC data");

    vector<double> baseline_scope(4*3);
    vector<double> underline_scope(4*3);
    vector<double> signal_scope(4);
    vector<double> t_scope(4*10);
    input_scopeTree->SetBranchAddress("baseline", baseline_scope.data());
    input_scopeTree->SetBranchAddress("underline", underline_scope.data());
    input_scopeTree->SetBranchAddress("signal", signal_scope.data());
    input_scopeTree->SetBranchAddress("t", t_scope.data());

    const int Neventi = input_scopeTree->GetEntries();

    TCanvas *c1 = new TCanvas("c1", "Signal Histograms", 1600, 1200);
    c1->cd();
    c1->Divide(2, 2);
    TH1D *histograms[4];

    int pixels[4] = {5, 6, 9, 10}; // Pixels to analyze

    // Capire come scegliere il binning
    for(int i = 0; i < 4; i++) 
        histograms[i] = new TH1D(Form("histogram_%d", pixels[i]), Form("Segnale pixel %d", pixels[i]), 256, 0.1, 100);

    for (int i = 0; i < Neventi; i++) 
    {
        input_scopeTree->GetEntry(i);
        for (int j = 0; j < 4; j++) 
            histograms[j]->Fill(signal_scope[j]);
    }


    for(int i = 0; i < 4; i++) 
    {
        c1->cd(i + 1);
        histograms[i]->Draw();
    }

}