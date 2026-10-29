#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TCanvas.h>
#include <TH2D.h>
#include <TStyle.h>

void p2p_time_correlation(){
    gStyle->SetOptStat(0);

    TFile *inputfile = TFile::Open("data.root");
    if (!inputfile || inputfile->IsZombie()) {
        std::cerr << "Errore nell'apertura di data.root" << std::endl;
        return;
    }

    TCanvas *c1 = new TCanvas("c1", "p2p_corr_t_50", 1600, 1200);
    c1->Divide(2, 2);

    TCanvas *c2 = new TCanvas("c2", "p2p_corr_t_90", 1600, 1200);
    c2->Divide(2, 2);

    TH2D *histo_t50[4];
    TH2D *histo_t90[4];
    TTree *inputTree[4];

    const int p[4] = {5, 6, 9, 10};

    double signal[4];
    double t_10[4];
    double t_50[4];
    double t_90[4];
    
    double minBinV = 1.0; // Taglio eventi sotto i 1 mV (rumore/baseline)
    double maxBinV = 90.0;
    int NbinsV = 100;

    int Nbinst = 100;

    for(int i = 0; i < 4; i++) {
        // Istogrammi per t50 - t10 (valori positivi: ~0.3 ns)
        histo_t50[i] = new TH2D(
            Form("h_pixel%d_t50", p[i]), 
            Form("Pixel %d: t_{50} - t_{10};Signal [mV];t_{50} - t_{10} [ns]", p[i]), 
            NbinsV, minBinV, maxBinV, 
            Nbinst, 0.0, 2.0
        );
        histo_t50[i]->SetDirectory(0);

        // Istogrammi per t90 - t10 (valori negativi: ~ -0.07 ns)
        histo_t90[i] = new TH2D(
            Form("h_pixel%d_t90", p[i]), 
            Form("Pixel %d: t_{90} - t_{10};Signal [mV];t_{90} - t_{10} [ns]", p[i]), 
            NbinsV, minBinV, maxBinV, 
            Nbinst, -1.0, 1.0
        );
        histo_t90[i]->SetDirectory(0);

        // Connessione ai Tree
        inputTree[i] = (TTree*)inputfile->Get(Form("pixel%d", p[i]));
        inputTree[i]->SetBranchAddress("signal", &signal[i]);
        inputTree[i]->SetBranchAddress("CFD_times.t10", &t_10[i]);
        inputTree[i]->SetBranchAddress("CFD_times.t50", &t_50[i]);
        inputTree[i]->SetBranchAddress("CFD_times.t90", &t_90[i]);
    }

    int Nevents = inputTree[0]->GetEntries();

    // UNICO loop sugli eventi e sui 4 pixel
    for(int i = 0; i < Nevents; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            inputTree[j]->GetEntry(i); // Carica tutti i rami del pixel per l'evento i

            if (signal[j] > minBinV) {
                // Riempiamo il grafico con t50 - t10
                if (t_10[j] > -900 && t_50[j] > -900) {
                    histo_t50[j]->Fill(signal[j], -t_10[j] + t_50[j]);
                }
                // Riempiamo il grafico con t90 - t10
                if (t_10[j] > -900 && t_90[j] > -900) {
                    histo_t90[j]->Fill(signal[j], -t_10[j] + t_90[j]);
                }
            }
        }
    }

    // Disegno Canvas 1 (t50)
    for(int i = 0; i < 4; i++)
    {
        c1->cd(i + 1);
        gPad->SetRightMargin(0.13);
        histo_t50[i]->Draw("COLZ");
    }

    // Disegno Canvas 2 (t90)
    for(int i = 0; i < 4; i++)
    {
        c2->cd(i + 1);
        gPad->SetRightMargin(0.13);
        histo_t90[i]->Draw("COLZ");
    }

    inputfile->Close();
}