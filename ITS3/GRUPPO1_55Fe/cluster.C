#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <iostream>
#include <vector>

using namespace std;

void cluster() {
    
    //------------------------------------------------------------------------Mappa Pixel 2x2------------------------

    TFile *inputfile = TFile::Open("data.root", "UPDATE");

    // I tuoi 4 pixel fisici mappati su una griglia 2x2
    // p[0] -> indice 0 (es. x=0, y=0)
    // p[1] -> indice 1 (es. x=1, y=0)
    // p[2] -> indice 2 (es. x=0, y=1)
    // p[3] -> indice 3 (es. x=1, y=1)
    
    const int p[4] = {5, 6, 9, 10};

    TCanvas *c1 = new TCanvas("c1", "Mappa Pixel 2x2", 800, 600);    
    TH2D *h_map = new TH2D("h_map", "Stato Pixel (2x2);Asse X;Asse Y", 2, 0, 2, 2, 0, 2);
    gStyle->SetOptStat(0);
    
    std::vector<TTree*> input_pixels(4);
    double baseline[4] = {-999., -999., -999., -999.};
    double signal[4] = {};

    for (int i = 0; i < 4; i++) {
        TString treename = TString::Format("pixel%d", p[i]);
        input_pixels[i] = (TTree*)inputfile->Get(treename);

        if (input_pixels[i] && input_pixels[i]->GetBranch("baseline")) {
            input_pixels[i]->SetBranchAddress("baseline", &baseline[i]);
        } else {
            std::cerr << "Attenzione: impossibile trovare il tree o il branch baseline di "
                      << treename << std::endl;
            input_pixels[i] = nullptr;
        }
        if (input_pixels[i] && input_pixels[i]->GetBranch("signal")) {
            input_pixels[i]->SetBranchAddress("signal", &signal[i]);
        } else {
            std::cerr << "Attenzione: impossibile trovare il branch signal di "
                      << treename << std::endl;
            input_pixels[i] = nullptr;
        }
    }


    Long64_t nEntries = input_pixels[0]->GetEntries();
    std::cout << "Numero totale di eventi da analizzare: " << nEntries << std::endl;

    // 4. Ciclo sugli eventi 
    for (Long64_t ev = 0; ev < nEntries; ev++) {
        
        int pixel_accesi_in_questo_evento = 0;

        for (int i = 0; i < 4; i++) {
            if (!input_pixels[i]) continue;

            input_pixels[i]->GetEntry(ev);

            // Condizione: se la baseline è diversa da -999, il pixel ha registrato un evento
            if (baseline[i] != -999.) {
                int x = i % 2;         // Colonna (0 o 1)
                int y = 1 - (i / 2);   // Riga invertita: 0 e 1 vanno in alto (y=1), 2 e 3 in basso (y=0)
                
                h_map->Fill(x + 0.5, y + 0.5, 1);
                
                pixel_accesi_in_questo_evento++;
            }
        }
    }
    
    // 5. Disegno della mappa finale
    h_map->GetXaxis()->SetNdivisions(2);
    h_map->GetYaxis()->SetNdivisions(2);
    h_map->SetMinimum(0);
    
    c1->cd();
    h_map->Draw("COLZ TEXT");
    c1->cd();
    h_map->SetMarkerSize(1.8);        // dimensione del testo con il numero di eventi
    h_map->Draw("COLZ TEXT");

    // etichetta con il numero del pixel in ogni cella
    for (int i = 0; i < 4; i++) {
        int x = i % 2;
        int y = 1 - (i / 2);          // stessa mappatura usata nel Fill

        TLatex *label = new TLatex(x + 0.5, y + 0.7, Form("Pixel %d", p[i]));
        label->SetTextAlign(22);      // centrato orizzontalmente e verticalmente
        label->SetTextSize(0.03);
        label->SetTextColor(kBlack);
        label->Draw();
    }

    c1->Modified();
    c1->Update();


    //---------------------------------------------------------------------------------------------------------------------------------------------

    //aggiungo cluster a data.root
    inputfile->cd();   // il tree viene associato al file
    TTree *t_cluster = new TTree("cluster", "Cluster dei pixel interni");

    Int_t    b_size = 0;                 // numero di pixel accesi (cluster size)
    Double_t b_totsig = 0.;              // segnale totale del cluster
    Int_t    b_on[4] = {0, 0, 0, 0};     // 1 se il pixel p[i] e' acceso
    Double_t b_sig[4] = {0., 0., 0., 0.};// segnale del singolo pixel (0 se spento)

    t_cluster->Branch("cluster_size", &b_size,   "cluster_size/I");
    t_cluster->Branch("total_signal", &b_totsig, "total_signal/D");
    t_cluster->Branch("pixel_on",     b_on,      "pixel_on[4]/I");
    t_cluster->Branch("pixel_signal", b_sig,     "pixel_signal[4]/D");

    // creazione istogramma per contare quanti eventi hanno registrato un cluster e quanti pixel si sono accesi in quell'evento
    TCanvas *c2 = new TCanvas("c2", "Cluster Count", 800, 600);
    TH1D *h_cluster = new TH1D("h_cluster", "Numero di pixel accesi per evento;Numero di pixel accesi;Numero di eventi", 5, -0.5, 4.5);
    h_cluster->GetXaxis()->SetBinLabel(1, "0");
    h_cluster->GetXaxis()->SetBinLabel(2, "1");
    h_cluster->GetXaxis()->SetBinLabel(3, "2");
    h_cluster->GetXaxis()->SetBinLabel(4, "3");
    h_cluster->GetXaxis()->SetBinLabel(5, "4");

    // un unico canvas diviso in 4 pad (2x2) per il segnale totale in funzione del cluster size
    TCanvas *c4 = new TCanvas("c4", "Signal vs Cluster Size (1D)", 1200, 900);
    c4->Divide(2, 2);

    TH1D *h_signal_cluster_1 = new TH1D("h_signal_cluster_1", "Segnale per cluster size = 1;Segnale [mV];Numero di eventi", 100, 0, 100);
    TH1D *h_signal_cluster_2 = new TH1D("h_signal_cluster_2", "Segnale per cluster size = 2;Segnale [mV];Numero di eventi", 100, 0, 100);
    TH1D *h_signal_cluster_3 = new TH1D("h_signal_cluster_3", "Segnale per cluster size = 3;Segnale [mV];Numero di eventi", 100, 0, 100);
    TH1D *h_signal_cluster_4 = new TH1D("h_signal_cluster_4", "Segnale per cluster size = 4;Segnale [mV];Numero di eventi", 100, 0, 100);

    // array comodo: signal_histograms[n-1] e' l'istogramma per cluster size = n
    TH1D *signal_histograms[4] = {
        h_signal_cluster_1, h_signal_cluster_2,
        h_signal_cluster_3, h_signal_cluster_4
    };

    int no_events = 0;
    int one_pixel_on = 0;
    int two_pixels_on = 0;
    int three_pixels_on = 0;
    int four_pixels_on = 0;

    for (Long64_t ev = 0; ev < nEntries; ev++) {
        int pixel_accesi_in_questo_evento = 0;
        double total_signal = 0.;

        for (int i = 0; i < 4; i++) { b_on[i] = 0; b_sig[i] = 0.; }

        for (int i = 0; i < 4; i++) {
            if (!input_pixels[i]) continue;
            input_pixels[i]->GetEntry(ev);

            if (baseline[i] != -999.) {
                pixel_accesi_in_questo_evento++;
                total_signal += signal[i];
                b_on[i]  = 1;
                b_sig[i] = signal[i];
            }
        }

        // riempimento del nuovo tree
        b_size   = pixel_accesi_in_questo_evento;
        b_totsig = total_signal;
        t_cluster->Fill();

        // istogrammi generali dell'evento
        h_cluster->Fill(pixel_accesi_in_questo_evento);

        if (pixel_accesi_in_questo_evento > 0) {
            // segnale totale in funzione del cluster size: UN solo Fill per evento
            signal_histograms[pixel_accesi_in_questo_evento - 1]->Fill(total_signal);
        }

        if (pixel_accesi_in_questo_evento == 0) no_events++;
        if (pixel_accesi_in_questo_evento == 1) one_pixel_on++;
        if (pixel_accesi_in_questo_evento == 2) two_pixels_on++;
        if (pixel_accesi_in_questo_evento == 3) three_pixels_on++;
        if (pixel_accesi_in_questo_evento == 4) four_pixels_on++;
    }   // fine ciclo sugli eventi

    inputfile->cd();
    t_cluster->Write("", TObject::kOverwrite);   // sovrascrive se rilanci la macro

    int cluster_count = one_pixel_on + two_pixels_on + three_pixels_on + four_pixels_on;

    std::cout << "Numero di eventi con un solo pixel acceso: " << one_pixel_on << std::endl;
    std::cout << "Numero di eventi con due pixel accesi: " << two_pixels_on << std::endl;
    std::cout << "Numero di eventi con tre pixel accesi: " << three_pixels_on << std::endl;
    std::cout << "Numero di eventi con quattro pixel accesi: " << four_pixels_on << std::endl;
    std::cout << "Numero di eventi senza cluster: " << no_events << std::endl;
    std::cout << "Numero di eventi con cluster: " <<  cluster_count << std::endl;

    c2->cd();
    h_cluster->Draw();
    c2->Modified();
    c2->Update();

    // 4 plot in un solo canvas (conteggi reali, non normalizzati)
    for (int k = 0; k < 4; k++) {
        c4->cd(k + 1);
        signal_histograms[k]->Draw();
    }
    c4->Modified();
    c4->Update();

    //plot di sovrapposizione dei segnali dei cluster in funzione della dimensione del cluster
    TCanvas *c8 = new TCanvas("c8", "Segnale vs Dimensione del Cluster", 800, 600);

    // cloni normalizzati: gli originali (disegnati in c4) restano intatti
    TH1D *h_norm[4];
    for (int k = 0; k < 4; k++) {
        h_norm[k] = (TH1D*)signal_histograms[k]->Clone(Form("h_norm_%d", k + 1));
        h_norm[k]->SetDirectory(0);   
        const double integral = h_norm[k]->Integral();
        if (integral > 0.) h_norm[k]->Scale(1.0 / integral);
    }

    h_norm[0]->SetLineColor(kRed);
    h_norm[1]->SetLineColor(kGreen + 2);
    h_norm[2]->SetLineColor(kBlue);
    h_norm[3]->SetLineColor(kBlack);

    // massimo tra tutti gli istogrammi normalizzati
    double ymax = 0.;
    for (int k = 0; k < 4; k++) {
        if (h_norm[k]->GetMaximum() > ymax) ymax = h_norm[k]->GetMaximum();
    }

    // il primo istogramma disegnato definisce gli assi: massimo globale + 15% di margine
    h_norm[0]->SetTitle("Segnale totale normalizzato per cluster size;Segnale [mV];Frazione di eventi");
    h_norm[0]->SetMaximum(1.15 * ymax);
    h_norm[0]->SetMinimum(0.);

    c8->cd();
    h_norm[0]->Draw("HIST");
    h_norm[1]->Draw("HIST SAME");
    h_norm[2]->Draw("HIST SAME");
    h_norm[3]->Draw("HIST SAME");

    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(h_norm[0], "1 Pixel", "l");
    legend->AddEntry(h_norm[1], "2 Pixels", "l");
    legend->AddEntry(h_norm[2], "3 Pixels", "l");
    legend->AddEntry(h_norm[3], "4 Pixels", "l");
    legend->Draw();

    c8->Modified();
    c8->Update();


    //----------------------------------------------------------------------------------------------------------------------------------------------

    
    // =====================================================================
    // Distribuzioni di signal (pixel interni) e p2p (pixel esterni)
    // in funzione del cluster size (numero di pixel interni accesi)
    // =====================================================================
    {
        const int extPx[12] = {0, 1, 2, 3, 4, 7, 8, 11, 12, 13, 14, 15};
        const int cenPx[4]  = {5, 6, 9, 10};
        const int nCS = 4;   // cluster size 1..4

        // ---- canvas: uno per cluster size ----
        TCanvas *cCS[nCS];
        for (int s = 0; s < nCS; s++) {
            cCS[s] = new TCanvas(Form("cCS%d", s + 1), Form("P2P / Signal - %d pixel interni accesi", s + 1), 1200, 800);
            cCS[s]->Divide(4, 4);
        }

        // ---- tutto indicizzato per NUMERO DI PIXEL (0..15) ----
        TTree *tr[16] = {};
        double v_p2p[16], v_sig[16], v_base[16];
        for (int i = 0; i < 16; i++) { v_p2p[i] = 0; v_sig[i] = 0; v_base[i] = 0; }

        TH1D *h_cs[nCS][16] = {};   // p2p (esterni) o signal (centrali) per cluster size

        // ---- tree + istogrammi pixel esterni ----
        for (int k = 0; k < 12; k++) {
            int px = extPx[k];
            tr[px] = (TTree*)inputfile->Get(Form("pixel%d", px));
            if (!tr[px] || !tr[px]->GetBranch("p2p")) {
                std::cerr << "Manca tree/branch p2p per pixel " << px << std::endl;
                tr[px] = nullptr;
                continue;
            }
            tr[px]->SetBranchAddress("p2p", &v_p2p[px]);

            for (int s = 0; s < nCS; s++) {
                h_cs[s][px] = new TH1D(Form("h_cs_p2p_s%d_pixel%d", s + 1, px), Form("Pixel %d, cluster size = %d;P2P [ADC];Eventi", px, s + 1), 40, 0, 8.5);
                h_cs[s][px]->SetDirectory(nullptr);
            }
        }

        // ---- tree + istogrammi pixel centrali ----
        for (int k = 0; k < 4; k++) {
            int px = cenPx[k];
            tr[px] = (TTree*)inputfile->Get(Form("pixel%d", px));
            if (!tr[px] || !tr[px]->GetBranch("signal") || !tr[px]->GetBranch("baseline")) {
                std::cerr << "Manca tree/branch signal o baseline per pixel " << px << std::endl;
                tr[px] = nullptr;
                continue;
            }
            tr[px]->SetBranchAddress("signal",   &v_sig[px]);
            tr[px]->SetBranchAddress("baseline", &v_base[px]);

            for (int s = 0; s < nCS; s++) {
                h_cs[s][px] = new TH1D(Form("h_cs_sig_s%d_pixel%d", s + 1, px), Form("Pixel %d, cluster size = %d;signal [SCOPE];Eventi", px, s + 1), 100, 0, 100);
                h_cs[s][px]->SetDirectory(nullptr);
            }
        }

        // ---- un solo ciclo sugli eventi ----
        for (Long64_t ev = 0; ev < nEntries; ev++) {

            for (int px = 0; px < 16; px++)
                if (tr[px]) tr[px]->GetEntry(ev);

            // cluster size = numero di pixel interni accesi
            int csize = 0;
            for (int k = 0; k < 4; k++) {
                int px = cenPx[k];
                if (tr[px] && v_base[px] != -999.) csize++;
            }
            if (csize < 1 || csize > nCS) continue;   // nessun pixel acceso: scarta
            int sIdx = csize - 1;

            // pixel esterni: p2p
            for (int k = 0; k < 12; k++) {
                int px = extPx[k];
                if (!tr[px]) continue;
                h_cs[sIdx][px]->Fill(v_p2p[px]);
            }

            // pixel centrali: signal (solo se acceso)
            for (int k = 0; k < 4; k++) {
                int px = cenPx[k];
                if (!tr[px] || v_base[px] == -999.) continue;
                h_cs[sIdx][px]->Fill(v_sig[px]);
            }
        }

        // ---- disegno ----
        for (int s = 0; s < nCS; s++) {
            for (int px = 0; px < 16; px++) {
                if (!h_cs[s][px]) continue;
                cCS[s]->cd(px + 1);
                h_cs[s][px]->Draw();
            }
            cCS[s]->Modified();
            cCS[s]->Update();
        }
    }
        
}