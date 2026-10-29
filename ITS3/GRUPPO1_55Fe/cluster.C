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
    TFile *inputfile = TFile::Open("data.root", "READ");
    if (!inputfile || inputfile->IsZombie()) {
        std::cerr << "Errore nell'apertura del file!" << std::endl;
        return;
    }

    // I tuoi 4 pixel fisici mappati su una griglia 2x2
    // p[0] -> indice 0 (es. x=0, y=0)
    // p[1] -> indice 1 (es. x=1, y=0)
    // p[2] -> indice 2 (es. x=0, y=1)
    // p[3] -> indice 3 (es. x=1, y=1)
    const int p[4] = {5, 6, 9, 10};

    // 1. Crea un canvas quadrato
    TCanvas *c1 = new TCanvas("c1", "Mappa Pixel 2x2", 800, 600);
    
    // 2. Crea un istogramma 2D 2x2
    TH2D *h_map = new TH2D("h_map", "Stato Pixel (2x2);Asse X;Asse Y", 2, 0, 2, 2, 0, 2);
    
    gStyle->SetOptStat(0);
    //gStyle->SetPalette(kViridis);

    // 3. Carica i 4 TTree e imposta i branch della baseline
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

    // Assumiamo che tutti i tree abbiano lo stesso numero di entry (eventi)
    // Prendiamo come riferimento il primo albero valido.
    TTree *reference_tree = nullptr;
    for (TTree *tree : input_pixels) {
        if (tree) {
            reference_tree = tree;
            break;
        }
    }
    if (!reference_tree) {
        std::cerr << "Errore: nessun tree pixel valido da analizzare." << std::endl;
        return;
    }
    Long64_t nEntries = reference_tree->GetEntries();

    std::cout << "Numero totale di eventi da analizzare: " << nEntries << std::endl;

    // 4. Ciclo sugli eventi (entry per entry)
    for (Long64_t ev = 0; ev < nEntries; ev++) {
        
        // Contatore per capire quanti pixel si sono accesi in questo singolo evento (utile per i cluster)
        int pixel_accesi_in_questo_evento = 0;

        for (int i = 0; i < 4; i++) {
            if (!input_pixels[i]) continue;

            // Legge l'evento corrente per l'albero i-esimo
            input_pixels[i]->GetEntry(ev);

            // Condizione: se la baseline è diversa da -999, il pixel ha registrato un evento
            if (baseline[i] != -999.) {
                int x = i % 2; // Colonna (0 o 1)
                int y = i / 2; // Riga (0 o 1)
                
                // Riempie la mappa incrementando il conteggio per quel pixel
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
    c1->Modified();
    c1->Update();


    //---------------------------------------------------------------------------------------------------------------------------------------------

    // creazione istogramma per contare quanti eventi hanno registrato un cluster e quanti pixel si sono accesi in quell'evento
    TCanvas *c2 = new TCanvas("c2", "Cluster Count", 800, 600);
    TH1D *h_cluster = new TH1D("h_cluster", "Numero di pixel accesi per evento;Numero di pixel accesi;Numero di eventi", 5, -0.5, 4.5);
    h_cluster->GetXaxis()->SetBinLabel(1, "0");
    h_cluster->GetXaxis()->SetBinLabel(2, "1");
    h_cluster->GetXaxis()->SetBinLabel(3, "2");
    h_cluster->GetXaxis()->SetBinLabel(4, "3");
    h_cluster->GetXaxis()->SetBinLabel(5, "4"); 

    TCanvas *c3 = new TCanvas("c3", "Cluster vs Size Signal", 800, 600);
    TH2D *h_signal_cluster = new TH2D("h_signal_cluster", "Cluster Size vs Signal;Segnale [mV];Numero di pixel accesi;", 100, 0, 100, 4, 0.5, 4.5);
    h_signal_cluster->GetYaxis()->SetBinLabel(1, "1");
    h_signal_cluster->GetYaxis()->SetBinLabel(2, "2");
    h_signal_cluster->GetYaxis()->SetBinLabel(3, "3");
    h_signal_cluster->GetYaxis()->SetBinLabel(4, "4");

    //histogrammi per vedere segnale totale in funzione del numero di pixel accesi (cluster size)
    TCanvas *c4 = new TCanvas("c4", "Signal vs Cluster Size (1D)", 800, 600);
    TCanvas *c5 = new TCanvas("c5", "Signal vs Cluster Size (1D)", 800, 600);
    TCanvas *c6 = new TCanvas("c6", "Signal vs Cluster Size (1D)", 800, 600);
    TCanvas *c7 = new TCanvas("c7", "Signal vs Cluster Size (1D)", 800, 600);
    TH1D *h_signal_cluster_1 = new TH1D("h_signal_cluster_1", "Segnale per cluster size = 1;Segnale [mV];Numero di eventi", 100, 0, 100);
    TH1D *h_signal_cluster_2 = new TH1D("h_signal_cluster_2", "Segnale per cluster size = 2;Segnale [mV];Numero di eventi", 100, 0, 100);
    TH1D *h_signal_cluster_3 = new TH1D("h_signal_cluster_3", "Segnale per cluster size = 3;Segnale [mV];Numero di eventi", 100, 0, 100);
    TH1D *h_signal_cluster_4 = new TH1D("h_signal_cluster_4", "Segnale per cluster size = 4;Segnale [mV];Numero di eventi", 100, 0, 100);
    

    int no_events = 0; 
    int one_pixel_on = 0;
    int two_pixels_on = 0;
    int three_pixels_on = 0;
    int four_pixels_on = 0;
    
    for (Long64_t ev = 0; ev < nEntries; ev++) {
        int pixel_accesi_in_questo_evento = 0;  
        double total_signal = 0.; 


        // --- PASSO 1: Conta quanti pixel si sono accesi in questo evento ---
        for (int i = 0; i < 4; i++) {
            if (!input_pixels[i]) continue;
            
            input_pixels[i]->GetEntry(ev); // Legge l'evento
            
            if (baseline[i] != -999.) {
                pixel_accesi_in_questo_evento++;
                total_signal += signal[i];
            }
        }
        
        // Riempie gli istogrammi generali dell'evento
        h_cluster->Fill(pixel_accesi_in_questo_evento);

        if(pixel_accesi_in_questo_evento > 0)h_signal_cluster->Fill(total_signal, pixel_accesi_in_questo_evento);
        
        if(!pixel_accesi_in_questo_evento) no_events++;
        if(pixel_accesi_in_questo_evento == 1) one_pixel_on++;
        if(pixel_accesi_in_questo_evento == 2) two_pixels_on++;
        if(pixel_accesi_in_questo_evento == 3) three_pixels_on++;
        if(pixel_accesi_in_questo_evento == 4) four_pixels_on++;


        // --- PASSO 2: Rileggi i pixel e riempi l'istogramma 1D con il *singolo* segnale ---
        // (Lo facciamo solo se l'evento ha una dimensione di cluster valida)
        for (int i = 0; i < 4; i++) {
            if (!input_pixels[i]) continue;
            
            input_pixels[i]->GetEntry(ev); // Rileggiamo lo stesso evento in modo sicuro
            
            if (baseline[i] != -999.) {
                // Ora sappiamo con certezza qual è la dimensione totale del cluster (pixel_accesi_in_questo_evento)
                // e prendiamo esattamente il signal[i] di questo specifico pixel.
                
                if (pixel_accesi_in_questo_evento == 1) {
                    h_signal_cluster_1->Fill(signal[i]);
                } 
                else if (pixel_accesi_in_questo_evento == 2) {
                    h_signal_cluster_2->Fill(signal[i]);
                } 
                else if (pixel_accesi_in_questo_evento == 3) {
                    h_signal_cluster_3->Fill(signal[i]);
                } 
                else if (pixel_accesi_in_questo_evento == 4) {
                    h_signal_cluster_4->Fill(signal[i]);
                }
            }
        }
    }

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

    c3->cd();
    h_signal_cluster->Draw("COLZ");
    c3->Modified();
    c3->Update();

    c4->cd();
    h_signal_cluster_1->Draw();
    c4->Modified();
    c4->Update();

    c5->cd();
    h_signal_cluster_2->Draw();
    c5->Modified();
    c5->Update();

    c6->cd();
    h_signal_cluster_3->Draw();
    c6->Modified();
    c6->Update();

    c7->cd();
    h_signal_cluster_4->Draw();
    c7->Modified();
    c7->Update();      

    //plot di sovrapposizione dei segnali dei cluster in funzione della dimensione del cluster
    TCanvas *c8 = new TCanvas("c8", "Segnale vs Dimensione del Cluster", 800, 600);
    //normalizzo gli istogrammi per poterli confrontare meglio
    TH1D *signal_histograms[] = {
        h_signal_cluster_1, h_signal_cluster_2,
        h_signal_cluster_3, h_signal_cluster_4
    };
    for (TH1D *histogram : signal_histograms) {
        const double integral = histogram->Integral();
        if (integral > 0.) {
            histogram->Scale(1.0 / integral);
        }
    }

    h_signal_cluster_1->SetLineColor(kRed);
    h_signal_cluster_2->SetLineColor(kGreen);
    h_signal_cluster_3->SetLineColor(kBlue);
    h_signal_cluster_4->SetLineColor(kBlack);

    h_signal_cluster_1->Draw("HIST");
    h_signal_cluster_2->Draw("HIST SAME");
    h_signal_cluster_3->Draw("HIST SAME");
    h_signal_cluster_4->Draw("HIST SAME");

    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(h_signal_cluster_1, "1 Pixel", "l");
    legend->AddEntry(h_signal_cluster_2, "2 Pixels", "l");
    legend->AddEntry(h_signal_cluster_3, "3 Pixels", "l");
    legend->AddEntry(h_signal_cluster_4, "4 Pixels", "l");
    legend->Draw();

    c8->Modified();
    c8->Update();


    //----------------------------------------------------------------------------------------------------------------------------------------------

    /*
    //esploriamo gli altri pixel 0,1,2,3,4,7,8,11,12,13,14,15
    //p2p è il segnale di ampiezza ADC
    //voglio mettere su un unico canvas tutti i p2p dei pixel in plot diversi.

    const int pixel[12] = {0, 1, 2, 3, 4, 7, 8, 11, 12, 13, 14, 15};
    double p2p[12] = {};
    std::vector<TTree*> p2p_trees(12, nullptr);
    std::vector<TH1D*> h_p2p(12, nullptr);
    std::vector<TH1D*> h_p2p_range1(12, nullptr);
    std::vector<TH1D*> h_p2p_range2(12, nullptr);
    std::vector<TH1D*> h_p2p_range3(12, nullptr);

    for (int i = 0; i < 12; i++) {
        TString treename = TString::Format("pixel%d", pixel[i]);
        p2p_trees[i] = (TTree*)inputfile->Get(treename);

        if (p2p_trees[i] && p2p_trees[i]->GetBranch("p2p")) {
            p2p_trees[i]->SetBranchAddress("p2p", &p2p[i]);
        } else {
            std::cerr << "Attenzione: impossibile trovare il tree o il branch p2p di "
                      << treename << std::endl;
            p2p_trees[i] = nullptr;
        }

    }

    TCanvas *c9 = new TCanvas("c9", "P2P Signals", 1200, 800);
    c9->Divide(4, 4);
    TCanvas *c10 = new TCanvas("c10", "P2P Signals - Range 1", 1200, 800);
    c10->Divide(4, 4);
    TCanvas *c11 = new TCanvas("c11", "P2P Signals - Range 2", 1200, 800);
    c11->Divide(4, 4);
    TCanvas *c12 = new TCanvas("c12", "P2P Signals - Range 3", 1200, 800);
    c12->Divide(4, 4);

    double ref_p2p = 0;

    //pixel esterni
    for (int i = 0; i < 12; i++) {
        c9->cd(pixel[i] + 1);
        h_p2p[i] = new TH1D(TString::Format("h_p2p_pixel%d", pixel[i]), TString::Format("Pixel %d;P2P [ADC];Numero di eventi", pixel[i]), 40, 2, 8.5);
        h_p2p[i]->SetDirectory(nullptr);
        h_p2p_range1[i] = new TH1D(TString::Format("h_p2p_range1_pixel%d", pixel[i]), TString::Format("Pixel %d;P2P [ADC];Numero di eventi", pixel[i]), 40, 2, 8.5);
        h_p2p_range1[i]->SetDirectory(nullptr);
        h_p2p_range2[i] = new TH1D(TString::Format("h_p2p_range2_pixel%d", pixel[i]), TString::Format("Pixel %d;P2P [ADC];Numero di eventi", pixel[i]), 40, 2, 8.5);
        h_p2p_range2[i]->SetDirectory(nullptr);
        h_p2p_range3[i] = new TH1D(TString::Format("h_p2p_range3_pixel%d", pixel[i]), TString::Format("Pixel %d;P2P [ADC];Numero di eventi", pixel[i]), 40, 2, 8.5);
        h_p2p_range3[i]->SetDirectory(nullptr);

        if (p2p_trees[i] && p2p_trees[i]->GetBranch("p2p")) {
            for (Long64_t ev = 0; ev < p2p_trees[i]->GetEntries(); ev++) {
                p2p_trees[i]->GetEntry(ev);
                h_p2p[i]->Fill(p2p[i]);
                if (p2p[i] >= 0 && p2p[i] < 4) {
                    c10->cd(pixel[i] + 1);
                    h_p2p_range1[i]->Fill(p2p[i]);
                }
                if (p2p[i] >= 4 && p2p[i] < 5.5) {
                    c11->cd(pixel[i] + 1);
                    h_p2p_range2[i]->Fill(p2p[i]);
                }
                if (p2p[i] >= 5.5 && p2p[i] < 8) {
                    c12->cd(pixel[i] + 1);
                    h_p2p_range3[i]->Fill(p2p[i]);
                }
                c9->cd(pixel[i] + 1);
            }
        }
        h_p2p[i]->Draw();
        if(i == 0) ref_p2p = 
    }

    //pixel interni
    for (int i = 0; i < 4; i++){
        c9->cd(p[i] + 1);
        h_p2p[i] = new TH1D(TString::Format("h_p2p_pixel%d", p[i]), TString::Format("Pixel %d;signal [SCOPE];Numero di eventi", p[i]), 100, 0, 100);
        h_p2p[i]->SetDirectory(nullptr);
        h_p2p_range1[i] = new TH1D(TString::Format("h_p2p_range1_pixel%d", p[i]), TString::Format("Pixel %d;P2P [ADC];Numero di eventi", p[i]), 40, 2, 8.5);
        h_p2p_range1[i]->SetDirectory(nullptr);
        h_p2p_range2[i] = new TH1D(TString::Format("h_p2p_range2_pixel%d", p[i]), TString::Format("Pixel %d;P2P [ADC];Numero di eventi", p[i]), 40, 2, 8.5);
        h_p2p_range2[i]->SetDirectory(nullptr);
        h_p2p_range3[i] = new TH1D(TString::Format("h_p2p_range3_pixel%d", p[i]), TString::Format("Pixel %d;P2P [ADC];Numero di eventi", p[i]), 40, 2, 8.5);
        h_p2p_range3[i]->SetDirectory(nullptr);

        if (input_pixels[i] &&
            input_pixels[i]->GetBranch("signal") &&
            input_pixels[i]->GetBranch("baseline")) {
            for (Long64_t ev = 0; ev < input_pixels[i]->GetEntries(); ev++) {
                input_pixels[i]->GetEntry(ev);
                
                double ref_p2p[12] = {};
                input_pixels[0]->GetEntry(ev);
                input_pixels[0]->SetBranchAddress("p2p", &ref_p2p[1]);

                if (baseline[i] != -999.) h_p2p[i]->Fill(signal[i]);
                if (baseline[i] != -999. && ref_p2p[1] >= 0. && ref_p2p[1] < 4.) {
                    c10->cd(p[i] + 1);
                    h_p2p_range1[i]->Fill(ref_p2p[1]);
                }
                if (baseline[i] != -999. && ref_p2p[1] >= 4. && ref_p2p[1] < 5.5) {
                    c11->cd(p[i] + 1);
                    h_p2p_range2[i]->Fill(ref_p2p[1]);
                }
                if (baseline[i] != -999. && ref_p2p[1] >= 5.5 && ref_p2p[1] < 8.) {
                    c12->cd(p[i] + 1);
                    h_p2p_range3[i]->Fill(ref_p2p[1]);
                }
                c9->cd(p[i] + 1);
            }
        }
        c9->cd(p[i] + 1);
        h_p2p[i]->Draw();
        c10->cd(p[i] + 1);
        h_p2p_range1[i]->Draw();
        c11->cd(p[i] + 1);
        h_p2p_range2[i]->Draw();
        c12->cd(p[i] + 1);
        h_p2p_range3[i]->Draw();

    }
    */


    const int ext[12] = {0, 1, 2, 3, 4, 7, 8, 11, 12, 13, 14, 15};
    const int cen[4]  = {5, 6, 9, 10};
    const int refPixel = 0;          // pixel usato per selezionare i range di p2p

    const int nR = 3;
    const double lo[nR] = {0., 4.,  5.5};
    const double hi[nR] = {4., 5.5, 8.};

    // ---- canvas ----
    TCanvas *c9 = new TCanvas("c9", "P2P / Signal", 1200, 800);
    c9->Divide(4, 4);
    TCanvas *cR[nR];
    for (int r = 0; r < nR; r++) {
        cR[r] = new TCanvas(Form("c%d", 10 + r),
                            Form("p2p(pixel %d) in [%.1f, %.1f)", refPixel, lo[r], hi[r]), 1200, 800);
        cR[r]->Divide(4, 4);
    }

    // ---- tutto indicizzato per NUMERO DI PIXEL (0..15) ----
    TTree *tree[16] = {};
    double p2p[16], sig[16], base[16];
    for (int i = 0; i < 16; i++) { p2p[i] = 0; sig[i] = 0; base[i] = 0; }

    TH1D *h_tot[16]     = {};          // p2p (esterni) o signal (centrali)
    TH1D *h_rng[nR][16] = {};          // idem, selezionato sul range di p2p di refPixel

    // ---- tree + istogrammi pixel esterni ----
    for (int k = 0; k < 12; k++) {
        int px = ext[k];
        tree[px] = (TTree*)inputfile->Get(Form("pixel%d", px));
        if (!tree[px] || !tree[px]->GetBranch("p2p")) {
            std::cerr << "Manca tree/branch p2p per pixel " << px << std::endl;
            tree[px] = nullptr;
            continue;
        }
        tree[px]->SetBranchAddress("p2p", &p2p[px]);

        h_tot[px] = new TH1D(Form("h_p2p_pixel%d", px),
                            Form("Pixel %d;P2P [ADC];Eventi", px), 40, 0, 8.5);
        h_tot[px]->SetDirectory(nullptr);
        for (int r = 0; r < nR; r++) {
            h_rng[r][px] = new TH1D(Form("h_p2p_r%d_pixel%d", r, px),
                                    Form("Pixel %d, p2p ref [%.1f,%.1f);P2P [ADC];Eventi", px, lo[r], hi[r]),
                                    40, 0, 8.5);
            h_rng[r][px]->SetDirectory(nullptr);
        }
    }

    // ---- tree + istogrammi pixel centrali ----
    for (int k = 0; k < 4; k++) {
        int px = cen[k];
        tree[px] = (TTree*)inputfile->Get(Form("pixel%d", px));
        if (!tree[px] || !tree[px]->GetBranch("signal") || !tree[px]->GetBranch("baseline")) {
            std::cerr << "Manca tree/branch signal o baseline per pixel " << px << std::endl;
            tree[px] = nullptr;
            continue;
        }
        tree[px]->SetBranchAddress("signal",   &sig[px]);
        tree[px]->SetBranchAddress("baseline", &base[px]);

        h_tot[px] = new TH1D(Form("h_sig_pixel%d", px),
                            Form("Pixel %d;signal [SCOPE];Eventi", px), 100, 0, 100);
        h_tot[px]->SetDirectory(nullptr);
        for (int r = 0; r < nR; r++) {
            h_rng[r][px] = new TH1D(Form("h_sig_r%d_pixel%d", r, px),
                                    Form("Pixel %d, p2p ref [%.1f,%.1f);signal [SCOPE];Eventi", px, lo[r], hi[r]),
                                    100, 0, 100);
            h_rng[r][px]->SetDirectory(nullptr);
        }
    }

    // ---- controlli ----
    if (!tree[refPixel]) { std::cerr << "Tree di riferimento mancante" << std::endl; return; }
    Long64_t nEv = tree[refPixel]->GetEntries();
    for (int px = 0; px < 16; px++)
        if (tree[px] && tree[px]->GetEntries() != nEv)
            std::cerr << "ATTENZIONE: pixel " << px << " ha " << tree[px]->GetEntries()
                    << " entry invece di " << nEv << std::endl;

    // ================= UN SOLO CICLO SUGLI EVENTI =================
    for (Long64_t ev = 0; ev < nEv; ev++) {

        // leggo l'evento in tutti i tree
        for (int px = 0; px < 16; px++)
            if (tree[px]) tree[px]->GetEntry(ev);

        double ref = p2p[refPixel];

        // indice del range in cui cade ref (-1 se fuori da tutti)
        int rIdx = -1;
        for (int r = 0; r < nR; r++)
            if (ref >= lo[r] && ref < hi[r]) { rIdx = r; break; }

        // pixel esterni: p2p
        for (int k = 0; k < 12; k++) {
            int px = ext[k];
            if (!tree[px]) continue;
            h_tot[px]->Fill(p2p[px]);
            if (rIdx >= 0) h_rng[rIdx][px]->Fill(p2p[px]);
        }

        // pixel centrali: signal (solo se acceso)
        for (int k = 0; k < 4; k++) {
            int px = cen[k];
            if (!tree[px] || base[px] == -999.) continue;
            h_tot[px]->Fill(sig[px]);
            if (rIdx >= 0) h_rng[rIdx][px]->Fill(sig[px]);
        }
    }

    // ================= DISEGNO =================
    for (int px = 0; px < 16; px++) {
        if (!h_tot[px]) continue;
        c9->cd(px + 1);
        h_tot[px]->Draw();
        for (int r = 0; r < nR; r++) {
            cR[r]->cd(px + 1);
            h_rng[r][px]->Draw();
        }
    }
    c9->Update();
    for (int r = 0; r < nR; r++) cR[r]->Update();

}