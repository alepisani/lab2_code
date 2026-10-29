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
    TFile *inputfile = TFile::Open("data.root", "UPDATE");
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
                
                /*
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
                */
               // --- Segnale totale (somma dei pixel accesi) in funzione del cluster size ---
                if (pixel_accesi_in_questo_evento == 1)      h_signal_cluster_1->Fill(total_signal);
                else if (pixel_accesi_in_questo_evento == 2) h_signal_cluster_2->Fill(total_signal);
                else if (pixel_accesi_in_questo_evento == 3) h_signal_cluster_3->Fill(total_signal);
                else if (pixel_accesi_in_questo_evento == 4) h_signal_cluster_4->Fill(total_signal);
            }   // fine ciclo sugli eventi
            
        }
    }

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

    // massimo tra tutti gli istogrammi normalizzati
    double ymax = 0.;
    for (TH1D *histogram : signal_histograms) {
        if (histogram->GetMaximum() > ymax) ymax = histogram->GetMaximum();
    }

    // l'istogramma disegnato per primo definisce gli assi: gli do il massimo globale + 15% di margine
    h_signal_cluster_1->SetMaximum(1.15 * ymax);
    h_signal_cluster_1->SetMinimum(0.);

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
        cR[r] = new TCanvas(Form("c%d", 10 + r), Form("p2p(pixel %d) in [%.1f, %.1f)", refPixel, lo[r], hi[r]), 1200, 800);
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

        h_tot[px] = new TH1D(Form("h_p2p_pixel%d", px), Form("Pixel %d;P2P [ADC];Eventi", px), 40, 0, 8.5);
        h_tot[px]->SetDirectory(nullptr);
        for (int r = 0; r < nR; r++) {
            h_rng[r][px] = new TH1D(Form("h_p2p_r%d_pixel%d", r, px), Form("Pixel %d, p2p ref [%.1f,%.1f);P2P [ADC];Eventi", px, lo[r], hi[r]), 40, 0, 8.5);
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

        h_tot[px] = new TH1D(Form("h_sig_pixel%d", px), Form("Pixel %d;signal [SCOPE];Eventi", px), 100, 0, 100);
        h_tot[px]->SetDirectory(nullptr);
        for (int r = 0; r < nR; r++) {
            h_rng[r][px] = new TH1D(Form("h_sig_r%d_pixel%d", r, px), Form("Pixel %d, p2p ref [%.1f,%.1f);signal [SCOPE];Eventi", px, lo[r], hi[r]), 100, 0, 100);
            h_rng[r][px]->SetDirectory(nullptr);
        }
    }


    
    // ================= UN SOLO CICLO SUGLI EVENTI =================
    
    Long64_t nEv = tree[refPixel]->GetEntries();
    
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

    */

    //------------------------------------------------------------------------------------------------------------------------------

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
            cCS[s] = new TCanvas(Form("cCS%d", s + 1),
                                 Form("P2P / Signal - %d pixel interni accesi", s + 1),
                                 1200, 800);
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
                h_cs[s][px] = new TH1D(Form("h_cs_p2p_s%d_pixel%d", s + 1, px),
                                       Form("Pixel %d, cluster size = %d;P2P [ADC];Eventi", px, s + 1),
                                       40, 0, 8.5);
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
                h_cs[s][px] = new TH1D(Form("h_cs_sig_s%d_pixel%d", s + 1, px),
                                       Form("Pixel %d, cluster size = %d;signal [SCOPE];Eventi", px, s + 1),
                                       100, 0, 100);
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