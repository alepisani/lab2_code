#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2D.h>
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
    double baseline[4];
    double signal[4];

    for (int i = 0; i < 4; i++) {
        TString treename = TString::Format("pixel%d", p[i]);
        input_pixels[i] = (TTree*)inputfile->Get(treename);

        if (input_pixels[i]) {
            input_pixels[i]->SetBranchAddress("baseline", &baseline[i]);
        } else {
            std::cerr << "Attenzione: impossibile trovare il tree " << treename << std::endl;
        }
        if (input_pixels[i]) {
            input_pixels[i]->SetBranchAddress("signal", &signal[i]);
        } else {
            std::cerr << "Attenzione: impossibile trovare il tree " << treename << std::endl;
        }
    }

    // Assumiamo che tutti i tree abbiano lo stesso numero di entry (eventi)
    // Prendiamo come riferimento il primo albero valido trovato
    Long64_t nEntries = 0;
    for (int i = 0; i < 4; i++) {
        if (input_pixels[i]) {
            nEntries = input_pixels[i]->GetEntries();
            break;
        }
    }

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

        // Qui puoi aggiungere la logica per il controllo dei cluster (es. se pixel_accesi_in_questo_evento > 1)
        if (pixel_accesi_in_questo_evento > 1) {
            // Evento con cluster maggiore di 1!
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


    //-----------------------------------------------------
    //controllo se ci sono eventi in cui si sono accesi più di un pixel (cluster)
    int cluster_count = 0;
    
    // creazione istogramma per contare quanti eventi hanno registrato un cluster e quanti pixel si sono accesi in quell'evento
    TCanvas *c2 = new TCanvas("c2", "Cluster Count", 800, 600);
    TH1D *h_cluster = new TH1D("h_cluster", "Numero di pixel accesi per evento;Numero di pixel accesi;Numero di eventi", 5, -0.5, 4.5);
    h_cluster->GetXaxis()->SetBinLabel(1, "0");
    h_cluster->GetXaxis()->SetBinLabel(2, "1");
    h_cluster->GetXaxis()->SetBinLabel(3, "2");
    h_cluster->GetXaxis()->SetBinLabel(4, "3");
    h_cluster->GetXaxis()->SetBinLabel(5, "4"); 

    TCanvas *c3 = new TCanvas("c3", "Signal vs Cluster Size", 800, 600);
    TH2D *h_signal_cluster = new TH2D("h_signal_cluster", "Signal vs Cluster Size;Numero di pixel accesi;Segnale [mV]", 5, -0.5, 4.5, 100, 0, 175);

    //histogrammi per vedere segnale totale in funzione del numero di pixel accesi (cluster size)
    TCanvas *c4 = new TCanvas("c4", "Signal vs Cluster Size (1D)", 800, 600);
    TCanvas *c5 = new TCanvas("c5", "Signal vs Cluster Size (1D)", 800, 600);
    TCanvas *c6 = new TCanvas("c6", "Signal vs Cluster Size (1D)", 800, 600);
    TCanvas *c7 = new TCanvas("c7", "Signal vs Cluster Size (1D)", 800, 600);
    TH1D *h_signal_cluster_1 = new TH1D("h_signal_cluster_1", "Segnale per cluster size = 1;Segnale [mV];Numero di eventi", 100, 0, 100);
    TH1D *h_signal_cluster_2 = new TH1D("h_signal_cluster_2", "Segnale per cluster size = 2;Segnale [mV];Numero di eventi", 100, 0, 100);
    TH1D *h_signal_cluster_3 = new TH1D("h_signal_cluster_3", "Segnale per cluster size = 3;Segnale [mV];Numero di eventi", 100, 0, 100);
    TH1D *h_signal_cluster_4 = new TH1D("h_signal_cluster_4", "Segnale per cluster size = 4;Segnale [mV];Numero di eventi", 100, 0, 100);
    

    

    
    std::vector<int> cluster_sizes; // per memorizzare le dimensioni dei cluster
    for (Long64_t ev = 0; ev < nEntries; ev++) {
        int pixel_accesi_in_questo_evento = 0;  
        double total_signal = 0.; 

        for (int i = 0; i < 4; i++) {
            if (!input_pixels[i]) continue;
            

            input_pixels[i]->GetEntry(ev);
            if (baseline[i] != -999.) {
                pixel_accesi_in_questo_evento++;
                total_signal += signal[i];

                // riempio gli histogrammi per ogni dimensione di cluster
                if (pixel_accesi_in_questo_evento == 1) {
                    h_signal_cluster_1->Fill(signal[i]);
                } else if (pixel_accesi_in_questo_evento == 2) {
                    h_signal_cluster_2->Fill(signal[i]);
                } else if (pixel_accesi_in_questo_evento == 3) {
                    h_signal_cluster_3->Fill(signal[i]);
                } else if (pixel_accesi_in_questo_evento == 4) {
                    h_signal_cluster_4->Fill(signal[i]);
                }
            }
        }
        
        //considero la somma dei segnali dei pixel accesi in questo evento come il segnale totale del cluster
        h_signal_cluster->Fill(pixel_accesi_in_questo_evento, total_signal);
        h_cluster->Fill(pixel_accesi_in_questo_evento);

        if (pixel_accesi_in_questo_evento > 1) {
            cluster_count++;
            //std::cout << "Evento " << ev << " ha un cluster con " << pixel_accesi_in_questo_evento << " pixel accesi." << std::endl;
        }

    }

    std::cout << "Numero di eventi con cluster: " << cluster_count << std::endl;
    
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
    h_signal_cluster_1->Scale(1.0 / h_signal_cluster_1->Integral());
    h_signal_cluster_2->Scale(1.0 / h_signal_cluster_2->Integral());
    h_signal_cluster_3->Scale(1.0 / h_signal_cluster_3->Integral());
    h_signal_cluster_4->Scale(1.0 / h_signal_cluster_4->Integral());

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
}