#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <TProfile.h>
#include <TPad.h>

using namespace std;

// =====================================================================
// COSTANTI
// =====================================================================
const double PITCH = 10.;                   // [micron] dimensione del pixel 
const int    NPIX = 16;
const int    NCEN = 4;
const int    NEXT = 12;
const int    CEN_PX[NCEN] = {5, 6, 9, 10};
const int    EXT_PX[NEXT] = {0, 1, 2, 3, 4, 7, 8, 11, 12, 13, 14, 15};
const double NOHIT = -999.;     // baseline = NOHIT -> pixel spento
const char  *DATAFILE = "data2.root";
const double SEG_MAX = 300.;    // range istogrammi di "segnale" dei pixel esterni [ADC]
double       SEG_THR = 8.;     // soglia sul segnale dei pixel esterni [ADC]: si usa solo segnale > SEG_THR
double       SEG_MIN_PLOT = 10.;   // limite inferiore asse x dei segnali esterni
double       SEG_MAX_PLOT = 100.;  // limite superiore asse x dei segnali esterni

bool IsCentralPx(int px) {
    for (int k = 0; k < NCEN; k++) if (CEN_PX[k] == px) return true;
    return false;
}
int CenIndex(int px) {          // posizione in CEN_PX, -1 se esterno
    for (int k = 0; k < NCEN; k++) if (CEN_PX[k] == px) return k;
    return -1;
}

// =====================================================================
// CLASSE PIXEL: contiene tutti i dati di un pixel, evento per evento
//   - centrali: baseline + signal
//   - esterni:  segnale (ampiezza), presente solo se sopra soglia
// Dopo Expand() tutti i vettori hanno lunghezza nEv (indice = evento originale)
// =====================================================================
class Pixel {
    public:
    int  id      = -1;
    bool central = false;
    bool hasBaseline = false, hasSignal = false, hasSegnale = false;

    std::vector<double> baseline, signal, segnale;
    std::vector<char>   present;        // 1 se l'evento e' presente nel tree del pixel

    Long64_t maxEv = 0;                 // (massimo indice evento + 1)

    // dati grezzi come stanno nel tree (prima di Expand)
    std::vector<Long64_t> evRaw;
    std::vector<double>   bRaw, sRaw, gRaw;

    // legge il tree "pixelN" una sola volta; ritorna il numero di entries
    Long64_t Load(TFile *f, int pxId) {
        id      = pxId;
        central = IsCentralPx(pxId);
        evRaw.clear(); bRaw.clear(); sRaw.clear(); gRaw.clear();

        TTree *t = (TTree*)f->Get(Form("pixel%d", id));
        if (!t) { std::cerr << "Manca il tree pixel" << id << std::endl; return 0; }

        double b = 0., s = 0., g = 0.;
        Long64_t e = 0;
        bool hasEvent = t->GetBranch("event") != nullptr;
        hasBaseline = t->GetBranch("baseline") != nullptr;
        hasSignal   = t->GetBranch("signal")   != nullptr;
        hasSegnale  = t->GetBranch("segnale")  != nullptr;

        if (hasEvent)    t->SetBranchAddress("event",    &e);
        if (hasBaseline) t->SetBranchAddress("baseline", &b);
        if (hasSignal)   t->SetBranchAddress("signal",   &s);
        if (hasSegnale)  t->SetBranchAddress("segnale",  &g);

        Long64_t n = t->GetEntries();
        evRaw.reserve(n);
        if (hasBaseline) bRaw.reserve(n);
        if (hasSignal)   sRaw.reserve(n);
        if (hasSegnale)  gRaw.reserve(n);

        maxEv = 0;
        for (Long64_t i = 0; i < n; i++) {
            t->GetEntry(i);
            Long64_t idx = hasEvent ? e : i;     // senza "event" si assume allineato
            evRaw.push_back(idx);
            maxEv = std::max(maxEv, idx + 1);
            if (hasBaseline) bRaw.push_back(b);
            if (hasSignal)   sRaw.push_back(s);
            if (hasSegnale)  gRaw.push_back(g);
        }
        t->ResetBranchAddresses();   // b, s, g, e sono variabili locali
        return n;
    }

    // porta i dati su vettori lunghi nEv (indice = evento originale)
    void Expand(Long64_t nEv) {
        baseline.assign(nEv, NOHIT);
        signal.assign(nEv, 0.);
        segnale.assign(nEv, 0.);
        present.assign(nEv, 0);
        for (size_t i = 0; i < evRaw.size(); i++) {
            Long64_t ev = evRaw[i];
            if (ev < 0 || ev >= nEv) continue;
            present[ev] = 1;
            if (hasBaseline) baseline[ev] = bRaw[i];
            if (hasSignal)   signal[ev]   = sRaw[i];
            if (hasSegnale)  segnale[ev]  = gRaw[i];
        }
        evRaw.clear(); bRaw.clear(); sRaw.clear(); gRaw.clear();
    }

    // il pixel ha i dati che servono alle analisi?
    bool Valid() const { return central ? (hasBaseline && hasSignal) : hasSegnale; }

    // pixel centrale acceso
    bool IsOn(Long64_t ev) const {
        return central && hasBaseline && present[ev] && baseline[ev] != NOHIT;
    }

    // pixel esterno con segnale sopra soglia nell'evento
    // bool Fired(Long64_t ev) const { return !central && hasSegnale && present[ev]; }
    // pixel esterno con segnale sopra soglia nell'evento
    bool Fired(Long64_t ev) const {
        return !central && hasSegnale && present[ev] && segnale[ev] > SEG_THR;
    }

    // grandezza "da mostrare": signal per i centrali, segnale per gli esterni
    double Value(Long64_t ev) const { return central ? signal[ev] : segnale[ev]; }
};

// =====================================================================
// CLASSE DATASET: 16 pixel + informazioni di cluster calcolate in memoria
// =====================================================================
class Dataset {
    public:
    Pixel    px[NPIX];
    Long64_t nEv = 0;

    // info di cluster (una entry per evento)
    std::vector<int>                csize;    // numero di pixel centrali accesi
    std::vector<std::array<int,4>>  on;       // on[ev][i] = 1 se CEN_PX[i] acceso
    std::vector<double>             totSig;   // segnale totale del cluster

    bool Load(TFile *f) {
        bool any = false;
        for (int i = 0; i < NPIX; i++) {
            Long64_t n = px[i].Load(f, i);
            if (!px[i].Valid()) {
                std::cerr << "Pixel " << i << ": dati mancanti, escluso dalle analisi" << std::endl;
                continue;
            }
            any = true;
            nEv = std::max(nEv, px[i].maxEv);
            // std::cout << "Pixel " << i << ": " << n << " entries" << std::endl;
        }
        if (!any) { std::cerr << "Nessun pixel valido" << std::endl; return false; }

        for (int i = 0; i < NPIX; i++) if (px[i].Valid()) px[i].Expand(nEv);

        std::cout << "Numero totale di eventi da analizzare: " << nEv << std::endl;
        BuildClusters();
        return true;
    }

    void BuildClusters() {
        csize.assign(nEv, 0);
        on.assign(nEv, std::array<int,4>{0, 0, 0, 0});
        totSig.assign(nEv, 0.);
        for (Long64_t ev = 0; ev < nEv; ev++) {
            for (int i = 0; i < NCEN; i++) {
                const Pixel &p = px[CEN_PX[i]];
                if (p.Valid() && p.IsOn(ev)) {
                    on[ev][i] = 1;
                    csize[ev]++;
                    totSig[ev] += p.signal[ev];
                }
            }
        }
    }

    // indice (0..3) dell'unico pixel centrale acceso, -1 se csize != 1
    int SeedIndex(Long64_t ev) const {
        if (csize[ev] != 1) return -1;
        for (int i = 0; i < NCEN; i++) if (on[ev][i]) return i;
        return -1;
    }
};


// =====================================================================
// OPZIONI (valori di default) E DATASET GLOBALE
// =====================================================================
Dataset d;                                  // dati letti una sola volta
bool   optWriteClusterTree = false;
bool   optDoMap16          = false;
bool   optDoClusterCount   = false;
bool   optDoPerClusterSz   = false;
bool   optDoSeed           = false;
bool   optDoSegSelection   = false;
bool   optDoCentroid       = false;
bool   optDoOpposite       = false;

// legge data2.root una sola volta (le funzioni Plot* lo chiamano da sole se i dati non ci sono)
bool LoadData() {
    if (d.nEv > 0) return true;               // gia' caricato
    TFile *f = TFile::Open(DATAFILE, "READ");
    if (!f || f->IsZombie()) { std::cerr << "Impossibile aprire " << DATAFILE << std::endl; return false; }
    bool ok = d.Load(f);
    f->Close();   // da qui in poi si lavora solo sui dati in memoria
    return ok;
}

// =====================================================================
// FUNZIONI DI SUPPORTO (istogrammi per tutti i 16 pixel)
// =====================================================================

// crea un istogramma per ogni pixel valido: signal [0,100] per i centrali, segnale [0,SEG_MAX] per gli esterni
void MakeHists16(const char *name, const char *title, TH1D *h[NPIX]) {
    for (int i = 0; i < NPIX; i++) {
        h[i] = nullptr;
        if (!d.px[i].Valid()) continue;
        if (d.px[i].central)
            h[i] = new TH1D(Form("%s_pixel%d", name, i),
                            Form("Pixel %d (%s);signal [SCOPE];Eventi", i, title), 50, 0, 100);
        else
            h[i] = new TH1D(Form("%s_pixel%d", name, i),
                            Form("Pixel %d (%s);Segnale [ADC];Eventi", i, title), 50, 0, 100);
        h[i]->SetDirectory(nullptr);
    }
}

// riempie gli istogrammi con i valori dell'evento ev:
// esterni solo se sopra soglia (segnale), centrali solo se accesi (signal)
void FillEvent16(Long64_t ev, TH1D *h[NPIX]) {
    for (int i = 0; i < NPIX; i++) {
        if (!h[i]) continue;
        if (d.px[i].central) { if (!d.px[i].IsOn(ev)) continue; }
        else                 { if (!d.px[i].Fired(ev)) continue; }
        h[i]->Fill(d.px[i].Value(ev));
    }
}

// disegna i 16 istogrammi in un canvas 4x4 (pad = numero pixel + 1)
void Draw16(TCanvas *c, TH1D *h[NPIX]) {
    for (int i = 0; i < NPIX; i++) {
        if (!h[i]) continue;
        c->cd(i + 1);
        h[i]->Draw();
    }
    c->Modified();
    c->Update();
}

// =====================================================================
// ANALISI 1b: mappa 4x4 di tutti i pixel
//   centrali: acceso (baseline != NOHIT)
//   esterni : segnale > SEG_THR
// =====================================================================
void PlotMap16() {
    if (!LoadData()) return;

    TCanvas *c = new TCanvas("c1_16", "Mappa Pixel 4x4", 1000, 800);
    TH2D *h = new TH2D("h_map16",
                       Form("Stato Pixel (4x4), esterni con segnale > %g;Colonna;Riga", SEG_THR),
                       4, 0, 4, 4, 0, 4);
    h->SetDirectory(nullptr);
    gStyle->SetOptStat(0);

    for (Long64_t ev = 0; ev < d.nEv; ev++) {
        for (int i = 0; i < NPIX; i++) {
            if (!d.px[i].Valid()) continue;
            bool hit = d.px[i].central ? d.px[i].IsOn(ev) : d.px[i].Fired(ev);
            if (!hit) continue;
            // riga 0 in alto, come la numerazione dei pixel
            h->Fill(i % 4 + 0.5, 3 - i / 4 + 0.5);
        }
    }

    h->GetXaxis()->SetNdivisions(4);
    h->GetYaxis()->SetNdivisions(4);
    h->GetXaxis()->SetLabelSize(0);     // gli assi non hanno senso come numeri: ci sono le etichette
    h->GetYaxis()->SetLabelSize(0);
    h->SetMinimum(0);
    h->SetMarkerSize(1.6);
    c->cd();
    c->SetRightMargin(0.15);   // spazio per la palette (default 0.1)
    c->SetLeftMargin(0.10);
    h->Draw("COLZ TEXT");

    // etichetta con il numero del pixel in ogni cella
    for (int i = 0; i < NPIX; i++) {
        TLatex *label = new TLatex(i % 4 + 0.5, 3 - i / 4 + 0.78,
                                   Form("Pixel %d", i));
        label->SetTextAlign(22);
        label->SetTextSize(0.025);
        label->SetTextColor(kBlack);
        label->Draw();
    }
    c->Modified();
    c->Update();
    c->SaveAs("plots/mappa_4x4.pdf");
    c->SaveAs("plots/mappa_4x4.tex");
}

// =====================================================================
// ANALISI 2: numero di pixel accesi per evento + segnale totale vs cluster size
// =====================================================================
void PlotClusterCount() {
    if (!LoadData()) return;
    // ---- conteggi ----
    int count[NCEN + 1] = {0, 0, 0, 0, 0};
    for (Long64_t ev = 0; ev < d.nEv; ev++) count[d.csize[ev]]++;

    const char *nomi[NCEN + 1] = {"senza cluster", "un solo pixel acceso", "due pixel accesi",
                                  "tre pixel accesi", "quattro pixel accesi"};
    for (int k = 1; k <= NCEN; k++) std::cout << "Numero di eventi con " << nomi[k] << ": " << count[k] << std::endl;
    std::cout << "Numero di eventi " << nomi[0] << ": " << count[0] << std::endl;
    std::cout << "Numero di eventi con cluster: " << (d.nEv - count[0]) << std::endl;

    // ---- istogramma numero di pixel accesi ----
    TCanvas *c2 = new TCanvas("c2", "Cluster Count", 800, 600);
    TH1D *h_cluster = new TH1D("h_cluster", "Numero di pixel accesi per evento;Cluster Size;Conteggi",
                               5, -0.5, 4.5);
    h_cluster->SetDirectory(nullptr);
    for (int k = 0; k <= NCEN; k++) h_cluster->GetXaxis()->SetBinLabel(k + 1, Form("%d", k));
    for (Long64_t ev = 0; ev < d.nEv; ev++) h_cluster->Fill(d.csize[ev]);
    c2->cd();
    h_cluster->SetMarkerSize(1.5);
    h_cluster->Draw("HIST TEXT0");
    c2->Modified();
    c2->Update();
    c2->SaveAs("plots/pixel_on_per_event.pdf");
    c2->SaveAs("plots/pixel_on_per_event.tex");

    // ---- segnale totale per cluster size (4 pad) ----
    TCanvas *c4 = new TCanvas("c4", "Signal vs Cluster Size (1D)", 1200, 900);
    c4->Divide(2, 2);
    TH1D *hs[NCEN];
    for (int k = 0; k < NCEN; k++) {
        hs[k] = new TH1D(Form("h_signal_cluster_%d", k + 1),
                         Form("Segnale per cluster size = %d;Segnale [mV];Conteggi", k + 1), 100, 0, 100);
        hs[k]->SetDirectory(nullptr);
    }
    for (Long64_t ev = 0; ev < d.nEv; ev++)
        if (d.csize[ev] >= 1) hs[d.csize[ev] - 1]->Fill(d.totSig[ev]);
    for (int k = 0; k < NCEN; k++) { c4->cd(k + 1); hs[k]->Draw(); }
    c4->Modified();
    c4->Update();
    c4->SaveAs("plots/signtot_per_cls.pdf");
    c4->SaveAs("plots/signtot_per_cls.tex");

    // ---- sovrapposizione normalizzata ----
    TCanvas *c8 = new TCanvas("c8", "Segnale vs Dimensione del Cluster", 800, 600);
    TH1D *hn[NCEN];
    const int colori[NCEN] = {kRed, kGreen + 2, kBlue, kBlack};
    double ymax = 0.;
    for (int k = 0; k < NCEN; k++) {
        hn[k] = (TH1D*)hs[k]->Clone(Form("h_norm_%d", k + 1));
        hn[k]->SetDirectory(nullptr);
        double integral = hn[k]->Integral();
        if (integral > 0.) hn[k]->Scale(1.0 / integral);
        hn[k]->SetLineColor(colori[k]);
        ymax = std::max(ymax, hn[k]->GetMaximum());
    }
    hn[0]->SetTitle("Confronto segnale totale al variare di clsuster size;Segnale [mV];Conteggi / Conteggi totali");
    hn[0]->SetMaximum(1.15 * ymax);
    hn[0]->SetMinimum(0.);
    c8->cd();
    hn[0]->Draw("HIST");
    for (int k = 1; k < NCEN; k++) hn[k]->Draw("HIST SAME");

    TLegend *leg = new TLegend(0.1, 0.7, 0.3, 0.9);
    for (int k = 0; k < NCEN; k++) leg->AddEntry(hn[k], Form("%d Pixel%s", k + 1, k ? "s" : ""), "l");
    leg->Draw();
    c8->Modified();
    c8->Update();
    c8->SaveAs("plots/overlay_sign_vs_clsize.pdf");
    c8->SaveAs("plots/overlay_sign_vs_clsize.tex");
}

// =====================================================================
// ANALISI 2b: cluster size "completa" su tutti i 16 pixel
//   centrali: acceso (baseline != NOHIT)
//   esterni : segnale > SEG_THR
//   Produce l'istogramma del numero di pixel colpiti per evento (0..16)
//   e, opzionalmente, il segnale totale (centrali: signal, esterni: segnale)
//   per ogni cluster size.
// =====================================================================
void PlotClusterCountFull(bool plotSignal = false) {
    if (!LoadData()) return;

    // ---- calcolo cluster size completa ----
    std::vector<int>    csizeFull(d.nEv, 0);
    std::vector<double> totFull(d.nEv, 0.);

    for (Long64_t ev = 0; ev < d.nEv; ev++) {
        for (int i = 0; i < NPIX; i++) {
            if (!d.px[i].Valid()) continue;
            bool hit = d.px[i].central ? d.px[i].IsOn(ev) : d.px[i].Fired(ev);
            if (!hit) continue;
            csizeFull[ev]++;
            totFull[ev] += d.px[i].Value(ev);
        }
    }

    // ---- conteggi ----
    std::vector<int> count(NPIX + 1, 0);
    for (Long64_t ev = 0; ev < d.nEv; ev++) count[csizeFull[ev]]++;

    for (int k = 1; k <= NPIX; k++)
        if (count[k] > 0)
            std::cout << "Numero di eventi con cluster size " << k << ": " << count[k] << std::endl;
    std::cout << "Numero di eventi senza pixel colpiti: " << count[0] << std::endl;
    std::cout << "Numero di eventi con cluster (esterni sopra soglia " << SEG_THR << "): "
              << (d.nEv - count[0]) << std::endl;

    // ---- istogramma cluster size ----
    int kMax = 0;
    for (int k = NPIX; k >= 0; k--) if (count[k] > 0) { kMax = k; break; }

    TCanvas *c = new TCanvas("c2_full", "Cluster Count (16 pixel)", 900, 600);
    TH1D *h = new TH1D("h_cluster_full",
                       Form("Numero di pixel colpiti per evento (centrali accesi + esterni con segnale > %g);Pixel;Conteggi", SEG_THR),
                       kMax + 1, -0.5, kMax + 0.5);
    h->SetDirectory(nullptr);
    for (int k = 0; k <= kMax; k++) h->GetXaxis()->SetBinLabel(k + 1, Form("%d", k));
    for (Long64_t ev = 0; ev < d.nEv; ev++) h->Fill(csizeFull[ev]);

    h->SetMaximum(1.15 * h->GetMaximum());    // spazio per i numeri sopra i bin

    c->cd();
    h->SetMarkerSize(1.5);
    h->Draw("HIST TEXT0");                    // stile ROOT di default + conteggio sopra ogni bin
    c->Modified();
    c->Update();
    c->SaveAs("plots/pixel_on_per_event_full.pdf");
    c->SaveAs("plots/pixel_on_per_event_full.tex");

    if (!plotSignal) return;

    // ---- segnale totale per cluster size (raggruppato: 1,2,3,4,5,6+) ----
    const int NG = 6;
    TCanvas *c2 = new TCanvas("c2_full_sig", "Segnale totale vs cluster size (16 pixel)", 1200, 900);
    c2->Divide(3, 2);
    TH1D *hs[NG];
    for (int k = 0; k < NG; k++) {
        TString lab = (k < NG - 1) ? Form("= %d", k + 1) : Form(">= %d", k + 1);
        hs[k] = new TH1D(Form("h_signal_full_%d", k + 1),
                         Form("Segnale totale per cluster size %s;Segnale;Conteggi", lab.Data()),
                         100, 0, 300);
        hs[k]->SetDirectory(nullptr);
    }
    for (Long64_t ev = 0; ev < d.nEv; ev++) {
        int cs = csizeFull[ev];
        if (cs < 1) continue;
        hs[std::min(cs, NG) - 1]->Fill(totFull[ev]);
    }
    for (int k = 0; k < NG; k++) { c2->cd(k + 1); hs[k]->Draw(); }
    c2->Modified();
    c2->Update();
    c2->SaveAs("plots/signtot_per_cls_full.pdf");
    c2->SaveAs("plots/signtot_per_cls_full.tex");
}

// =====================================================================
// ANALISI 3: tutti i 16 pixel in funzione del cluster size (1..4)
// =====================================================================
void PlotPerClusterSize() {
    if (!LoadData()) return;
    TCanvas *cCS[NCEN];
    TH1D    *h[NCEN][NPIX];
    for (int s = 0; s < NCEN; s++) {
        cCS[s] = new TCanvas(Form("cCS%d", s + 1),
                             Form("Segnale - %d pixel interni accesi", s + 1), 1600, 1200);
        cCS[s]->Divide(4, 4);
        MakeHists16(Form("h_cs_s%d", s + 1), Form("cluster size = %d", s + 1), h[s]);
    }
    for (Long64_t ev = 0; ev < d.nEv; ev++) {
        int cs = d.csize[ev];
        if (cs < 1 || cs > NCEN) continue;
        FillEvent16(ev, h[cs - 1]);
    }
    for (int s = 0; s < NCEN; s++) {
        Draw16(cCS[s], h[s]);
        cCS[s]->SaveAs(Form("plots/segnale_clusize_%d.pdf", s + 1));   
        cCS[s]->SaveAs(Form("plots/segnale_clusize_%d.tex", s + 1));   
    }
    
}

// =====================================================================
// ANALISI 4: cluster size = 1, per ogni pixel centrale come seed.
// Taglio opzionale: signal del pixel acceso in [sigCut_min, sigCut_max]
// =====================================================================
void PlotSeed(bool applyCut = true, double sigCut_min = 0., double sigCut_max = 150.) {
    if (!LoadData()) return;

    TString cutTitle = "", cutName = "";
    if (applyCut) {
        cutTitle = Form(", signal in [%.0f, %.0f]", sigCut_min, sigCut_max);
        cutName  = Form("_sig%.0f_%.0f", sigCut_min, sigCut_max);
    }

    TCanvas *c[NCEN];
    TH1D    *h[NCEN][NPIX];
    int      nSel[NCEN] = {0, 0, 0, 0};

    for (int s = 0; s < NCEN; s++) {
        c[s] = new TCanvas(Form("c5_seed%d%s", CEN_PX[s], cutName.Data()),
                           Form("Solo pixel %d acceso%s", CEN_PX[s], cutTitle.Data()), 1200, 900);
        c[s]->Divide(4, 4);
        MakeHists16(Form("h_seed%d%s", CEN_PX[s], cutName.Data()),
                    Form("solo pixel %d acceso%s", CEN_PX[s], cutTitle.Data()), h[s]);
    }

    for (Long64_t ev = 0; ev < d.nEv; ev++) {
        int s = d.SeedIndex(ev);
        if (s < 0) continue;
        if (applyCut) {
            double sig = d.px[CEN_PX[s]].signal[ev];
            if (sig < sigCut_min || sig > sigCut_max) continue;
        }
        nSel[s]++;
        FillEvent16(ev, h[s]);
    }

    for (int s = 0; s < NCEN; s++) {
        std::cout << "Eventi con cluster size 1 (solo pixel " << CEN_PX[s] << ")"
                  << cutTitle << ": " << nSel[s] << std::endl;
        Draw16(c[s], h[s]);
    }
}

// =====================================================================
// ANALISI 5: eventi in cui il segnale di un pixel ESTERNO (selPx) e' in [segMin, segMax]
// (il pixel deve essere presente, cioe' sopra soglia)
// =====================================================================
void PlotSegnaleSelection(int selPx = 0, double segMin = 60., double segMax = 150.) {
    if (!LoadData()) return;

    if (selPx < 0 || selPx >= NPIX || IsCentralPx(selPx) || !d.px[selPx].Valid()) {
        std::cerr << "selPx deve essere un pixel esterno valido: 0,1,2,3,4,7,8,11,12,13,14,15" << std::endl;
        return;
    }

    TString selTitle = Form("segnale pixel %d in [%g, %g]", selPx, segMin, segMax);
    TString selName  = Form("_sel%d_seg_%g_%g", selPx, segMin, segMax);
    selName.ReplaceAll(".", "p");
    selName.ReplaceAll("-", "m");

    TCanvas *c = new TCanvas(Form("c6%s", selName.Data()), Form("Eventi con %s", selTitle.Data()), 1200, 900);
    c->Divide(4, 4);
    TH1D *h[NPIX];
    MakeHists16(Form("h_sel%s", selName.Data()), selTitle.Data(), h);

    int nSel = 0;
    for (Long64_t ev = 0; ev < d.nEv; ev++) {
        if (!d.px[selPx].Fired(ev)) continue;
        double g = d.px[selPx].segnale[ev];
        if (g < segMin || g > segMax) continue;
        nSel++;
        FillEvent16(ev, h);
    }
    std::cout << "Eventi con " << selTitle << ": " << nSel << std::endl;
    Draw16(c, h);
}

/*
// =====================================================================
// SCRITTURA DEL TREE "cluster" in data2.root (dai dati gia' in memoria)
// =====================================================================
void WriteClusterTree() {
    if (!LoadData()) return;
    TFile *f = TFile::Open(DATAFILE, "UPDATE");
    if (!f || f->IsZombie()) { std::cerr << "Impossibile aprire " << DATAFILE << " in UPDATE" << std::endl; return; }
    f->cd();
    TTree *t = new TTree("cluster", "Cluster dei pixel interni");

    Int_t    b_size = 0;
    Double_t b_totsig = 0.;
    Int_t    b_on[4]  = {0, 0, 0, 0};
    Double_t b_sig[4] = {0., 0., 0., 0.};

    t->Branch("cluster_size", &b_size,   "cluster_size/I");
    t->Branch("total_signal", &b_totsig, "total_signal/D");
    t->Branch("pixel_on",     b_on,      "pixel_on[4]/I");
    t->Branch("pixel_signal", b_sig,     "pixel_signal[4]/D");

    for (Long64_t ev = 0; ev < d.nEv; ev++) {
        b_size   = d.csize[ev];
        b_totsig = d.totSig[ev];
        for (int i = 0; i < NCEN; i++) {
            b_on[i]  = d.on[ev][i];
            b_sig[i] = d.on[ev][i] ? d.px[CEN_PX[i]].signal[ev] : 0.;
        }
        t->Fill();
    }
    t->Write("", TObject::kOverwrite);
    f->Close();
}
*/


// =====================================================================
// ANALISI 6: eventi con i due pixel centrali OPPOSTI (diagonale) accesi
//   coppia A: (5,10)   coppia B: (6,9)
//   exclusive = true  -> solo quei due pixel accesi (cluster size = 2),
//                        gli altri due centrali spenti
//   exclusive = false -> basta che i due pixel della coppia siano accesi
//                        (gli altri centrali possono essere accesi o no)
//   Un canvas 4x4 per coppia, con signal (centrali) / segnale (esterni)
//   di tutti i 16 pixel.
// =====================================================================
void PlotOppositePairs() {
    bool exclusive = true;
    if (!LoadData()) return;

    // indici in CEN_PX: 5->0, 6->1, 9->2, 10->3
    const int NPAIR = 2;
    const int pairs[NPAIR][2] = { {0, 3},     // (5,10)
                                  {1, 2} };   // (6,9)

    TCanvas *c[NPAIR];
    TH1D    *h[NPAIR][NPIX];
    int      nSel[NPAIR] = {0, 0};

    for (int p = 0; p < NPAIR; p++) {
        int pa = CEN_PX[pairs[p][0]], pb = CEN_PX[pairs[p][1]];
        TString suffix = exclusive ? "solo" : "almeno";
        c[p] = new TCanvas(Form("c7_opp%d_%d_%s", pa, pb, suffix.Data()),
                           Form("Pixel opposti %d e %d accesi (%s)", pa, pb, suffix.Data()),
                           1200, 900);
        c[p]->Divide(4, 4);
        MakeHists16(Form("h_opp%d_%d_%s", pa, pb, suffix.Data()),
                    Form("px %d e %d accesi%s", pa, pb, exclusive ? ", solo questi" : ""),
                    h[p]);
    }

    for (Long64_t ev = 0; ev < d.nEv; ev++) {
        if (exclusive && d.csize[ev] != 2) continue;
        for (int p = 0; p < NPAIR; p++) {
            if (!d.on[ev][pairs[p][0]] || !d.on[ev][pairs[p][1]]) continue;
            nSel[p]++;
            FillEvent16(ev, h[p]);
        }
    }

    for (int p = 0; p < NPAIR; p++) {
        int pa = CEN_PX[pairs[p][0]], pb = CEN_PX[pairs[p][1]];
        std::cout << "Eventi con pixel " << CEN_PX[pairs[p][0]] << " e " << CEN_PX[pairs[p][1]]
                  << " accesi" << (exclusive ? " (solo questi)" : "") << ": " << nSel[p] << std::endl;
        Draw16(c[p], h[p]);
        c[p]->Update();
        c[p]->SaveAs(Form("plots/opposti_%d_%d.pdf", pa, pb));      
        c[p]->SaveAs(Form("plots/opposti_%d_%d.tex", pa, pb));       
    }
}

// =====================================================================
// ANALISI 7: segnali di tutti i 16 pixel su tutti gli eventi (nessuna selezione)
//   centrali: signal, solo se il pixel e' acceso
//   esterni : segnale, solo se > SEG_THR
// =====================================================================
void PlotAllSignals(bool logY = false) {
    if (!LoadData()) return;

    TCanvas *c = new TCanvas("c_all_signals", "Segnali di tutti i pixel", 1600, 1200);
    c->Divide(4, 4);

    TH1D *h[NPIX];
    MakeHists16("h_all", "tutti gli eventi", h);

    for (Long64_t ev = 0; ev < d.nEv; ev++) FillEvent16(ev, h);

    Draw16(c, h);

    if (logY) {
        for (int i = 0; i < NPIX; i++) c->cd(i + 1)->SetLogy();
        c->Modified();
        c->Update();
    }

    c->SaveAs("plots/signal_all.pdf");
    c->SaveAs("plots/signal_all.tex");
}

// =====================================================================
// MACRO PRINCIPALE: esegue tutte le analisi con gli argomenti di default
// (le singole funzioni si possono richiamare da ROOT anche con altri argomenti)
// =====================================================================
void cls_analysis() {
    gStyle->SetLineScalePS(1);

    bool   optWriteClusterTree = false;
    if (!LoadData()) return;
    //if (optWriteClusterTree) WriteClusterTree();
    bool   optDoSeed           = false;      // check what happen selezionando un range nei pixel interni
    bool   optDoSegSelection   = false;      // check what happen selezionando un range nei pixel esterni
    if (optDoSeed)         PlotSeed(true, 50., 150);          // bool applyCut, double sigCut_min, double sigCut_max
    if (optDoSegSelection) PlotSegnaleSelection(0, 0., 150.); // pixel esterno 0, segnale in [60, 150]


    bool   optDoMap16          = true;      // mostra i conteggi di tutti i pixel 4x4
    bool   optDoClusterCount   = true;      // hist dimensione cls (2x2) + segnale sui cls
    bool   optDoClusterCountFull = true;     // cluster size su tutti i 16 pixel
    bool   optDoPerClusterSz   = true;      // segnale 16pixel divisi per cluster size (dei pixel interni)
    bool   optDoOpposite       = true;      // check eventi in cui pixel interni opposti si accendono
    bool   optDoAllSignals     = true;       // segnali dei 16 pixel su tutti gli eventi
    
    if (optDoMap16)        PlotMap16();
    if (optDoClusterCount) PlotClusterCount();
    if (optDoClusterCountFull) PlotClusterCountFull(false);   // true -> anche segnale totale
    if (optDoPerClusterSz) PlotPerClusterSize();
    if (optDoOpposite)     PlotOppositePairs();
    if (optDoAllSignals)   PlotAllSignals();   // true -> asse y logaritmico

    
}