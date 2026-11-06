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
#include <TProfile.h>
#include <TPad.h>

using namespace std;

// =====================================================================
// COSTANTI
// =====================================================================
const int    NPIX = 16;
const int    NCEN = 4;
const int    NEXT = 12;
const int    CEN_PX[NCEN] = {5, 6, 9, 10};
const int    EXT_PX[NEXT] = {0, 1, 2, 3, 4, 7, 8, 11, 12, 13, 14, 15};
const double NOHIT = -999.;     // baseline = NOHIT -> pixel spento

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
//   - esterni:  p2p
// =====================================================================
class Pixel {
    public:
    int  id      = -1;
    bool central = false;
    bool hasBaseline = false, hasSignal = false, hasP2P = false;
    
    std::vector<double> baseline, signal, p2p;
    
    // legge il tree "pixelN" una sola volta; ritorna il numero di entries
    Long64_t Load(TFile *f, int pxId) {
        id      = pxId;
        central = IsCentralPx(pxId);
        baseline.clear(); signal.clear(); p2p.clear();
        
        TTree *t = (TTree*)f->Get(Form("pixel%d", id));
        if (!t) { std::cerr << "Manca il tree pixel" << id << std::endl; return 0; }
        
        double b = 0., s = 0., p = 0.;
        hasBaseline = t->GetBranch("baseline") != nullptr;
        hasSignal   = t->GetBranch("signal")   != nullptr;
        hasP2P      = t->GetBranch("p2p")      != nullptr;
        
        if (hasBaseline) t->SetBranchAddress("baseline", &b);
        if (hasSignal)   t->SetBranchAddress("signal",   &s);
        if (hasP2P)      t->SetBranchAddress("p2p",      &p);
        
        Long64_t n = t->GetEntries();
        if (hasBaseline) baseline.reserve(n);
        if (hasSignal)   signal.reserve(n);
        if (hasP2P)      p2p.reserve(n);
        
        for (Long64_t i = 0; i < n; i++) {
            t->GetEntry(i);
            if (hasBaseline) baseline.push_back(b);
            if (hasSignal)   signal.push_back(s);
            if (hasP2P)      p2p.push_back(p);
        }
        t->ResetBranchAddresses();   // b, s, p sono variabili locali
        return n;
    }
    
    // il pixel ha i dati che servono alle analisi?
    bool Valid() const { return central ? (hasBaseline && hasSignal) : hasP2P; }
    
    // pixel acceso (solo centrali)
    bool IsOn(Long64_t ev) const { return central && hasBaseline && baseline[ev] != NOHIT; }
    
    // grandezza "da mostrare": signal per i centrali, p2p per gli esterni
    double Value(Long64_t ev) const { return central ? signal[ev] : p2p[ev]; }
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
        bool first = true;
        for (int i = 0; i < NPIX; i++) {
            Long64_t n = px[i].Load(f, i);
            if (!px[i].Valid()) {
                std::cerr << "Pixel " << i << ": dati mancanti, escluso dalle analisi" << std::endl;
                continue;
            }
            if (first) { nEv = n; first = false; }
            else if (n != nEv) {
                std::cerr << "Attenzione: pixel " << i << " ha " << n << " entries (altri: " << nEv
                << "), uso il minimo" << std::endl;
                nEv = std::min(nEv, n);
            }
        }
        if (first) { std::cerr << "Nessun pixel valido" << std::endl; return false; }
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
const bool   optWriteClusterTree = false;   
const bool   optDoMap          = false;
const bool   optDoClusterCount = false;
const bool   optDoPerClusterSz = false;
const bool   optDoSeed         = false;
const bool   optDoP2PSelection = false;
const bool   optDoCentroid = false;
const bool   optDoOpposite = true;

// legge data.root una sola volta (le funzioni Plot* lo chiamano da sole se i dati non ci sono)
bool LoadData() {
    if (d.nEv > 0) return true;               // gia' caricato
    TFile *f = TFile::Open("data.root", "READ");
    if (!f || f->IsZombie()) { std::cerr << "Impossibile aprire data.root" << std::endl; return false; }
    bool ok = d.Load(f);
    f->Close();   // da qui in poi si lavora solo sui dati in memoria
    return ok;
}

// =====================================================================
// FUNZIONI DI SUPPORTO (istogrammi per tutti i 16 pixel)
// =====================================================================

// crea un istogramma per ogni pixel valido: signal [0,100] per i centrali, p2p [0,8.5] per gli esterni
void MakeHists16(const char *name, const char *title, TH1D *h[NPIX]) {
    for (int i = 0; i < NPIX; i++) {
        h[i] = nullptr;
        if (!d.px[i].Valid()) continue;
        if (d.px[i].central)
            h[i] = new TH1D(Form("%s_pixel%d", name, i),
                            Form("Pixel %d (%s);signal [SCOPE];Eventi", i, title), 100, 0, 100);
        else
            h[i] = new TH1D(Form("%s_pixel%d", name, i),
                            Form("Pixel %d (%s);P2P [ADC];Eventi", i, title), 40, 0, 8.5);
        h[i]->SetDirectory(nullptr);
    }
}

// riempie gli istogrammi con i valori dell'evento ev:
// esterni sempre (p2p), centrali solo se accesi (signal)
void FillEvent16(Long64_t ev, TH1D *h[NPIX]) {
    for (int i = 0; i < NPIX; i++) {
        if (!h[i]) continue;
        if (d.px[i].central && !d.px[i].IsOn(ev)) continue;
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
// ANALISI 1: mappa 2x2 dei pixel centrali
// =====================================================================
void PlotMap() {
    if (!LoadData()) return;
    TCanvas *c = new TCanvas("c1", "Mappa Pixel 2x2", 800, 600);
    TH2D *h = new TH2D("h_map", "Stato Pixel (2x2);Asse X;Asse Y", 2, 0, 2, 2, 0, 2);
    h->SetDirectory(nullptr);
    gStyle->SetOptStat(0);

    for (Long64_t ev = 0; ev < d.nEv; ev++)
        for (int i = 0; i < NCEN; i++)
            if (d.on[ev][i]) h->Fill(i % 2 + 0.5, (1 - i / 2) + 0.5);

    h->GetXaxis()->SetNdivisions(2);
    h->GetYaxis()->SetNdivisions(2);
    h->SetMinimum(0);
    h->SetMarkerSize(1.8);
    c->cd();
    h->Draw("COLZ TEXT");

    for (int i = 0; i < NCEN; i++) {
        int x = i % 2, y = 1 - (i / 2);
        TLatex *label = new TLatex(x + 0.5, y + 0.7, Form("Pixel %d", CEN_PX[i]));
        label->SetTextAlign(22);
        label->SetTextSize(0.03);
        label->SetTextColor(kBlack);
        label->Draw();
    }
    c->Modified();
    c->Update();
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
    TH1D *h_cluster = new TH1D("h_cluster", "Numero di pixel accesi per evento;Numero di pixel accesi;Numero di eventi",
                               5, -0.5, 4.5);
    h_cluster->SetDirectory(nullptr);
    for (int k = 0; k <= NCEN; k++) h_cluster->GetXaxis()->SetBinLabel(k + 1, Form("%d", k));
    for (Long64_t ev = 0; ev < d.nEv; ev++) h_cluster->Fill(d.csize[ev]);
    c2->cd();
    h_cluster->Draw();
    c2->Modified();
    c2->Update();

    // ---- segnale totale per cluster size (4 pad) ----
    TCanvas *c4 = new TCanvas("c4", "Signal vs Cluster Size (1D)", 1200, 900);
    c4->Divide(2, 2);
    TH1D *hs[NCEN];
    for (int k = 0; k < NCEN; k++) {
        hs[k] = new TH1D(Form("h_signal_cluster_%d", k + 1),
                         Form("Segnale per cluster size = %d;Segnale [mV];Numero di eventi", k + 1), 100, 0, 100);
        hs[k]->SetDirectory(nullptr);
    }
    for (Long64_t ev = 0; ev < d.nEv; ev++)
        if (d.csize[ev] >= 1) hs[d.csize[ev] - 1]->Fill(d.totSig[ev]);
    for (int k = 0; k < NCEN; k++) { c4->cd(k + 1); hs[k]->Draw(); }
    c4->Modified();
    c4->Update();

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
    hn[0]->SetTitle("Segnale totale normalizzato per cluster size;Segnale [mV];Frazione di eventi");
    hn[0]->SetMaximum(1.15 * ymax);
    hn[0]->SetMinimum(0.);
    c8->cd();
    hn[0]->Draw("HIST");
    for (int k = 1; k < NCEN; k++) hn[k]->Draw("HIST SAME");

    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    for (int k = 0; k < NCEN; k++) leg->AddEntry(hn[k], Form("%d Pixel%s", k + 1, k ? "s" : ""), "l");
    leg->Draw();
    c8->Modified();
    c8->Update();
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
                             Form("P2P / Signal - %d pixel interni accesi", s + 1), 1200, 800);
        cCS[s]->Divide(4, 4);
        MakeHists16(Form("h_cs_s%d", s + 1), Form("cluster size = %d", s + 1), h[s]);
    }
    for (Long64_t ev = 0; ev < d.nEv; ev++) {
        int cs = d.csize[ev];
        if (cs < 1 || cs > NCEN) continue;
        FillEvent16(ev, h[cs - 1]);
    }
    for (int s = 0; s < NCEN; s++) Draw16(cCS[s], h[s]);
}

// =====================================================================
// ANALISI 4: cluster size = 1, per ogni pixel centrale come seed.
// Taglio opzionale: signal del pixel acceso > sigCut
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
// ANALISI 5: eventi in cui il p2p di un pixel ESTERNO (selPx) e' in [p2pMin, p2pMax]
// =====================================================================
void PlotP2PSelection(int selPx = 0, double p2pMin = 1.0, double p2pMax = 3.0) {
    if (!LoadData()) return;

    if (selPx < 0 || selPx >= NPIX || IsCentralPx(selPx) || !d.px[selPx].Valid()) {
        std::cerr << "selPx deve essere un pixel esterno valido: 0,1,2,3,4,7,8,11,12,13,14,15" << std::endl;
        return;
    }

    TString selTitle = Form("p2p pixel %d in [%g, %g]", selPx, p2pMin, p2pMax);
    TString selName  = Form("_sel%d_p2p_%g_%g", selPx, p2pMin, p2pMax);
    selName.ReplaceAll(".", "p");
    selName.ReplaceAll("-", "m");

    TCanvas *c = new TCanvas(Form("c6%s", selName.Data()), Form("Eventi con %s", selTitle.Data()), 1200, 900);
    c->Divide(4, 4);
    TH1D *h[NPIX];
    MakeHists16(Form("h_sel%s", selName.Data()), selTitle.Data(), h);

    int nSel = 0;
    for (Long64_t ev = 0; ev < d.nEv; ev++) {
        double p = d.px[selPx].p2p[ev];
        if (p < p2pMin || p > p2pMax) continue;
        nSel++;
        FillEvent16(ev, h);
    }
    std::cout << "Eventi con " << selTitle << ": " << nSel << std::endl;
    Draw16(c, h);
}

// =====================================================================
// SCRITTURA DEL TREE "cluster" in data.root (dai dati gia' in memoria)
// =====================================================================
void WriteClusterTree() {
    if (!LoadData()) return;
    TFile *f = TFile::Open("data.root", "UPDATE");
    if (!f || f->IsZombie()) { std::cerr << "Impossibile aprire data.root in UPDATE" << std::endl; return; }
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
// =====================================================================
// COSTANTI DA AGGIUNGERE IN ALTO (dopo le altre costanti / opzioni)
// =====================================================================
const double PITCH = 10.;                   // [micron] dimensione del pixel 20x20
const int    centroidSelPx = 0;             // pixel esterno da studiare

// =====================================================================
// ANALISI 6: segnale dei pixel interni quando si considera un evento
//   del pixel esterno selPx, in funzione della distanza (geometria 4x4)
//   Eventi selezionati: p2p di selPx in [p2pMin, p2pMax]
//   Tutto in un unico canvas: griglia 4x4 (pad = posizione del pixel)
//   + due riepiloghi a destra. Distanze in micron.
// =====================================================================
inline double PxX(int id) { return (id % 4) * PITCH; }   // colonna [um]
inline double PxY(int id) { return (id / 4) * PITCH; }   // riga [um]

// distanza tra due pixel [um]
inline double PxDist(int a, int b) {
    double dx = PxX(a) - PxX(b), dy = PxY(a) - PxY(b);
    return std::sqrt(dx*dx + dy*dy);
}

// calcola il centroide interno dell'evento [um]; false se non definito
bool Centroid(Long64_t ev, double &cx, double &cy) {
    if (d.csize[ev] < 1 || d.totSig[ev] <= 0.) return false;
    cx = 0.; cy = 0.;
    for (int i = 0; i < NCEN; i++) {
        if (!d.on[ev][i]) continue;
        double s = d.px[CEN_PX[i]].signal[ev];
        cx += s * PxX(CEN_PX[i]);
        cy += s * PxY(CEN_PX[i]);
    }
    cx /= d.totSig[ev];
    cy /= d.totSig[ev];
    return true;
}

void PlotCentroid(int selPx = 0, double p2pMin = 0., double p2pMax = 1e9) {
    if (!LoadData()) return;

    if (selPx < 0 || selPx >= NPIX || IsCentralPx(selPx) || !d.px[selPx].Valid()) {
        std::cerr << "selPx deve essere un pixel esterno valido: 0,1,2,3,4,7,8,11,12,13,14,15" << std::endl;
        return;
    }

    const double dMax   = 4.5 * PITCH;   // 0 - 90 um
    const double sMax   = 300.;          // range segnale totale
    const double sMaxPx = 100.;          // range segnale singolo pixel

    // ---- istogrammi ----
    // segnale totale interno vs distanza centroide - selPx
    TH2D *hTot = new TH2D(Form("h_intTot_sel%d", selPx),
                          Form("Segnale totale interni vs distanza dal pixel %d;Distanza centroide - pixel %d [#mum];Segnale totale [SCOPE]",
                               selPx, selPx),
                          45, 0, dMax, 100, 0, sMax);
    hTot->SetDirectory(nullptr);

    // segnale di ogni pixel centrale vs distanza (fissa) pixel centrale - selPx
    TH2D *hAll = new TH2D(Form("h_intAll_sel%d", selPx),
                          Form("Segnale pixel interni vs distanza dal pixel %d;Distanza pixel interno - pixel %d [#mum];Signal [SCOPE]",
                               selPx, selPx),
                          45, 0, dMax, 100, 0, sMaxPx);
    hAll->SetDirectory(nullptr);

    // un istogramma per ogni pixel centrale: signal vs distanza centroide - selPx
    TH2D *hEach[NCEN];
    for (int k = 0; k < NCEN; k++) {
        hEach[k] = new TH2D(Form("h_intPx%d_sel%d", CEN_PX[k], selPx),
                            Form("Pixel %d (d=%.0f #mum);Dist. centroide-px %d [#mum];Signal [SCOPE]",
                                 CEN_PX[k], PxDist(CEN_PX[k], selPx), selPx),
                            30, 0, dMax, 50, 0, sMaxPx);
        hEach[k]->SetDirectory(nullptr);
    }

    // p2p del pixel selezionato (eventi scelti)
    TH1D *hSel = new TH1D(Form("h_selP2P_%d", selPx),
                          Form("Pixel %d: P2P;P2P [ADC];Eventi", selPx), 40, 0, 8.5);
    hSel->SetDirectory(nullptr);

    // ---- riempimento ----
    int nUsed = 0;
    for (Long64_t ev = 0; ev < d.nEv; ev++) {
        double p0 = d.px[selPx].p2p[ev];
        if (p0 < p2pMin || p0 > p2pMax) continue;

        double cx, cy;
        if (!Centroid(ev, cx, cy)) continue;
        nUsed++;

        double dx = cx - PxX(selPx);
        double dy = cy - PxY(selPx);
        double distC = std::sqrt(dx*dx + dy*dy);

        hSel->Fill(p0);
        hTot->Fill(distC, d.totSig[ev]);

        for (int k = 0; k < NCEN; k++) {
            if (!d.on[ev][k]) continue;
            double s = d.px[CEN_PX[k]].signal[ev];
            hAll->Fill(PxDist(CEN_PX[k], selPx), s);
            hEach[k]->Fill(distC, s);
        }
    }
    std::cout << "Eventi con p2p pixel " << selPx << " in [" << p2pMin << ", " << p2pMax
              << "] e centroide definito: " << nUsed << std::endl;

    // ---- canvas unico ----
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.07);
    TCanvas *c = new TCanvas(Form("c_intSig_sel%d", selPx),
                             Form("Segnale pixel interni - eventi pixel %d", selPx), 1600, 900);

    const double gridW = 0.70;           // frazione di larghezza per la griglia 4x4
    const double w = gridW / 4., h = 1. / 4.;

    for (int id = 0; id < NPIX; id++) {
        int r = id / 4, col = id % 4;
        double x1 = col * w, x2 = (col + 1) * w;
        double y2 = 1. - r * h, y1 = y2 - h;

        c->cd();
        TPad *pad = new TPad(Form("pad_px%d_sel%d", id, selPx), "", x1, y1, x2, y2);
        pad->SetMargin(0.17, 0.14, 0.17, 0.12);
        pad->Draw();
        pad->cd();

        int k = CenIndex(id);
        if (k >= 0) {
            // pixel centrale: signal vs distanza centroide
            hEach[k]->GetXaxis()->SetLabelSize(0.06);
            hEach[k]->GetYaxis()->SetLabelSize(0.06);
            hEach[k]->GetXaxis()->SetTitleSize(0.06);
            hEach[k]->GetYaxis()->SetTitleSize(0.06);
            hEach[k]->GetXaxis()->SetTitleOffset(1.0);
            hEach[k]->GetYaxis()->SetTitleOffset(1.2);
            hEach[k]->Draw("COLZ");
        } else if (id == selPx) {
            // pixel selezionato: distribuzione p2p
            pad->SetFillColor(kYellow - 9);
            hSel->GetXaxis()->SetLabelSize(0.06);
            hSel->GetYaxis()->SetLabelSize(0.06);
            hSel->GetXaxis()->SetTitleSize(0.06);
            hSel->GetYaxis()->SetTitleSize(0.06);
            hSel->GetYaxis()->SetTitleOffset(1.2);
            hSel->Draw();
        } else {
            // altri esterni: solo etichetta
            pad->SetFillColor(kGray);
            TLatex t;
            t.SetNDC();
            t.SetTextAlign(22);
            t.SetTextSize(0.12);
            t.DrawLatex(0.5, 0.58, Form("Pixel %d", id));
            t.SetTextSize(0.08);
            t.DrawLatex(0.5, 0.38, Form("d = %.0f #mum da px %d", PxDist(id, selPx), selPx));
        }
    }

    // ---- riepiloghi a destra ----
    c->cd();
    TPad *pR1 = new TPad(Form("pad_tot_sel%d", selPx), "", gridW, 0.5, 1., 1.);
    pR1->SetMargin(0.14, 0.14, 0.12, 0.10);
    pR1->Draw();
    pR1->cd();
    hTot->Draw("COLZ");
    TProfile *pTot = hTot->ProfileX(Form("p_intTot_sel%d", selPx));
    pTot->SetLineColor(kRed);
    pTot->SetMarkerColor(kRed);
    pTot->SetMarkerStyle(20);
    pTot->SetMarkerSize(0.7);
    pTot->Draw("SAME");

    c->cd();
    TPad *pR2 = new TPad(Form("pad_all_sel%d", selPx), "", gridW, 0., 1., 0.5);
    pR2->SetMargin(0.14, 0.14, 0.12, 0.10);
    pR2->Draw();
    pR2->cd();
    hAll->Draw("COLZ");

    c->Modified();
    c->Update();
}

// =====================================================================
// ANALISI 7: eventi con i due pixel centrali OPPOSTI (diagonale) accesi
//   coppia A: (5,10)   coppia B: (6,9)
//   exclusive = true  -> solo quei due pixel accesi (cluster size = 2),
//                        gli altri due centrali spenti
//   exclusive = false -> basta che i due pixel della coppia siano accesi
//                        (gli altri centrali possono essere accesi o no)
//   Un canvas 4x4 per coppia, con signal (centrali) / p2p (esterni)
//   di tutti i 16 pixel.
// =====================================================================
void PlotOppositePairs(bool exclusive = true) {
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
        std::cout << "Eventi con pixel " << CEN_PX[pairs[p][0]] << " e " << CEN_PX[pairs[p][1]]
                  << " accesi" << (exclusive ? " (solo questi)" : "") << ": " << nSel[p] << std::endl;
        Draw16(c[p], h[p]);
    }
}

// =====================================================================
// MACRO PRINCIPALE: esegue tutte le analisi con gli argomenti di default
// (le singole funzioni si possono richiamare da ROOT anche con altri argomenti)
// =====================================================================
void cls_analysis() {
    if (!LoadData()) return;

    if (optWriteClusterTree) WriteClusterTree();

    if (optDoMap)          PlotMap();
    if (optDoClusterCount) PlotClusterCount();
    if (optDoPerClusterSz) PlotPerClusterSize();
    if (optDoSeed)         PlotSeed(true, 50., 150);    //bool applyCut = true, double sigCut_min = 0., double sigCut_max = 150.
    if (optDoP2PSelection) PlotP2PSelection(0, 5.5, 6); // pixel esterno 4, range [2, 5]
    if (optDoCentroid)     PlotCentroid(centroidSelPx, 0., 8.5);   // eventi con p2p pixel 0 in range [xmin, xmax]
    if (optDoOpposite)     PlotOppositePairs(true);   // true = solo i due pixel opposti accesi

}