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

const bool   optWriteClusterTree = false;   // true: scrive/aggiorna il tree "cluster" in data.root

// quali analisi eseguire
const bool   optDoMap          = true;
const bool   optDoClusterCount = true;
const bool   optDoPerClusterSz = true;
const bool   optDoSeed         = true;
const bool   optDoP2PSelection = true;

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
void PlotSeed(bool applyCut = true, double sigCut = 60.) {
    if (!LoadData()) return;

    TString cutTitle = "", cutName = "";
    if (applyCut) {
        cutTitle = Form(", signal > %.0f", sigCut);
        cutName  = Form("_gt%.0f", sigCut);
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
        if (applyCut && d.px[CEN_PX[s]].signal[ev] <= sigCut) continue;
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
// MACRO PRINCIPALE: esegue tutte le analisi con gli argomenti di default
// (le singole funzioni si possono richiamare da ROOT anche con altri argomenti)
// =====================================================================
void cluster() {
    if (!LoadData()) return;

    if (optWriteClusterTree) WriteClusterTree();

    if (optDoMap)          PlotMap();
    if (optDoClusterCount) PlotClusterCount();
    if (optDoPerClusterSz) PlotPerClusterSize();
    if (optDoSeed)         PlotSeed();
    if (optDoP2PSelection) PlotP2PSelection();
}