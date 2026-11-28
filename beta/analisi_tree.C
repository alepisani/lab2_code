#include <iostream>
#include <fstream>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TF1.h>
#include <TLine.h>

using namespace std;

#define TIME_CUT true

void analisi_tree()
{
    double pmax_thr = 12;


    gStyle->SetOptStat(0);
    int biasVoltage = 600;
    int Nsigma = 3;

    int nBins_pmax = 801;
    double binmax_pmax = 200.5;

    string inputFolder = "FBK/Run4_PRE_LF-FBK_QC-TS_W16_S52_LGAD-B/";
    string inputName   = "stats_Sr_Run4_Ch0-" + to_string(biasVoltage) + "V_trig2750V.root";

    TFile *inputFile = new TFile(("Misure/" + inputFolder + inputName).c_str(), "READ");
    TTree *tree = (TTree*)inputFile->Get("Analysis");
    
    vector<float> *tree_pmax=nullptr, *tree_tmax=nullptr, *tree_area_new=nullptr, *tree_dvdt=nullptr, *tree_rms=nullptr;
    vector<vector<float>> *tree_cfd=nullptr;
    tree->SetBranchAddress("pmax", &tree_pmax);
    tree->SetBranchAddress("tmax", &tree_tmax);
    tree->SetBranchAddress("area_new", &tree_area_new);
    tree->SetBranchAddress("dvdt", &tree_dvdt);
    tree->SetBranchAddress("rms", &tree_rms);
    tree->SetBranchAddress("cfd", &tree_cfd);

    TH1F *h_pmax = new TH1F("h_pmax", "Distribuzione dell'ampiezza massima;Ampiezza [mV];Conteggi", nBins_pmax, -0.5, binmax_pmax);
    TH1F *h_tmax = new TH1F("h_tmax", "Tempo del massimo del segnale;t_{max} [ns];Conteggi", 101, -1, 4);

    vector<float> pmax, tmax, area_new, dvdt, cfd20_DUT, cfd70_DUT, cfd80_DUT, cfd20_MCP, rms;

    for(int i=0; i<tree->GetEntries(); i++)
    {
        tree->GetEntry(i);
        h_pmax->Fill((*tree_pmax)[0]);
        h_tmax->Fill((*tree_tmax)[0]);
        pmax.push_back((*tree_pmax)[0]);
        tmax.push_back((*tree_tmax)[0]);
        rms.push_back((*tree_rms)[0]);
        area_new.push_back((*tree_area_new)[0]);
        dvdt.push_back((*tree_dvdt)[0]);
        cfd20_DUT.push_back((*tree_cfd)[0][1]);
        cfd70_DUT.push_back((*tree_cfd)[0][6]);
        cfd80_DUT.push_back((*tree_cfd)[0][7]);
        cfd20_MCP.push_back((*tree_cfd)[1][1]);
    }

    TCanvas *c1 = new TCanvas("c1", "Peak Amplitude Distribution", 800, 1200);
    TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    c1->Divide(1, 2);
    c1->cd(1);
    leg->AddEntry(h_pmax, "Scartati (tmax)", "f");
    h_pmax->Draw();
    gPad->SetLogy();
    h_pmax->SetLineColor(kRed);
    h_pmax->SetLineWidth(2);
    h_pmax->SetFillColor(kRed);
    c1->cd(2);
    h_tmax->SetLineColor(kRed);
    h_tmax->SetLineWidth(2);
    h_tmax->SetFillColor(kRed);
    h_tmax->Draw();

    TF1 *fit_tmax = new TF1("fit_tmax", "gaus(0)+pol0(3)", 0, 2.2);
    fit_tmax->SetParameters(300,0.8,0.1,40);
    fit_tmax->SetParLimits(2, 0, 100);
    fit_tmax->SetNpx(400);
    fit_tmax->SetLineColor(kBlack);
    fit_tmax->SetLineWidth(3);
    c1->cd(2);
    h_tmax->Fit(fit_tmax, "RQ");
    cout << "Norm: " << fit_tmax->GetParameter(0) << endl;
    cout << "Mu: " << fit_tmax->GetParameter(1) << endl;
    cout << "Sigma: " << fit_tmax->GetParameter(2) << endl;
    cout << "Fit result: p-value = " << fit_tmax->GetProb() << endl;
    double tmean = fit_tmax->GetParameter(1);
    double tsigma = fit_tmax->GetParameter(2);

    vector<float> pmax_tclean, tmax_tclean, area_new_tclean, dvdt_tclean, cfd20_DUT_tclean, cfd70_DUT_tclean, cfd80_DUT_tclean, cfd20_MCP_tclean, rms_tclean;
    for(int i=0; i<pmax.size(); i++)
    {
        if((TIME_CUT && tmax[i]>tmean-Nsigma*tsigma && tmax[i]<tmean+Nsigma*tsigma) || !TIME_CUT)
        {
            pmax_tclean.push_back(pmax[i]);
            tmax_tclean.push_back(tmax[i]);
            area_new_tclean.push_back(area_new[i]);
            dvdt_tclean.push_back(dvdt[i]);
            cfd20_DUT_tclean.push_back(cfd20_DUT[i]);
            cfd70_DUT_tclean.push_back(cfd70_DUT[i]);
            cfd80_DUT_tclean.push_back(cfd80_DUT[i]);
            cfd20_MCP_tclean.push_back(cfd20_MCP[i]);
            rms_tclean.push_back(rms[i]);
        }
    }
    TH1F *h_pmax_tclean = new TH1F("h_pmax_tclean", "Distribuzione dell'ampiezza massima;Ampiezza [mV];Conteggi", nBins_pmax, -0.5, binmax_pmax);
    for(int i=0; i<pmax_tclean.size(); i++)
        h_pmax_tclean->Fill(pmax_tclean[i]);
    c1->cd(1);
    leg->AddEntry(h_pmax_tclean, "Scartati (pmax)", "f");
    h_pmax_tclean->SetLineColor(kAzure+7);
    h_pmax_tclean->SetLineWidth(2);
    h_pmax_tclean->SetFillColor(kAzure+7);
    h_pmax_tclean->Draw("SAME");
    leg->Draw();
    leg->SetTextSize(0.05);

    TH1F *h_tmax_tclean = new TH1F("h_tmax_tclean", "Tempo del massimo del segnale;t_{max} [ns];Conteggi", 101, -1, 4);
    for(int i=0; i<tmax_tclean.size(); i++)
        h_tmax_tclean->Fill(tmax_tclean[i]);
    c1->cd(2);
    h_tmax_tclean->SetLineColor(kAzure+7);
    h_tmax_tclean->SetLineWidth(2);
    h_tmax_tclean->SetFillColor(kAzure+7);
    h_tmax_tclean->Draw("SAME");

    vector<float> pmax_clean, tmax_clean, area_new_clean, dvdt_clean, cfd20_DUT_clean, cfd70_DUT_clean, cfd80_DUT_clean, cfd20_MCP_clean, rms_clean;
    for(int i=0; i<pmax_tclean.size(); i++)
    {
        if(pmax_tclean[i]>pmax_thr)
        {
            pmax_clean.push_back(pmax_tclean[i]);
            tmax_clean.push_back(tmax_tclean[i]);
            area_new_clean.push_back(area_new_tclean[i]);
            dvdt_clean.push_back(dvdt_tclean[i]);
            cfd20_DUT_clean.push_back(cfd20_DUT_tclean[i]);
            cfd70_DUT_clean.push_back(cfd70_DUT_tclean[i]);
            cfd80_DUT_clean.push_back(cfd80_DUT_tclean[i]);
            cfd20_MCP_clean.push_back(cfd20_MCP_tclean[i]);
            rms_clean.push_back(rms_tclean[i]);
        }
    }
    TH1F *h_pmax_clean = new TH1F("h_pmax_clean", "Distribuzione dell'ampiezza massima;Ampiezza [mV];Conteggi", nBins_pmax, -0.5, binmax_pmax);
    for(int i=0; i<pmax_clean.size(); i++)
        h_pmax_clean->Fill(pmax_clean[i]);
    c1->cd(1);
    leg->AddEntry(h_pmax_clean, "Dati finali", "f");
    h_pmax_clean->SetLineColor(kGreen+2);
    h_pmax_clean->SetLineWidth(2);
    h_pmax_clean->SetFillColor(kGreen+2);
    h_pmax_clean->Draw("SAME");
    leg->Draw();
    leg->SetTextSize(0.05);

    TH1F *h_tmax_clean = new TH1F("h_tmax_clean", "Tempo del massimo del segnale;t_{max} [ns];Conteggi", 101, -1, 4);
    for(int i=0; i<tmax_clean.size(); i++)
        h_tmax_clean->Fill(tmax_clean[i]);
    c1->cd(2);
    h_tmax_clean->SetLineColor(kGreen+2);
    h_tmax_clean->SetLineWidth(2);
    h_tmax_clean->SetFillColor(kGreen+2);
    h_tmax_clean->Draw("SAME");

    c1->cd(1);
    TLine *line2 = new TLine(pmax_thr, h_pmax->GetYaxis()->GetBinLowEdge(1), pmax_thr, h_pmax->GetMaximum());
    line2->SetLineColor(kBlack);
    line2->SetLineWidth(2);
    line2->SetLineStyle(2);
    line2->Draw("SAME");

    c1->cd(2);
    fit_tmax->Draw("SAME");

    TFile *outputFile = new TFile(("Misure_clean/" + inputFolder + inputName).c_str(), "RECREATE");
    TTree *outputTree = new TTree("Analysis", "Analysis");

    double out_pmax_clean, out_tmax_clean, out_area_new_clean, out_dvdt_clean, out_cfd20_DUT_clean, out_cfd70_DUT_clean, out_cfd80_DUT_clean, out_cfd20_MCP_clean, out_rms_clean;
    double out_mutime, out_sigmatime, out_pthr;
    bool time_cut;

    TBranch *b_pmax_clean = outputTree->Branch("pmax", &out_pmax_clean, "pmax_clean/D");
    TBranch *b_tmax_clean = outputTree->Branch("tmax", &out_tmax_clean, "tmax_clean/D");
    TBranch *b_area_new_clean = outputTree->Branch("area_new", &out_area_new_clean, "area_new_clean/D");
    TBranch *b_dvdt_clean = outputTree->Branch("dvdt", &out_dvdt_clean, "dvdt_clean/D");
    TBranch *b_cfd20_DUT_clean = outputTree->Branch("cfd20_DUT", &out_cfd20_DUT_clean, "cfd20_DUT_clean/D");
    TBranch *b_cfd70_DUT_clean = outputTree->Branch("cfd70_DUT", &out_cfd70_DUT_clean, "cfd70_DUT_clean/D");
    TBranch *b_cfd80_DUT_clean = outputTree->Branch("cfd80_DUT", &out_cfd80_DUT_clean, "cfd80_DUT_clean/D");
    TBranch *b_cfd20_MCP_clean = outputTree->Branch("cfd20_MCP", &out_cfd20_MCP_clean, "cfd20_MCP_clean/D");
    TBranch *b_rms_clean = outputTree->Branch("rms", &out_rms_clean, "rms_clean/D");
    TBranch *b_mutime = outputTree->Branch("mutime", &out_mutime, "mutime/D");
    TBranch *b_sigmatime = outputTree->Branch("sigmatime", &out_sigmatime, "sigmatime/D");
    TBranch *b_pthr = outputTree->Branch("pthr", &out_pthr, "pthr/D");
    TBranch *b_time_cut = outputTree->Branch("time_cut", &time_cut, "time_cut/O");

    for(int i=0; i<pmax_clean.size(); i++)
    {
        out_pmax_clean = pmax_clean[i];
        out_tmax_clean = tmax_clean[i];
        out_area_new_clean = area_new_clean[i];
        out_dvdt_clean = dvdt_clean[i];
        out_cfd20_DUT_clean = cfd20_DUT_clean[i];
        out_cfd70_DUT_clean = cfd70_DUT_clean[i];
        out_cfd80_DUT_clean = cfd80_DUT_clean[i];
        out_cfd20_MCP_clean = cfd20_MCP_clean[i];
        out_rms_clean = rms_clean[i];
        out_mutime = tmean;
        out_sigmatime = tsigma;
        out_pthr = pmax_thr;
        time_cut = TIME_CUT;
        outputTree->Fill();
    }
    outputTree->Write();
    outputFile->Close();
}
