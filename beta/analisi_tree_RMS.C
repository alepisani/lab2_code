#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TInterpreter.h>

using namespace std;

#define TIME_CUT false

void analisi_tree_RMS()
{
    int biasVoltage = 200;

    string inputFolder = "HPK/Batch2/Run3_PRE_HPK_QC-TS_W81_S2_LGAD-B/";
    string inputName   = "stats_Sr_Run3_Ch0-" + to_string(biasVoltage) + "V_trig2750V.root";


    TFile *inputFile2 = new TFile(("Misure_temp/" + inputFolder + inputName).c_str(), "READ");
    TTree *tree2 = (TTree*)inputFile2->Get("Analysis");
    double pmax_thr, tmean, tsigma;

    tree2->SetBranchAddress("pthr", &pmax_thr);
    tree2->SetBranchAddress("mutime", &tmean);
    tree2->SetBranchAddress("sigmatime", &tsigma);

    tree2->GetEntry(0);

    double tmin = tmean - 3*tsigma;
    double tmax = tmean + 3*tsigma;

    inputFile2->Close();

    TFile *inputFile1 = new TFile(("Misure/" + inputFolder + inputName).c_str(), "READ");
    TTree *tree1 = (TTree*)inputFile1->Get("Analysis");

    gInterpreter->GenerateDictionary("vector<vector<double> >", "vector");
    
    vector<float> *tree_pmax=nullptr, *tree_tmax=nullptr, *tree_area_new=nullptr, *tree_dvdt=nullptr, *tree_rms=nullptr;
    vector<vector<double>> *tree_cfd=nullptr;
    tree1->SetBranchAddress("pmax", &tree_pmax);
    tree1->SetBranchAddress("tmax", &tree_tmax);
    tree1->SetBranchAddress("area_new", &tree_area_new);
    tree1->SetBranchAddress("dvdt", &tree_dvdt);
    tree1->SetBranchAddress("rms", &tree_rms);
    tree1->SetBranchAddress("cfd", &tree_cfd);
    
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

    for(int i=0; i<tree1->GetEntries(); i++)
    {
        tree1->GetEntry(i);
        if(TIME_CUT && ((*tree_tmax)[0] <= tmin || (*tree_tmax)[0] >= tmax))
            continue;
        else if((*tree_pmax)[0] > pmax_thr)
        {
            out_pmax_clean = (*tree_pmax)[0];
            out_tmax_clean = (*tree_tmax)[0];
            out_area_new_clean = (*tree_area_new)[0];
            out_dvdt_clean = (*tree_dvdt)[0];
            out_cfd20_DUT_clean = (*tree_cfd)[0][1];
            out_cfd70_DUT_clean = (*tree_cfd)[0][6];
            out_cfd80_DUT_clean = (*tree_cfd)[0][7];
            out_cfd20_MCP_clean = (*tree_cfd)[1][1];
            out_rms_clean = (*tree_rms)[0];
            out_mutime = tmean;
            out_sigmatime = tsigma;
            out_pthr = pmax_thr;
            time_cut = TIME_CUT;
            outputTree->Fill();
        }

    }
    outputTree->Write();
    outputFile->Close();

}



void test()
{
    TFile *f = TFile::Open("Misure/FBK/Run4_PRE_LF-FBK_QC-TS_W16_S52_LGAD-B/stats_Sr_Run4_Ch0-120V_trig2750V.root");
    TTree *t = (TTree*)f->Get("Analysis");
    t->Print();
    cout << t->GetBranch("cfd")->GetClassName() << endl;
    f->ShowStreamerInfo();

    vector<vector<float>> *cfd = nullptr;
    t->SetBranchAddress("cfd", &cfd);
    t->GetEntry(0);
    cout << cfd->size() << " " << (*cfd)[0].size() << endl;
}