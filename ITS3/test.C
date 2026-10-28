#include "TFile.h"
#include "TTree.h"
#include <iostream>

void test() {
    TFile *f = TFile::Open("GRUPPO1_55Fe/data.root");
    TTree *t = (TTree*)f->Get("pixel5");
    double t0, t30, t50;
    t->SetBranchAddress("CFD_times.t0",  &t0);
    t->SetBranchAddress("CFD_times.t30", &t30);
    t->SetBranchAddress("CFD_times.t50", &t50);

    for (Long64_t i = 0; i < 5; i++) {
        t->GetEntry(i);
        std::cout << t0 << " " << t30 << " " << t50 << std::endl;
    }
}