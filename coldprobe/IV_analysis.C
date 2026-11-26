#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TF1.h>

using namespace std;

#define Kb 8.617333262145e-5 // eV/K

int find_100V(const vector<double>&);

void IV_analysis()
{
    double err_temperature = 0.2;

    int ordJ[] = {8,6,4,2,0,1,3,5,7}; //Mi dispiace

    string filename = "IV_data.root";
    string tree_names[4] = {"New", "1e14", "8e14", "15e14"};
    string fluences[4] = {"0.0", "1e14", "8e14", "15e14"};
    TFile *file = new TFile(filename.c_str(), "READ");

    int colors[] = {kBlack, kRed, kBlue, kGreen+2, kMagenta, kCyan+2, kOrange+7, kViolet+1, kTeal+1, kAzure+1};
    
    vector<double> *V=nullptr, *I=nullptr, *V_err=nullptr, *I_err=nullptr;

    vector<vector<double>> vI100(4), err_vI100(4), vT100(4);
    double T;

    TCanvas *c1 = new TCanvas("c1", "IV Curves", 2000, 1000);
    c1->Divide(2, 2);

    TCanvas *c2 = new TCanvas("c2", "Arrhenius Plot", 1800, 1400);
    c2->Divide(2, 2);

    TCanvas *c3 = new TCanvas("c3", "Arrhenius Plot (combined)", 900, 600);
    TMultiGraph *mg_arr = new TMultiGraph();
    TLegend *legend_arr = new TLegend(0.75, 0.7, 0.95, 0.9);

    for(int i=0;i<4;i++)
    {
        c1->cd(i+1);
        gPad->SetLeftMargin(0.1);
        gPad->SetBottomMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.10);
        TMultiGraph *mg = new TMultiGraph();
        TLegend *legend = new TLegend(0.7, 0.45, 0.9, 0.65);
        legend->SetName(Form("legend_%s", tree_names[i].c_str()));
        legend->SetTextSize(0.03);

        TTree *tree = (TTree*)file->Get(Form("IV_PiN_%s", tree_names[i].c_str()));
        tree->SetBranchAddress("T", &T);
        tree->SetBranchAddress("V", &V);
        tree->SetBranchAddress("I", &I);
        tree->SetBranchAddress("V_err", &V_err);
        tree->SetBranchAddress("I_err", &I_err);

        int counter = 0;

        vector<double> overT, err_overT, logIoverT2, err_logIoverT2;
        for(int j=0;j<tree->GetEntries();j++)
        {
            tree->GetEntry(ordJ[j]);
            double T_K = (double)T + 273.15;
            overT.push_back(1000./T_K);
            err_overT.push_back(err_temperature/(T_K*T_K));
            int index_100V = find_100V(*V);
            double I100 = (*I)[index_100V];
            double err_I100 = (*I_err)[index_100V];
            logIoverT2.push_back(log(I100) - 2.*log(T_K));

            vI100[i].push_back(I100);
            err_vI100[i].push_back(err_I100);
            vT100[i].push_back(T);

            double err = sqrt(pow(err_I100/I100, 2) + pow(2.*err_temperature/T_K, 2));
            err_logIoverT2.push_back(err);

            if(fabs(T) == 0 || fabs(T) == 10 || fabs(T) == 25)      //così non abbiamo troppe curve
            {
                TGraphErrors *gr = new TGraphErrors(V->size(), V->data(), I->data(), V_err->data(), I_err->data());
                gr->SetName(Form("IV_PiN_%s_T%.1fC", tree_names[i].c_str(), T));
                gr->SetTitle(Form("T = %.1f C", T));
                gr->SetMarkerStyle(20);
                gr->SetMarkerSize(0.2);
                gr->SetLineColor(colors[counter]);
                gr->SetMarkerColor(colors[counter]);

                mg->Add(gr, "P");
                legend->AddEntry(gr, Form("T = %.1f C", T), "pl");
                counter++;
            }
        }
        mg->SetTitle(Form("IV Curves for #Phi = %s n_{eq}/cm^{2}", fluences[i].c_str()));
        mg->Draw("AP");
        mg->GetXaxis()->SetTitle("V [V]");
        mg->GetYaxis()->SetTitle("I [nA]");
        mg->GetXaxis()->SetTitleSize(0.05);
        mg->GetYaxis()->SetTitleSize(0.05);
        mg->GetXaxis()->SetLabelSize(0.04);
        mg->GetYaxis()->SetLabelSize(0.04);
        legend->Draw();


        c2->cd(i+1);
        gPad->SetLeftMargin(0.1);
        gPad->SetBottomMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.10);
        TGraphErrors *agr = new TGraphErrors(overT.size(), overT.data(), logIoverT2.data(), err_overT.data(), err_logIoverT2.data());
        agr->SetName(Form("Arrhenius_PiN_%s", tree_names[i].c_str()));
        agr->SetTitle(Form("Arrhenius Plot for #Phi = %s n_{eq}/cm^{2}", fluences[i].c_str()));
        agr->SetMarkerStyle(20);
        agr->SetMarkerSize(0.5);
        agr->SetLineColor(colors[i]);
        agr->SetMarkerColor(colors[i]);
        agr->Draw("AP");
        agr->GetXaxis()->SetTitle("1000/T [1/K]");
        agr->GetYaxis()->SetTitle("log(I/T^{2}) [u.a.]");
        agr->GetXaxis()->SetTitleSize(0.05);
        agr->GetYaxis()->SetTitleSize(0.05);
        agr->GetXaxis()->SetLabelSize(0.04);
        agr->GetYaxis()->SetLabelSize(0.04);
        mg_arr->Add(agr);
        legend_arr->AddEntry(agr, Form("T = %.1f C", T), "pl");

        TF1 *arr_fit = new TF1(Form("arr_fit_%s", tree_names[i].c_str()), "[0]*x + [1]", 3.3, 4.1);
        arr_fit->SetLineColor(colors[i]);
        agr->Fit(arr_fit, "RQ");


        cout << "Fluence: " << fluences[i] << " n_{eq}/cm^{2}" << endl;
        cout << "Fit parameters: " << endl;
        cout << "Slope: " << arr_fit->GetParameter(0) << " +/- " << arr_fit->GetParError(0) << endl;
        cout << "Intercept: " << arr_fit->GetParameter(1) << " +/- " << arr_fit->GetParError(1) << endl;
        cout << "Chi2/NDF: " << arr_fit->GetChisquare() << " / " << arr_fit->GetNDF() << endl;
        cout << "P-value: " << arr_fit->GetProb() << endl;

        cout << "\nExpected energy gap: " << -arr_fit->GetParameter(0) * 2000. * Kb << " +/- " << arr_fit->GetParError(0) * 2000. * Kb << " eV" << endl;
        cout << "----------------------------------------" << endl;


    }

    c3->cd();
    gPad->SetLeftMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.10);
    mg_arr->SetTitle("Arrhenius Plot for all fluences");
    mg_arr->Draw("AP");
    mg_arr->GetXaxis()->SetTitle("1000/T [1/K]");
    mg_arr->GetYaxis()->SetTitle("log(I/T^{2}) [u.a.]");
    mg_arr->GetXaxis()->SetTitleSize(0.05);
    mg_arr->GetYaxis()->SetTitleSize(0.05);
    mg_arr->GetXaxis()->SetLabelSize(0.04);
    mg_arr->GetYaxis()->SetLabelSize(0.04);
    legend_arr->Draw();
}

int find_100V(const vector<double>& V)
{
    double mindist = 100.;
    int index = -1;
    for(int i=0;i<V.size();i++)
    {
        double dist = fabs(V[i] - 100.);
        if(dist < mindist)
        {
            mindist = dist;
            index = i;
        }
    }
    return index;
}