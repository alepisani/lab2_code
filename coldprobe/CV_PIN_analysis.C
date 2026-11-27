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

#define epsilonSi 1.05e-12  // F/cm
#define q 1.602e-19         // C
#define area 0.13*0.13      // cm^2

void CV_PIN_analysis()
{
    int colors[] = {kBlack, kRed, kBlue, kGreen+2, kMagenta, kCyan+2, kOrange+7, kViolet+1, kTeal+1, kAzure+1};
    string filename = "CV_data.root";
    string tree_names[4] = {"new", "1e14", "8e14", "15e14"};
    string fluences[4] = {"0.0", "1e14", "8e14", "15e14"};
    TFile *file = new TFile(filename.c_str(), "READ");

    vector<vector<double>> V(4), C(4), overC2(4), V_err(4), C_err(4), overC2_err(4);

    for(int i=0;i<4;i++)
    {
        double v,c, v_err, c_err;
        TTree *tree = (TTree*)file->Get(Form("CV_PiN_25C_%s", tree_names[i].c_str()));
        tree->SetBranchAddress("V", &v);
        tree->SetBranchAddress("C", &c);
        tree->SetBranchAddress("V_err", &v_err);
        tree->SetBranchAddress("C_err", &c_err);

        for(int j=0;j<tree->GetEntries();j++)
        {
            tree->GetEntry(j);
            V[i].push_back(v);
            C[i].push_back(c);
            overC2[i].push_back(1/(c*c));
            V_err[i].push_back(v_err);
            C_err[i].push_back(c_err);
            overC2_err[i].push_back(2*c_err/(c*c*c));
        }
    }
    file->Close();

    TCanvas *c1 = new TCanvas("c1", "C-V curves", 900, 600);
    TMultiGraph *mg = new TMultiGraph();
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->SetTextSize(0.03);

    for(int i=0;i<4;i++)
    {
        TGraphErrors *gr = new TGraphErrors(V[i].size(), &V[i][0], &C[i][0], &V_err[i][0], &C_err[i][0]);
        gr->SetLineColor(colors[i]);
        gr->SetMarkerColor(colors[i]);
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(0.3);
        mg->Add(gr, "P");
        legend->AddEntry(gr, Form("#Phi = %s n_{eq}/cm^{2}", fluences[i].c_str()), "lp");
    }
    mg->Draw("AP");
    mg->GetXaxis()->SetTitle("Voltage [V]");
    mg->GetYaxis()->SetTitle("Capacitance [pF]");
    mg->GetXaxis()->SetTitleSize(0.03);
    mg->GetYaxis()->SetTitleSize(0.03);
    mg->GetXaxis()->SetLabelSize(0.03);
    mg->GetYaxis()->SetLabelSize(0.03);
    mg->GetXaxis()->SetTitleOffset(1.2);
    mg->GetYaxis()->SetTitleOffset(1.2);
    legend->Draw();

    TCanvas *c2 = new TCanvas("c2", "1/C^2-V curves", 900, 600);
    TMultiGraph *mg2 = new TMultiGraph();
    TLegend *legend2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend2->SetTextSize(0.03);

    for(int i=0;i<4;i++)
    {
        TGraphErrors *gr = new TGraphErrors(V[i].size(), &V[i][0], &overC2[i][0], &V_err[i][0], &overC2_err[i][0]);
        gr->SetLineColor(colors[i]);
        gr->SetMarkerColor(colors[i]);
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(0.3);
        mg2->Add(gr, "P");
        legend2->AddEntry(gr, Form("#Phi = %s n_{eq}/cm^{2}", fluences[i].c_str()), "lp");
    }
    mg2->Draw("AP");
    mg2->GetXaxis()->SetTitle("Voltage [V]");
    mg2->GetYaxis()->SetTitle("1/C^2 [1/pF^2]");
    mg2->GetXaxis()->SetTitleSize(0.03);
    mg2->GetYaxis()->SetTitleSize(0.03);
    mg2->GetXaxis()->SetLabelSize(0.03);
    mg2->GetYaxis()->SetLabelSize(0.03);
    mg2->GetXaxis()->SetTitleOffset(1.2);
    mg2->GetYaxis()->SetTitleOffset(1.2);
    legend2->Draw(); 

    TCanvas *c3 = new TCanvas("c3", "C-V curves", 1800, 1200);
    c3->Divide(2,2);

    double xmin1[4]={1,1,4,4};
    double xmax1[4]={5,3.5,7,7.5};

    double xmin2[4]={20,20,20,30};
    double xmax2[4]={30,30,30,40};


    double Vfd[4], Vfd_err[4];

    for(int i=0;i<4;i++)
    {
        c3->cd(i+1);
        TGraphErrors *gr = new TGraphErrors(V[i].size(), &V[i][0], &overC2[i][0], &V_err[i][0], &overC2_err[i][0]);
        gr->SetName(Form("overC2_PiN_%s", tree_names[i].c_str()));
        gr->SetTitle(Form("#Phi = %s n_{eq}/cm^{2}", fluences[i].c_str()));
        gr->SetLineColor(colors[i]);
        gr->SetMarkerColor(colors[i]);
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(0.3);
        gr->Draw("AP");
        gr->GetXaxis()->SetTitle("Voltage [V]");
        gr->GetYaxis()->SetTitle("1/C^{2} [1/pF^{2}]");
        gr->GetXaxis()->SetTitleSize(0.03);
        gr->GetYaxis()->SetTitleSize(0.03);
        gr->GetXaxis()->SetLabelSize(0.03);
        gr->GetYaxis()->SetLabelSize(0.03);
        gr->GetXaxis()->SetTitleOffset(1.2);
        gr->GetYaxis()->SetTitleOffset(1.2);


        TF1 *fit_rise = new TF1(Form("fit_%s", tree_names[i].c_str()), "[0]*x + [1]", xmin1[i], xmax1[i]);
        fit_rise->SetLineColor(colors[i]);
        gr->Fit(fit_rise, "RQ");
        fit_rise->Draw("same");

        cout << " Rise fit for fluence " << fluences[i] << " n_{eq}/cm^{2}: slope = " << fit_rise->GetParameter(0) << " +/- " << fit_rise->GetParError(0) << endl;
        cout << "\t p value = " << fit_rise->GetProb() << endl;

        TF1 *fit_plateau = new TF1(Form("fit_plateau_%s", tree_names[i].c_str()), "[0]*x + [1]", xmin2[i], xmax2[i]);
        fit_plateau->SetLineColor(colors[i]);
        gr->Fit(fit_plateau, "RQ+");
        fit_plateau->Draw("same");

        cout << " Plateau fit for fluence " << fluences[i] << " n_{eq}/cm^{2}: slope = " << fit_plateau->GetParameter(0) << " +/- " << fit_plateau->GetParError(0) << endl;
        cout << "\t p value = " << fit_plateau->GetProb() << endl;

        Vfd[i] = (fit_plateau->GetParameter(1) - fit_rise->GetParameter(1)) / (fit_rise->GetParameter(0) - fit_plateau->GetParameter(0));
        Vfd_err[i] = sqrt(pow(fit_plateau->GetParError(1)/(fit_rise->GetParameter(0) - fit_plateau->GetParameter(0)), 2) +
                          pow(fit_rise->GetParError(1)/(fit_rise->GetParameter(0) - fit_plateau->GetParameter(0)), 2) +
                          pow((fit_plateau->GetParameter(1) - fit_rise->GetParameter(1)) * fit_rise->GetParError(0)/pow(fit_rise->GetParameter(0) - fit_plateau->GetParameter(0), 2), 2) +
                          pow((fit_plateau->GetParameter(1) - fit_rise->GetParameter(1)) * fit_plateau->GetParError(0)/pow(fit_rise->GetParameter(0) - fit_plateau->GetParameter(0), 2), 2));
    }

    TCanvas *c4 = new TCanvas("c4", "Neff vs fluence", 900, 600);

    double N_eff[4], N_eff_err[4];
    for(int i=0;i<4;i++)
    {
        N_eff[i] = 2*epsilonSi * Vfd[i] / (q * 0.003 * 0.003);
        N_eff_err[i] = 2*epsilonSi * Vfd_err[i] / (q * 0.003 * 0.003);
    }

    TGraphErrors *gr_Neff = new TGraphErrors(4, new double[4]{0, 1e14, 8e14, 15e14}, N_eff, new double[4]{0, 0, 0, 0}, N_eff_err);
    gr_Neff->SetTitle("N_{eff} vs fluence");
    gr_Neff->SetMarkerStyle(20);
    gr_Neff->SetMarkerSize(1);
    gr_Neff->GetXaxis()->SetTitle("#Phi [n_{eq}/cm^{2}]");
    gr_Neff->GetYaxis()->SetTitle("N_{eff} [cm^{-3}]");
    gr_Neff->GetXaxis()->SetTitleSize(0.03);
    gr_Neff->GetYaxis()->SetTitleSize(0.03);
    gr_Neff->GetXaxis()->SetLabelSize(0.03);
    gr_Neff->GetYaxis()->SetLabelSize(0.03);
    gr_Neff->GetXaxis()->SetTitleOffset(1.2);
    gr_Neff->GetYaxis()->SetTitleOffset(1.2);
    gr_Neff->Draw("AP");

    TF1 *fit_Neff = new TF1("fit_Neff", "[0]*exp(-[1]*x) + [2]*x", 0, 16e14);
    fit_Neff->SetLineColor(kRed);
    fit_Neff->FixParameter(0, N_eff[0]);
    fit_Neff->SetParameter(1, 1e-14);
    fit_Neff->SetParameter(2, 1e-2);
    gr_Neff->Fit(fit_Neff, "RQ");
    fit_Neff->Draw("same");


    vector<vector<double>> dCdV(4);
    vector<vector<double>> V_dCdV(4);
    vector<vector<double>> Neff_dyn(4);
    vector<vector<double>> depl_w(4);

    TCanvas *c5 = new TCanvas("c5", "dC/dV curves", 1800, 1200);
    c5->Divide(2,2);

    for(int i=0;i<4;i++)
    {
        dCdV[i].push_back((overC2[i][1]-overC2[i][0])/(V[i][1]-V[i][0]));
        for(int j=1;j<V[i].size()-2;j++)
            dCdV[i].push_back((overC2[i][j+1]-overC2[i][j-1])/(V[i][j+1]-V[i][j-1]));
        dCdV[i].push_back((overC2[i][V[i].size()-1]-overC2[i][V[i].size()-2])/(V[i][V[i].size()-1]-V[i][V[i].size()-2]));

        for(int j=0;j<dCdV[i].size();j++)
        {
            depl_w[i].push_back(epsilonSi * area / C[i][j]);
            Neff_dyn[i].push_back(2./(q * epsilonSi * area * area) * 1./dCdV[i][j]);
        }

        TGraph *gr_profile = new TGraph(depl_w[i].size(), &depl_w[i][0], &Neff_dyn[i][0]);
        gr_profile->SetName(Form("Neff_profile_%s", tree_names[i].c_str()));
        gr_profile->SetTitle(Form("Neff profile for #Phi = %s n_{eq}/cm^{2}", fluences[i].c_str()));
        gr_profile->SetLineColor(colors[i]);
        gr_profile->SetMarkerColor(colors[i]);
        gr_profile->SetMarkerStyle(20);
        c5->cd(i+1);
        gr_profile->Draw("AP");
        gr_profile->GetXaxis()->SetTitle("Depletion width [cm]");
        gr_profile->GetYaxis()->SetTitle("N_{eff} [cm^{-3}]");
        gr_profile->GetXaxis()->SetTitleSize(0.03);
        gr_profile->GetYaxis()->SetTitleSize(0.03);
        gr_profile->GetXaxis()->SetLabelSize(0.03);
        gr_profile->GetYaxis()->SetLabelSize(0.03);
        gr_profile->GetXaxis()->SetTitleOffset(1.2);
        gr_profile->GetYaxis()->SetTitleOffset(1.2);
    }
}