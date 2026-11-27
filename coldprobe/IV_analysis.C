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


/*
    Incertezze mancanti:
    - Valutare se 0.3C ci sta come incertezza sulla temperatura
    - Manca incertezza sulla tensione
    - Manca incertezza sulla corrente (nominale)
    - Manca incertezza sulla corrente dovuta all'incertezza sulla temperatura (ogni grado fa circa 10% di incertezza)

    Potremmo mettere 0.5 gradi per valutare il fatto che misuriamo la temperatura del backplane e non del pad, e il corrispondente 5% di incertezza relativa sulla corrente
*/

#define Kb 8.617333262145e-5 // eV/K

int find_100V(const vector<double>&);

void IV_analysis()
{
    double err_temperature = 0.5;

    double pixel_volume = 0.13 * 0.13 * 0.003; // cm^3

    int ordJ[] = {8,6,4,2,0,1,3,5,7}; //Mi dispiace

    string filename = "IV_data.root";
    string tree_names[4] = {"New", "1e14", "8e14", "15e14"};
    string fluences[4] = {"0.0", "1e14", "8e14", "15e14"};
    TFile *file = new TFile(filename.c_str(), "READ");

    int colors[] = {kBlack, kRed, kBlue, kGreen+2, kMagenta, kCyan+2, kOrange+7, kViolet+1, kTeal+1, kAzure+1};
    
    vector<double> *V=nullptr, *I=nullptr, *V_err=nullptr, *I_err=nullptr;

    vector<vector<double>> vI100(4), err_vI100(4), vT100(4), err_vT100(4);
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
            err_overT.push_back(1000.*err_temperature/(T_K*T_K));
            int index_100V = find_100V(*V);
            double I100 = (*I)[index_100V];
            double err_I100 = (*I_err)[index_100V];
            logIoverT2.push_back(log(I100) - 2.*log(T_K));

            vI100[i].push_back(I100);
            err_vI100[i].push_back(err_I100);
            vT100[i].push_back(T);
            err_vT100[i].push_back(err_temperature);

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
        gPad->SetLeftMargin(0.15);
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
        legend_arr->AddEntry(agr, Form("#Phi = %s n_{eq}/cm^{2}", fluences[i].c_str()), "pl");

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


    TCanvas *c4 = new TCanvas("c4", "Damage factor", 1100, 600);
    TMultiGraph *mg_df = new TMultiGraph();
    TLegend *legend_df = new TLegend(0.75, 0.7, 0.95, 0.9);

    double df_coeff[9], df_coeff_err[9];

    for(int i=0;i<9;i++)
    {
        double I2[] = {vI100[0][i], vI100[1][i], vI100[2][i], vI100[3][i]};
        double err_I2[] = {err_vI100[0][i], err_vI100[1][i], err_vI100[2][i], err_vI100[3][i]};
        double fluence2[] = {0., 1e14, 8e14, 15e14};
        double err_fluence2[] = {0., 1e13, 8e13, 1.5e14};
        TGraphErrors *gr2 = new TGraphErrors(4, fluence2, I2, err_fluence2, err_I2);
        gr2->SetName(Form("Damage_factor_T%.1fC", vT100[0][i]));
        gr2->SetTitle(Form("Damage factor for T = %.1f C", vT100[0][i]));
        gr2->SetMarkerStyle(20);
        gr2->SetMarkerSize(0.5);
        gr2->SetLineColor(colors[i]);
        gr2->SetMarkerColor(colors[i]);
        mg_df->Add(gr2);
        legend_df->AddEntry(gr2, Form("T = %.1f C", vT100[0][i]), "pl");

        TF1 *df_fit = new TF1(Form("df_fit_T%.1fC", vT100[0][i]), "[0]*x + [1]", 0.5e14, 2e15);
        df_fit->SetLineColor(colors[i]);
        df_fit->SetParameter(0, I2[3]/fluence2[3]);
        df_fit->SetParameter(1, 0);
        gr2->Fit(df_fit, "RQ");
        df_fit->Draw("same");
        cout << "T = " << vT100[0][i] << " C" << endl;
        cout << "Fit parameters: " << endl;
        cout << "Slope: " << df_fit->GetParameter(0) << " +/- " << df_fit->GetParError(0) << endl;
        cout << "Intercept: " << df_fit->GetParameter(1) << " +/- " << df_fit->GetParError(1) << endl;
        cout << "Chi2/NDF: " << df_fit->GetChisquare() << " / " << df_fit->GetNDF() << endl;
        cout << "P-value: " << df_fit->GetProb() << endl;
        cout << "----------------------------------------" << endl;

        df_coeff[i] = df_fit->GetParameter(0)/pixel_volume*1e9;
        df_coeff_err[i] = df_fit->GetParError(0)/pixel_volume*1e9;
    }

    c4->cd();
    gPad->SetLeftMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.10);
    mg_df->SetTitle("Damage factor for all temperatures");
    mg_df->Draw("AP");
    mg_df->GetXaxis()->SetTitle("#Phi [n_{eq}/cm^{2}]");
    mg_df->GetYaxis()->SetTitle("I_{100V} [nA]");
    mg_df->GetXaxis()->SetTitleSize(0.05);
    mg_df->GetYaxis()->SetTitleSize(0.05);
    mg_df->GetXaxis()->SetLabelSize(0.04);
    mg_df->GetYaxis()->SetLabelSize(0.04);
    legend_df->Draw();


    TCanvas *c5 = new TCanvas("c5", "Damage factor vs Temperature", 900, 600);
    gPad->SetLeftMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.10);
    TGraphErrors *gr_dfT = new TGraphErrors(9, vT100[0].data(), df_coeff, err_vT100[0].data(), df_coeff_err);
    gr_dfT->SetName("Damage_factor_vs_Temperature");
    gr_dfT->SetTitle("Damage factor vs Temperature");
    gr_dfT->SetMarkerStyle(20);
    gr_dfT->SetMarkerSize(0.5);
    gr_dfT->Draw("AP");
    gr_dfT->GetXaxis()->SetTitle("T [C]");
    gr_dfT->GetYaxis()->SetTitle("Damage factor [10^{-18} A/(n_{eq}*cm)]");
    gr_dfT->GetXaxis()->SetTitleSize(0.05);
    gr_dfT->GetYaxis()->SetTitleSize(0.05);
    gr_dfT->GetXaxis()->SetLabelSize(0.04);
    gr_dfT->GetYaxis()->SetLabelSize(0.04);

    TF1 *dfT_fit = new TF1("dfT_fit", "[0]*exp([1]*x)", -30, 30);
    dfT_fit->SetLineColor(kRed);
    dfT_fit->SetParameter(0, 5);
    dfT_fit->SetParameter(1, 0.10);
    gr_dfT->Fit(dfT_fit, "RQ");
    dfT_fit->Draw("same");

    cout << "----------------------------------------" << endl;
    cout << "Fit p value: " << dfT_fit->GetProb() << endl;
    double df_20C = dfT_fit->Eval(20);
    double df_20C_err = df_20C * sqrt(pow(dfT_fit->GetParError(0)/dfT_fit->GetParameter(0), 2) + pow(dfT_fit->GetParError(1)*20, 2));

    cout << "Damage factor at 20 C: (" << dfT_fit->Eval(20) << " +/- " << df_20C_err << ") * 10^{-18} A/(n_{eq}*cm)" << endl;

    TGraphErrors *point20C = new TGraphErrors(1, new double[1]{20}, new double[1]{df_20C}, new double[1]{0}, new double[1]{df_20C_err});
    point20C->SetMarkerStyle(21);
    point20C->SetMarkerSize(1);
    point20C->SetMarkerColor(kRed);
    point20C->Draw("P same");
   
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