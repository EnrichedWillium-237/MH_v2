# include "TCanvas.h"
# include "TDirectory.h"
# include "TFile.h"
# include "TGraphErrors.h"
# include "TH1.h"
# include "TH2.h"
# include "TLegend.h"
# include "TMath.h"
# include "TPaveText.h"
# include "Riostream.h"
# include <iostream>

using namespace std;

# include "src/HiEvtPlaneList.h"
# include "src/style.h"

using namespace hi;
static const int ncentbins = 13;
static const int cminCENT[] = {0,  5, 10, 15, 20, 25, 30, 35, 40, 50, 60,  0, 20,  60};
static const int cmaxCENT[] = {5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 20, 60, 100};
static const int netabins = 12;
static const double etabins[] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0,  0.4,  0.8,  1.2,  1.6,  2.0,  2.4};
static const int nbinsETA = 14;
static const double eminETA[] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0, -2.4,  0.0};
static const double emaxETA[] = {-2.0, -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0,  2.4,  0.0,  2.4};
static const double etaMid[] = {-2.2, -1.8, -1.4, -1.0, -0.6, -0.2,  0.2,  0.6,  1.0,  1.4,  1.8,  2.2};
TString etags[] = {"-24_-20", "-20_-16", "-16_-12", "-12_-8", "-8_-4", "-4_0", "0_4", "4_8", "8_12", "12_16", "16_20", "20_24",
                   "-24_0", "0_24"};
static const int nptbins = 18;
static const double ptbins[] = {0.30,  0.40,  0.50,  0.60,  0.80,  1.00,  1.25,  1.50,  2.00,  2.50,  3.00,
                     3.50,  4.00,  5.00,  6.00,  7.00,  8.00,  10.00,  12.00};
static int NANALS = 16;
TString ANAL[] = {"N1MCm22SUB2", "N1MCp22SUB2", "N1MCm22SUB3", "N1MCp22SUB3",
                  "N1SUB2",      "N1ASUB2",     "N1BSUB2",     "N112SUB2",    "N112ASUB2",   "N112BSUB2",
                  "N1SUB3",      "N1ASUB3",     "N1BSUB3"};

# include "src/getRatios.h"

TFile * finPt;
TFile * finEta;
TFile * fout;

void comparison()
{
    TH1::SetDefaultSumw2();

    finPt = new TFile("hists/MH_combined_Pt.root","read");

    for (int ebin = 0; ebin<nbinsETA; ebin++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            TString mtag = Form("MH_nominal/eta_%s/%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]);

            N1MCm22SUB2_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1MCm22SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MCp22SUB2_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1MCp22SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MCm22SUB3_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1MCm22SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MCp22SUB3_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1MCp22SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1SUB2_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1ASUB2_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1ASUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1BSUB2_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1BSUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112SUB2_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N112SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112ASUB2_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N112ASUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112BSUB2_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N112BSUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1SUB3_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1ASUB3_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1ASUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1BSUB3_pT[ebin][cbin] = (TH1D *) finPt->Get(Form("%s/N1BSUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
        }
    }

    finEta = new TFile("hists/MH_combined_Eta.root","read");

    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TString mtag = Form("MH_nominal/%d_%d",cminCENT[cbin],cmaxCENT[cbin]);

        N1MCm22SUB2_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MCm22SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp22SUB2_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MCp22SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MC22SUB2_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MC22SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCm22SUB3_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MCm22SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MCp22SUB3_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MCp22SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1MC22SUB3_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1MC22SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1SUB2_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1ASUB2_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1ASUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1BSUB2_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1BSUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112SUB2_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N112SUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112ASUB2_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N112ASUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N112BSUB2_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N112BSUB2_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1SUB3_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1SUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1ASUB3_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1ASUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
        N1BSUB3_eta[cbin] = (TH1D *) finEta->Get(Form("%s/N1BSUB3_%d_%d",mtag.Data(),cminCENT[cbin],cmaxCENT[cbin]));
    }

    # include "../../published_results/PHOBOS_AuAu.h"

    TH1D * PHOBOS_AuAu_Mixed_200GeV = new TH1D("PHOBOS_AuAu_Mixed_200GeV", "", PHOBOS_AuAu_nebins, PHOBOS_AuAu_eta);
    TH1D * PHOBOS_AuAu_2sub_200GeV = new TH1D("PHOBOS_AuAu_2sub_200GeV", "", PHOBOS_AuAu_nebins, PHOBOS_AuAu_eta);
    for (int i = 0; i<PHOBOS_AuAu_nebins; i++) {
        PHOBOS_AuAu_Mixed_200GeV->SetBinContent(i+1, PHOBOS_AuAu_mixed_v1_3[i]);
        PHOBOS_AuAu_Mixed_200GeV->SetBinError(i+1, PHOBOS_AuAu_mixed_v1Err_3[i]);

        PHOBOS_AuAu_Mixed_200GeV->SetMarkerColor(kBlack);
        PHOBOS_AuAu_Mixed_200GeV->SetLineColor(kBlack);
        PHOBOS_AuAu_Mixed_200GeV->SetMarkerStyle(24);
        PHOBOS_AuAu_Mixed_200GeV->SetMarkerSize(1.2);

        PHOBOS_AuAu_2sub_200GeV->SetBinContent(i+1, PHOBOS_AuAu_sym_v1_3[i]);
        PHOBOS_AuAu_2sub_200GeV->SetBinError(i+1, PHOBOS_AuAu_sym_v1Err_3[i]);

        PHOBOS_AuAu_2sub_200GeV->SetMarkerColor(kBlack);
        PHOBOS_AuAu_2sub_200GeV->SetLineColor(kBlack);
        PHOBOS_AuAu_2sub_200GeV->SetMarkerStyle(21);
        PHOBOS_AuAu_2sub_200GeV->SetMarkerSize(1.2);
    }


    TCanvas * c0 = new TCanvas("c0","c0",650,600);
    TPad * pad0 = (TPad *) c0->cd();
    pad0->SetGrid();
    TH1D * h0 = new TH1D("h0", "", 100, -2.6, 2.6);
    h0->SetXTitle("#eta");
    h0->SetYTitle("v_{1}^{odd}");
    h0->GetYaxis()->SetRangeUser(-0.04, 0.04);
    h0->Draw();
    PHOBOS_AuAu_2sub_200GeV->Draw("same");
    PHOBOS_AuAu_Mixed_200GeV->Draw("same");
    TLegend * leg0 = new TLegend(0.2, 0.2, 0.4, 0.4, "NDC");
    SetLegend(leg0, 18);
    leg0->AddEntry(PHOBOS_AuAu_2sub_200GeV,"PHOBOS 2SUB #sqrt{s_{NN}}=200 GeV (6-40%)","p");
    leg0->AddEntry(PHOBOS_AuAu_Mixed_200GeV,"PHOBOS Mixed #sqrt{s_{NN}}=200 GeV (6-40%)","p");
    leg0->Draw();



}
