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
static const int netabins = 14;
static const double eminETA[] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0, -2.4,  0.0};
static const double emaxETA[] = {-2.0, -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0,  2.4,  0.0,  2.4};
static const double etaMid[] = {-2.2, -1.8, -1.4, -1.0, -0.6, -0.2,  0.2,  0.6,  1.0,  1.4,  1.8,  2.2};
TString etags[] = {"-24_-20", "-20_-16", "-16_-12", "-12_-8", "-8_-4", "-4_0", "0_4", "4_8", "8_12", "12_16", "16_20", "20_24", "-24_0", "0_24"};
static const int nptbins = 18;
static const double ptbins[] = {0.30,  0.40,  0.50,  0.60,  0.80,  1.00,  1.25,  1.50,  2.00,  2.50,  3.00,
                     3.50,  4.00,  5.00,  6.00,  7.00,  8.00,  10.00,  12.00};
static int NANALS = 24;
TString ANAL[] = {"N1MCm22SUB2", "N1MCm18SUB2", "N1MCm14SUB2", "N1MCp22SUB2", "N1MCp18SUB2", "N1MCp14SUB2",
                  "N1MCm22SUB3", "N1MCm18SUB3", "N1MCm14SUB3", "N1MCp22SUB3", "N1MCp18SUB3", "N1MCp14SUB3",
                  "N1SUB2",      "N1ASUB2",     "N1BSUB2",     "N112SUB2",    "N112ASUB2",   "N112BSUB2",
                  "N1SUB3",      "N1ASUB3",     "N1BSUB3",     "N112SUB3",    "N112ASUB3",   "N112BSUB3"
};

TH1D * N1MCm22SUB2[netabins][ncentbins];
TH1D * N1MCp22SUB2[netabins][ncentbins];
TH1D * N1MC22SUB2[netabins][ncentbins];
TH1D * N1MCm18SUB2[netabins][ncentbins];
TH1D * N1MCp18SUB2[netabins][ncentbins];
TH1D * N1MC18SUB2[netabins][ncentbins];

TH1D * N1MCm22SUB3[netabins][ncentbins];
TH1D * N1MCp22SUB3[netabins][ncentbins];
TH1D * N1MC22SUB3[netabins][ncentbins];
TH1D * N1MCm18SUB3[netabins][ncentbins];
TH1D * N1MCp18SUB3[netabins][ncentbins];
TH1D * N1MC18SUB3[netabins][ncentbins];

TH1D * N1ASUB2[netabins][ncentbins];
TH1D * N1BSUB2[netabins][ncentbins];
TH1D * N1SUB2[netabins][ncentbins];
TH1D * N112ASUB2[netabins][ncentbins];
TH1D * N112BSUB2[netabins][ncentbins];
TH1D * N112SUB2[netabins][ncentbins];

TH1D * N1ASUB3[netabins][ncentbins];
TH1D * N1BSUB3[netabins][ncentbins];
TH1D * N1SUB3[netabins][ncentbins];
TH1D * N112ASUB3[netabins][ncentbins];
TH1D * N112BSUB3[netabins][ncentbins];
TH1D * N112SUB3[netabins][ncentbins];

TFile * tfin;

void plotPtDist()
{
    tfin = new TFile("results/MH_combined_Pt.root","read");

    for (int ebin = 0; ebin<netabins; ebin++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            TString mtag = Form("MH_nominal/eta_%s/%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]);

            N1MCm22SUB2[ebin][cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm22SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MCp22SUB2[ebin][cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp22SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MC22SUB2[ebin][cbin] = (TH1D *) N1MCp22SUB2[ebin][cbin]->Clone(Form("N1MC22SUB2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MC22SUB2[ebin][cbin]->Add(N1MCm22SUB2[ebin][cbin]);
            N1MC22SUB2[ebin][cbin]->Scale(0.5);

            N1MCm18SUB2[ebin][cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm18SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MCp18SUB2[ebin][cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp18SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MC18SUB2[ebin][cbin] = (TH1D *) N1MCp18SUB2[ebin][cbin]->Clone(Form("N1MC18SUB2_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MC18SUB2[ebin][cbin]->Add(N1MCm18SUB2[ebin][cbin]);
            N1MC18SUB2[ebin][cbin]->Scale(0.5);

            N1MCm22SUB3[ebin][cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm22SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MCp22SUB3[ebin][cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp22SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MC22SUB3[ebin][cbin] = (TH1D *) N1MCp22SUB3[ebin][cbin]->Clone(Form("N1MC22SUB3_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MC22SUB3[ebin][cbin]->Add(N1MCm22SUB3[ebin][cbin]);
            N1MC22SUB3[ebin][cbin]->Scale(0.5);

            N1MCm18SUB3[ebin][cbin] = (TH1D *) tfin->Get(Form("%s/N1MCm18SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MCp18SUB3[ebin][cbin] = (TH1D *) tfin->Get(Form("%s/N1MCp18SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MC18SUB3[ebin][cbin] = (TH1D *) N1MCp18SUB3[ebin][cbin]->Clone(Form("N1MC18SUB3_eta_%s_%d_%d",etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1MC18SUB3[ebin][cbin]->Add(N1MCm18SUB3[ebin][cbin]);
            N1MC18SUB3[ebin][cbin]->Scale(0.5);

            N1ASUB2[ebin][cbin] = (TH1D *) tfin->Get(Form("%s/N1ASUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1BSUB2[ebin][cbin] = (TH1D *) tfin->Get(Form("%s/N1BSUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1SUB2[ebin][cbin] = (TH1D *) tfin->Get(Form("%s/N1SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112ASUB2[ebin][cbin] = (TH1D *) tfin->Get(Form("%s/N112ASUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112BSUB2[ebin][cbin] = (TH1D *) tfin->Get(Form("%s/N112BSUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112SUB2[ebin][cbin] = (TH1D *) tfin->Get(Form("%s/N112SUB2_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));

            N1ASUB3[ebin][cbin] = (TH1D *) tfin->Get(Form("%s/N1ASUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1BSUB3[ebin][cbin] = (TH1D *) tfin->Get(Form("%s/N1BSUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N1SUB3[ebin][cbin] = (TH1D *) tfin->Get(Form("%s/N1SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112ASUB3[ebin][cbin] = (TH1D *) tfin->Get(Form("%s/N112ASUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112BSUB3[ebin][cbin] = (TH1D *) tfin->Get(Form("%s/N112BSUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
            N112SUB3[ebin][cbin] = (TH1D *) tfin->Get(Form("%s/N112SUB3_eta_%s_%d_%d",mtag.Data(),etags[ebin].Data(),cminCENT[cbin],cmaxCENT[cbin]));
        }
    }

    if (!fopen("results/Pt_Distributions","r")) system("mkdir results/Pt_Distributions");

    int erangeNeg = 0;
    int erangePos = 0;
    int crange = 0;
    TH1D * h0 = new TH1D("h0", "", 100, 0.0, 12.0);
    h0->SetStats(0);
    h0->SetXTitle("p_{T} (GeV/c)");
    h0->SetYTitle("v_{1}");
    h0->GetYaxis()->SetRangeUser(-0.03, 0.2);

    erangeNeg = 13;
    erangePos = 12;
    crange = 1;
    TCanvas * cN1MCpm22SUB2 = new TCanvas("cN1MCpm22SUB2","cN1MCpm22SUB2",650,600);
    TPad * padN1MCpm22SUB2 = (TPad *) cN1MCpm22SUB2->cd();
    padN1MCpm22SUB2->SetGrid();
    TH1D * hN1MCpm22SUB2 = (TH1D *) h0->Clone("hN1MCpm22SUB2");
    hN1MCpm22SUB2->SetYTitle("v_{1}^{even}");
    hN1MCpm22SUB2->GetYaxis()->SetRangeUser(-0.03, 0.16);
    hN1MCpm22SUB2->Draw();
    N1MCm22SUB2[erangeNeg][crange]->SetMarkerColor(kRed);
    N1MCm22SUB2[erangeNeg][crange]->SetLineColor(kRed);
    N1MCm22SUB2[erangeNeg][crange]->SetMarkerStyle(24);
    N1MCm22SUB2[erangeNeg][crange]->SetMarkerSize(1.2);
    N1MCm22SUB2[erangeNeg][crange]->Draw("same");
    N1MCp22SUB2[erangePos][crange]->SetMarkerColor(kBlue);
    N1MCp22SUB2[erangePos][crange]->SetLineColor(kBlue);
    N1MCp22SUB2[erangePos][crange]->SetMarkerStyle(24);
    N1MCp22SUB2[erangePos][crange]->SetMarkerSize(1.2);
    N1MCp22SUB2[erangePos][crange]->Draw("same");
    TLegend * legN1MCpm22SUB2 = new TLegend(0.49, 0.17, 0.78, 0.27);
    SetLegend(legN1MCpm22SUB2, 18);
    legN1MCpm22SUB2->AddEntry(N1MCm22SUB2[erangeNeg][crange],"#Psi_{1}^{trk}{2SUB} (-2.4 < #eta < -2.0)","p");
    legN1MCpm22SUB2->AddEntry(N1MCp22SUB2[erangePos][crange],"#Psi_{1}^{trk}{2SUB} (2.0 < #eta < 2.4)","p");
    legN1MCpm22SUB2->Draw();
    TPaveText * txN1MCpm22SUB2 = new TPaveText(0.19, 0.84, 0.44, 0.93, "NDC");
    SetTPaveTxt(txN1MCpm22SUB2, 18);
    txN1MCpm22SUB2->AddText(Form("POI: %1.1f < |#eta| < %1.1f",eminETA[erangeNeg],emaxETA[erangeNeg]));
    txN1MCpm22SUB2->AddText(Form("(%d-%d%%)",cminCENT[crange],cmaxCENT[crange]));
    txN1MCpm22SUB2->Draw();
    cN1MCpm22SUB2->Print(Form("results/Pt_Distributions/N1MCpm22SUB2_%d_%d.png",cminCENT[crange],cmaxCENT[crange]),"png");


    TCanvas * cratN1MCpm22SUB2 = new TCanvas("cratN1MCpm22SUB2","cratN1MCpm22SUB2",650,600);
    TPad * padratN1MCpm22SUB2 = (TPad *) cratN1MCpm22SUB2->cd();
    padratN1MCpm22SUB2->SetGrid();
    TH1D * hratN1MCpm22SUB2 = (TH1D *) h0->Clone("hratN1MCpm22SUB2");
    hratN1MCpm22SUB2->SetYTitle("v_{1}^{even}{+#eta} / v_{1}^{even}{-#eta}");
    hratN1MCpm22SUB2->GetYaxis()->SetRangeUser(0.6, 1.4);
    hratN1MCpm22SUB2->Draw();
    TH1D * ratN1MCpm22SUB2 = (TH1D *) N1MCp22SUB2[erangePos][crange]->Clone(Form("ratN1MCpm22SUB2_%d_%d",cminCENT[crange],cmaxCENT[crange]));
    ratN1MCpm22SUB2->Divide(N1MCm22SUB2[erangeNeg][crange]);
    for (int i = 1; i<=N1MCp22SUB2[erangeNeg][crange]->GetNbinsX(); i++) {
        double x = N1MCp22SUB2[erangePos][crange]->GetBinContent(i);
        double xe = N1MCp22SUB2[erangePos][crange]->GetBinError(i);
        double y = N1MCm22SUB2[erangeNeg][crange]->GetBinContent(i);
        double ye = N1MCm22SUB2[erangeNeg][crange]->GetBinError(i);
        //ratN1MCpm22SUB2->SetBinError(i, ErrRatCalc(x, y, xe, ye));
    }
    ratN1MCpm22SUB2->SetMarkerColor(kBlack);
    ratN1MCpm22SUB2->SetLineColor(kBlack);
    ratN1MCpm22SUB2->SetMarkerStyle(20);
    ratN1MCpm22SUB2->SetMarkerSize(1.1);
    ratN1MCpm22SUB2->Draw("same");
    TPaveText * txratN1MCpm22SUB2 = new TPaveText(0.19, 0.84, 0.44, 0.93, "NDC");
    SetTPaveTxt(txratN1MCpm22SUB2, 18);
    txratN1MCpm22SUB2->AddText(Form("POI: %1.1f < |#eta| < %1.1f",eminETA[erangeNeg],emaxETA[erangeNeg]));
    txratN1MCpm22SUB2->AddText(Form("(%d-%d%%)",cminCENT[crange],cmaxCENT[crange]));
    txratN1MCpm22SUB2->Draw();
    TLine * lnratN1MCpm22SUB2 = new TLine(0.0, 1.0, 12.0, 1.0);
    lnratN1MCpm22SUB2->SetLineWidth(2);
    lnratN1MCpm22SUB2->Draw();
    cratN1MCpm22SUB2->Print(Form("results/Pt_Distributions/N1MCpm22SUB2_ratio_%d_%d.png",cminCENT[crange],cmaxCENT[crange]),"png");


}
