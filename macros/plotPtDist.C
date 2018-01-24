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
static const int netabins = 16;
static const double eminETA[] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0, -2.4,  0.0, -2.4,  0.4};
static const double emaxETA[] = {-2.0, -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0,  2.4,  0.0,  2.4, -0.4,  2.4};
static const double etaMid[] = {-2.2, -1.8, -1.4, -1.0, -0.6, -0.2,  0.2,  0.6,  1.0,  1.4,  1.8,  2.2};
TString etags[] = {"-24_-20", "-20_-16", "-16_-12", "-12_-8", "-8_-4", "-4_0", "0_4", "4_8", "8_12", "12_16", "16_20", "20_24",
                   "-24_0", "0_24", "-24_-4", "4_24"};
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
    if (!fopen("results/Pt_Distributions/N1MC22","r")) system("mkdir results/Pt_Distributions/N1MC22");
    if (!fopen("results/Pt_Distributions/N1MC18","r")) system("mkdir results/Pt_Distributions/N1MC18");

    int erangeNeg = 0;
    int erangePos = 0;
    TH1D * h0 = new TH1D("h0", "", 100, 0.0, 12.0);
    h0->SetStats(0);
    h0->SetXTitle("p_{T} (GeV/c)");
    h0->SetYTitle("v_{1}");
    h0->GetYaxis()->SetRangeUser(-0.03, 0.2);

    erangeNeg = 13;
    erangePos = 12;
/*    for (int cbin = 0; cbin<ncentbins; cbin++) {
        TCanvas * cN1MCpm22SUB2 = new TCanvas(Form("cN1MCpm22SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]),"cN1MCpm22SUB2",650,600);
        TPad * padN1MCpm22SUB2 = (TPad *) cN1MCpm22SUB2->cd();
        padN1MCpm22SUB2->SetGrid();
        TH1D * hN1MCpm22SUB2 = (TH1D *) h0->Clone(Form("hN1MCpm22SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCpm22SUB2->SetYTitle("v_{1}^{even}");
        hN1MCpm22SUB2->GetYaxis()->SetRangeUser(-0.03, 0.16);
        hN1MCpm22SUB2->Draw();
        N1MCm22SUB2[erangeNeg][cbin]->SetMarkerColor(kRed);
        N1MCm22SUB2[erangeNeg][cbin]->SetLineColor(kRed);
        N1MCm22SUB2[erangeNeg][cbin]->SetMarkerStyle(24);
        N1MCm22SUB2[erangeNeg][cbin]->SetMarkerSize(1.2);
        N1MCm22SUB2[erangeNeg][cbin]->Draw("same");
        N1MCp22SUB2[erangePos][cbin]->SetMarkerColor(kBlue);
        N1MCp22SUB2[erangePos][cbin]->SetLineColor(kBlue);
        N1MCp22SUB2[erangePos][cbin]->SetMarkerStyle(24);
        N1MCp22SUB2[erangePos][cbin]->SetMarkerSize(1.2);
        N1MCp22SUB2[erangePos][cbin]->Draw("same");
        TLegend * legN1MCpm22SUB2 = new TLegend(0.49, 0.17, 0.78, 0.27);
        SetLegend(legN1MCpm22SUB2, 18);
        legN1MCpm22SUB2->AddEntry(N1MCm22SUB2[erangeNeg][cbin],"#Psi_{1}^{trk}{2SUB} (-2.4 < #eta < -2.0)","p");
        legN1MCpm22SUB2->AddEntry(N1MCp22SUB2[erangePos][cbin],"#Psi_{1}^{trk}{2SUB} (2.0 < #eta < 2.4)","p");
        legN1MCpm22SUB2->Draw();
        TPaveText * txN1MCpm22SUB2 = new TPaveText(0.19, 0.84, 0.44, 0.93, "NDC");
        SetTPaveTxt(txN1MCpm22SUB2, 18);
        txN1MCpm22SUB2->AddText(Form("POI: %1.1f < |#eta| < %1.1f",eminETA[erangeNeg],emaxETA[erangeNeg]));
        txN1MCpm22SUB2->AddText(Form("(%d-%d%%)",cminCENT[cbin],cmaxCENT[cbin]));
        txN1MCpm22SUB2->Draw();
        cN1MCpm22SUB2->Print(Form("results/Pt_Distributions/N1MC22/N1MCpm22SUB2_%d_%d.png",cminCENT[cbin],cmaxCENT[cbin]),"png");
        cN1MCpm22SUB2->Close();


        TCanvas * cratN1MCpm22SUB2 = new TCanvas(Form("cratN1MCpm22SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]),"cratN1MCpm22SUB2",650,600);
        TPad * padratN1MCpm22SUB2 = (TPad *) cratN1MCpm22SUB2->cd();
        padratN1MCpm22SUB2->SetGrid();
        TH1D * hratN1MCpm22SUB2 = (TH1D *) h0->Clone(Form("hratN1MCpm22SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hratN1MCpm22SUB2->SetYTitle("v_{1}^{even}{+#eta} / v_{1}^{even}{-#eta}");
        hratN1MCpm22SUB2->GetYaxis()->SetRangeUser(0.6, 1.4);
        hratN1MCpm22SUB2->Draw();
        TH1D * ratN1MCpm22SUB2 = (TH1D *) N1MCp22SUB2[erangePos][cbin]->Clone(Form("ratN1MCpm22SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        ratN1MCpm22SUB2->Divide(N1MCm22SUB2[erangeNeg][cbin]);
        for (int i = 1; i<=N1MCp22SUB2[erangeNeg][cbin]->GetNbinsX(); i++) {
            double x = N1MCp22SUB2[erangePos][cbin]->GetBinContent(i);
            double xe = N1MCp22SUB2[erangePos][cbin]->GetBinError(i);
            double y = N1MCm22SUB2[erangeNeg][cbin]->GetBinContent(i);
            double ye = N1MCm22SUB2[erangeNeg][cbin]->GetBinError(i);
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
        txratN1MCpm22SUB2->AddText(Form("(%d-%d%%)",cminCENT[cbin],cmaxCENT[cbin]));
        txratN1MCpm22SUB2->Draw();
        TLine * lnratN1MCpm22SUB2 = new TLine(0.0, 1.0, 12.0, 1.0);
        lnratN1MCpm22SUB2->SetLineWidth(2);
        lnratN1MCpm22SUB2->Draw();
        cratN1MCpm22SUB2->Print(Form("results/Pt_Distributions/N1MC22/N1MCpm22SUB2_ratio_%d_%d.png",cminCENT[cbin],cmaxCENT[cbin]),"png");
        cratN1MCpm22SUB2->Close();


        TCanvas * cN1MCpm22SUB3 = new TCanvas(Form("cN1MCpm22SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]),"cN1MCpm22SUB3",650,600);
        TPad * padN1MCpm22SUB3 = (TPad *) cN1MCpm22SUB3->cd();
        padN1MCpm22SUB3->SetGrid();
        TH1D * hN1MCpm22SUB3 = (TH1D *) h0->Clone(Form("hN1MCpm22SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCpm22SUB3->SetYTitle("v_{1}^{even}");
        hN1MCpm22SUB3->GetYaxis()->SetRangeUser(-0.03, 0.16);
        hN1MCpm22SUB3->Draw();
        N1MCm22SUB3[erangeNeg][cbin]->SetMarkerColor(kRed);
        N1MCm22SUB3[erangeNeg][cbin]->SetLineColor(kRed);
        N1MCm22SUB3[erangeNeg][cbin]->SetMarkerStyle(24);
        N1MCm22SUB3[erangeNeg][cbin]->SetMarkerSize(1.2);
        N1MCm22SUB3[erangeNeg][cbin]->Draw("same");
        N1MCp22SUB3[erangePos][cbin]->SetMarkerColor(kBlue);
        N1MCp22SUB3[erangePos][cbin]->SetLineColor(kBlue);
        N1MCp22SUB3[erangePos][cbin]->SetMarkerStyle(24);
        N1MCp22SUB3[erangePos][cbin]->SetMarkerSize(1.2);
        N1MCp22SUB3[erangePos][cbin]->Draw("same");
        TLegend * legN1MCpm22SUB3 = new TLegend(0.49, 0.17, 0.78, 0.27);
        SetLegend(legN1MCpm22SUB3, 18);
        legN1MCpm22SUB3->AddEntry(N1MCm22SUB3[erangeNeg][cbin],"#Psi_{1}^{trk}{3SUB} (-2.4 < #eta < -2.0)","p");
        legN1MCpm22SUB3->AddEntry(N1MCp22SUB3[erangePos][cbin],"#Psi_{1}^{trk}{3SUB} (2.0 < #eta < 2.4)","p");
        legN1MCpm22SUB3->Draw();
        TPaveText * txN1MCpm22SUB3 = new TPaveText(0.19, 0.84, 0.44, 0.93, "NDC");
        SetTPaveTxt(txN1MCpm22SUB3, 18);
        txN1MCpm22SUB3->AddText(Form("POI: %1.1f < |#eta| < %1.1f",eminETA[erangeNeg],emaxETA[erangeNeg]));
        txN1MCpm22SUB3->AddText(Form("(%d-%d%%)",cminCENT[cbin],cmaxCENT[cbin]));
        txN1MCpm22SUB3->Draw();
        cN1MCpm22SUB3->Print(Form("results/Pt_Distributions/N1MC22/N1MCpm22SUB3_%d_%d.png",cminCENT[cbin],cmaxCENT[cbin]),"png");
        cN1MCpm22SUB3->Close();


        TCanvas * cratN1MCpm22SUB3 = new TCanvas(Form("cratN1MCpm22SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]),"cratN1MCpm22SUB3",650,600);
        TPad * padratN1MCpm22SUB3 = (TPad *) cratN1MCpm22SUB3->cd();
        padratN1MCpm22SUB3->SetGrid();
        TH1D * hratN1MCpm22SUB3 = (TH1D *) h0->Clone(Form("hratN1MCpm22SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hratN1MCpm22SUB3->SetYTitle("v_{1}^{even}{+#eta} / v_{1}^{even}{-#eta}");
        hratN1MCpm22SUB3->GetYaxis()->SetRangeUser(0.6, 1.4);
        hratN1MCpm22SUB3->Draw();
        TH1D * ratN1MCpm22SUB3 = (TH1D *) N1MCp22SUB3[erangePos][cbin]->Clone(Form("ratN1MCpm22SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        ratN1MCpm22SUB3->Divide(N1MCm22SUB3[erangeNeg][cbin]);
        for (int i = 1; i<=N1MCp22SUB3[erangeNeg][cbin]->GetNbinsX(); i++) {
            double x = N1MCp22SUB3[erangePos][cbin]->GetBinContent(i);
            double xe = N1MCp22SUB3[erangePos][cbin]->GetBinError(i);
            double y = N1MCm22SUB3[erangeNeg][cbin]->GetBinContent(i);
            double ye = N1MCm22SUB3[erangeNeg][cbin]->GetBinError(i);
            //ratN1MCpm22SUB3->SetBinError(i, ErrRatCalc(x, y, xe, ye));
        }
        ratN1MCpm22SUB3->SetMarkerColor(kBlack);
        ratN1MCpm22SUB3->SetLineColor(kBlack);
        ratN1MCpm22SUB3->SetMarkerStyle(20);
        ratN1MCpm22SUB3->SetMarkerSize(1.1);
        ratN1MCpm22SUB3->Draw("same");
        TPaveText * txratN1MCpm22SUB3 = new TPaveText(0.19, 0.84, 0.44, 0.93, "NDC");
        SetTPaveTxt(txratN1MCpm22SUB3, 18);
        txratN1MCpm22SUB3->AddText(Form("POI: %1.1f < |#eta| < %1.1f",eminETA[erangeNeg],emaxETA[erangeNeg]));
        txratN1MCpm22SUB3->AddText(Form("(%d-%d%%)",cminCENT[cbin],cmaxCENT[cbin]));
        txratN1MCpm22SUB3->Draw();
        TLine * lnratN1MCpm22SUB3 = new TLine(0.0, 1.0, 12.0, 1.0);
        lnratN1MCpm22SUB3->SetLineWidth(2);
        lnratN1MCpm22SUB3->Draw();
        cratN1MCpm22SUB3->Print(Form("results/Pt_Distributions/N1MC22/N1MCpm22SUB3_ratio_%d_%d.png",cminCENT[cbin],cmaxCENT[cbin]),"png");
        cratN1MCpm22SUB3->Close();


        // N1MC18

        erangeNeg = 15;
        erangePos = 14;

        TCanvas * cN1MCpm18SUB2 = new TCanvas(Form("cN1MCpm18SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]),"cN1MCpm18SUB2",650,600);
        TPad * padN1MCpm18SUB2 = (TPad *) cN1MCpm18SUB2->cd();
        padN1MCpm18SUB2->SetGrid();
        TH1D * hN1MCpm18SUB2 = (TH1D *) h0->Clone(Form("hN1MCpm18SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCpm18SUB2->SetYTitle("v_{1}^{even}");
        hN1MCpm18SUB2->GetYaxis()->SetRangeUser(-0.03, 0.16);
        hN1MCpm18SUB2->Draw();
        N1MCm18SUB2[erangeNeg][cbin]->SetMarkerColor(kRed);
        N1MCm18SUB2[erangeNeg][cbin]->SetLineColor(kRed);
        N1MCm18SUB2[erangeNeg][cbin]->SetMarkerStyle(24);
        N1MCm18SUB2[erangeNeg][cbin]->SetMarkerSize(1.2);
        N1MCm18SUB2[erangeNeg][cbin]->Draw("same");
        N1MCp18SUB2[erangePos][cbin]->SetMarkerColor(kBlue);
        N1MCp18SUB2[erangePos][cbin]->SetLineColor(kBlue);
        N1MCp18SUB2[erangePos][cbin]->SetMarkerStyle(24);
        N1MCp18SUB2[erangePos][cbin]->SetMarkerSize(1.2);
        N1MCp18SUB2[erangePos][cbin]->Draw("same");
        TLegend * legN1MCpm18SUB2 = new TLegend(0.49, 0.17, 0.78, 0.27);
        SetLegend(legN1MCpm18SUB2, 18);
        legN1MCpm18SUB2->AddEntry(N1MCm18SUB2[erangeNeg][cbin],"#Psi_{1}^{trk}{2SUB} (-2.0 < #eta < -1.6)","p");
        legN1MCpm18SUB2->AddEntry(N1MCp18SUB2[erangePos][cbin],"#Psi_{1}^{trk}{2SUB} (1.6 < #eta < 2.0)","p");
        legN1MCpm18SUB2->Draw();
        TPaveText * txN1MCpm18SUB2 = new TPaveText(0.19, 0.84, 0.44, 0.93, "NDC");
        SetTPaveTxt(txN1MCpm18SUB2, 18);
        txN1MCpm18SUB2->AddText(Form("POI: %1.1f < |#eta| < %1.1f",eminETA[erangeNeg],emaxETA[erangeNeg]));
        txN1MCpm18SUB2->AddText(Form("(%d-%d%%)",cminCENT[cbin],cmaxCENT[cbin]));
        txN1MCpm18SUB2->Draw();
        cN1MCpm18SUB2->Print(Form("results/Pt_Distributions/N1MC18/N1MCpm18SUB2_%d_%d.png",cminCENT[cbin],cmaxCENT[cbin]),"png");
        cN1MCpm18SUB2->Close();


        TCanvas * cratN1MCpm18SUB2 = new TCanvas(Form("cratN1MCpm18SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]),"cratN1MCpm18SUB2",650,600);
        TPad * padratN1MCpm18SUB2 = (TPad *) cratN1MCpm18SUB2->cd();
        padratN1MCpm18SUB2->SetGrid();
        TH1D * hratN1MCpm18SUB2 = (TH1D *) h0->Clone(Form("hratN1MCpm18SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hratN1MCpm18SUB2->SetYTitle("v_{1}^{even}{+#eta} / v_{1}^{even}{-#eta}");
        hratN1MCpm18SUB2->GetYaxis()->SetRangeUser(0.6, 1.4);
        hratN1MCpm18SUB2->Draw();
        TH1D * ratN1MCpm18SUB2 = (TH1D *) N1MCp18SUB2[erangePos][cbin]->Clone(Form("ratN1MCpm18SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        ratN1MCpm18SUB2->Divide(N1MCm18SUB2[erangeNeg][cbin]);
        for (int i = 1; i<=N1MCp18SUB2[erangeNeg][cbin]->GetNbinsX(); i++) {
            double x = N1MCp18SUB2[erangePos][cbin]->GetBinContent(i);
            double xe = N1MCp18SUB2[erangePos][cbin]->GetBinError(i);
            double y = N1MCm18SUB2[erangeNeg][cbin]->GetBinContent(i);
            double ye = N1MCm18SUB2[erangeNeg][cbin]->GetBinError(i);
            //ratN1MCpm18SUB2->SetBinError(i, ErrRatCalc(x, y, xe, ye));
        }
        ratN1MCpm18SUB2->SetMarkerColor(kBlack);
        ratN1MCpm18SUB2->SetLineColor(kBlack);
        ratN1MCpm18SUB2->SetMarkerStyle(20);
        ratN1MCpm18SUB2->SetMarkerSize(1.1);
        ratN1MCpm18SUB2->Draw("same");
        TPaveText * txratN1MCpm18SUB2 = new TPaveText(0.19, 0.84, 0.44, 0.93, "NDC");
        SetTPaveTxt(txratN1MCpm18SUB2, 18);
        txratN1MCpm18SUB2->AddText(Form("POI: %1.1f < |#eta| < %1.1f",eminETA[erangeNeg],emaxETA[erangeNeg]));
        txratN1MCpm18SUB2->AddText(Form("(%d-%d%%)",cminCENT[cbin],cmaxCENT[cbin]));
        txratN1MCpm18SUB2->Draw();
        TLine * lnratN1MCpm18SUB2 = new TLine(0.0, 1.0, 12.0, 1.0);
        lnratN1MCpm18SUB2->SetLineWidth(2);
        lnratN1MCpm18SUB2->Draw();
        cratN1MCpm18SUB2->Print(Form("results/Pt_Distributions/N1MC18/N1MCpm18SUB2_ratio_%d_%d.png",cminCENT[cbin],cmaxCENT[cbin]),"png");
        cratN1MCpm18SUB2->Close();


        TCanvas * cN1MCpm18SUB3 = new TCanvas(Form("cN1MCpm18SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]),"cN1MCpm18SUB3",650,600);
        TPad * padN1MCpm18SUB3 = (TPad *) cN1MCpm18SUB3->cd();
        padN1MCpm18SUB3->SetGrid();
        TH1D * hN1MCpm18SUB3 = (TH1D *) h0->Clone(Form("hN1MCpm18SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hN1MCpm18SUB3->SetYTitle("v_{1}^{even}");
        hN1MCpm18SUB3->GetYaxis()->SetRangeUser(-0.03, 0.16);
        hN1MCpm18SUB3->Draw();
        N1MCm18SUB3[erangeNeg][cbin]->SetMarkerColor(kRed);
        N1MCm18SUB3[erangeNeg][cbin]->SetLineColor(kRed);
        N1MCm18SUB3[erangeNeg][cbin]->SetMarkerStyle(24);
        N1MCm18SUB3[erangeNeg][cbin]->SetMarkerSize(1.2);
        N1MCm18SUB3[erangeNeg][cbin]->Draw("same");
        N1MCp18SUB3[erangePos][cbin]->SetMarkerColor(kBlue);
        N1MCp18SUB3[erangePos][cbin]->SetLineColor(kBlue);
        N1MCp18SUB3[erangePos][cbin]->SetMarkerStyle(24);
        N1MCp18SUB3[erangePos][cbin]->SetMarkerSize(1.2);
        N1MCp18SUB3[erangePos][cbin]->Draw("same");
        TLegend * legN1MCpm18SUB3 = new TLegend(0.49, 0.17, 0.78, 0.27);
        SetLegend(legN1MCpm18SUB3, 18);
        legN1MCpm18SUB3->AddEntry(N1MCm18SUB3[erangeNeg][cbin],"#Psi_{1}^{trk}{3SUB} (-2.0 < #eta < -1.6)","p");
        legN1MCpm18SUB3->AddEntry(N1MCp18SUB3[erangePos][cbin],"#Psi_{1}^{trk}{3SUB} (1.6 < #eta < 2.0)","p");
        legN1MCpm18SUB3->Draw();
        TPaveText * txN1MCpm18SUB3 = new TPaveText(0.19, 0.84, 0.44, 0.93, "NDC");
        SetTPaveTxt(txN1MCpm18SUB3, 18);
        txN1MCpm18SUB3->AddText(Form("POI: %1.1f < |#eta| < %1.1f",eminETA[erangeNeg],emaxETA[erangeNeg]));
        txN1MCpm18SUB3->AddText(Form("(%d-%d%%)",cminCENT[cbin],cmaxCENT[cbin]));
        txN1MCpm18SUB3->Draw();
        cN1MCpm18SUB3->Print(Form("results/Pt_Distributions/N1MC18/N1MCpm18SUB3_%d_%d.png",cminCENT[cbin],cmaxCENT[cbin]),"png");
        cN1MCpm18SUB3->Close();


        TCanvas * cratN1MCpm18SUB3 = new TCanvas(Form("cratN1MCpm18SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]),"cratN1MCpm18SUB3",650,600);
        TPad * padratN1MCpm18SUB3 = (TPad *) cratN1MCpm18SUB3->cd();
        padratN1MCpm18SUB3->SetGrid();
        TH1D * hratN1MCpm18SUB3 = (TH1D *) h0->Clone(Form("hratN1MCpm18SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        hratN1MCpm18SUB3->SetYTitle("v_{1}^{even}{+#eta} / v_{1}^{even}{-#eta}");
        hratN1MCpm18SUB3->GetYaxis()->SetRangeUser(0.6, 1.4);
        hratN1MCpm18SUB3->Draw();
        TH1D * ratN1MCpm18SUB3 = (TH1D *) N1MCp18SUB3[erangePos][cbin]->Clone(Form("ratN1MCpm18SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        ratN1MCpm18SUB3->Divide(N1MCm18SUB3[erangeNeg][cbin]);
        for (int i = 1; i<=N1MCp18SUB3[erangeNeg][cbin]->GetNbinsX(); i++) {
            double x = N1MCp18SUB3[erangePos][cbin]->GetBinContent(i);
            double xe = N1MCp18SUB3[erangePos][cbin]->GetBinError(i);
            double y = N1MCm18SUB3[erangeNeg][cbin]->GetBinContent(i);
            double ye = N1MCm18SUB3[erangeNeg][cbin]->GetBinError(i);
            //ratN1MCpm18SUB3->SetBinError(i, ErrRatCalc(x, y, xe, ye));
        }
        ratN1MCpm18SUB3->SetMarkerColor(kBlack);
        ratN1MCpm18SUB3->SetLineColor(kBlack);
        ratN1MCpm18SUB3->SetMarkerStyle(20);
        ratN1MCpm18SUB3->SetMarkerSize(1.1);
        ratN1MCpm18SUB3->Draw("same");
        TPaveText * txratN1MCpm18SUB3 = new TPaveText(0.19, 0.84, 0.44, 0.93, "NDC");
        SetTPaveTxt(txratN1MCpm18SUB3, 18);
        txratN1MCpm18SUB3->AddText(Form("POI: %1.1f < |#eta| < %1.1f",eminETA[erangeNeg],emaxETA[erangeNeg]));
        txratN1MCpm18SUB3->AddText(Form("(%d-%d%%)",cminCENT[cbin],cmaxCENT[cbin]));
        txratN1MCpm18SUB3->Draw();
        TLine * lnratN1MCpm18SUB3 = new TLine(0.0, 1.0, 12.0, 1.0);
        lnratN1MCpm18SUB3->SetLineWidth(2);
        lnratN1MCpm18SUB3->Draw();
        cratN1MCpm18SUB3->Print(Form("results/Pt_Distributions/N1MC18/N1MCpm18SUB3_ratio_%d_%d.png",cminCENT[cbin],cmaxCENT[cbin]),"png");
        cratN1MCpm18SUB3->Close();
    }*/

    // cent scan (N1MC22)
    erangeNeg = 13;
    erangePos = 12;
    TCanvas * cN1MC22SUB2_centscan = new TCanvas("cN1MC22SUB2_centscan","cN1MC22SUB2_centscan",1100,600);
    cN1MC22SUB2_centscan->Divide(4,2,0,0);
    for (int cbin = 0; cbin<8; cbin++) {
        TPad * padN1MC22SUB2_centscan = (TPad *) cN1MC22SUB2_centscan->cd(cbin+1);
        padN1MC22SUB2_centscan->SetGrid();
        TH1D * hN1MC22SUB2_centscan = (TH1D *) h0->Clone(Form("hN1MC22SUB2_centscan_%d",cbin));
        hN1MC22SUB2_centscan->SetYTitle("v_{1}^{even}");
        hN1MC22SUB2_centscan->GetYaxis()->SetRangeUser(-0.03, 0.16);
        hN1MC22SUB2_centscan->Draw();
        N1MCm22SUB2[erangeNeg][cbin]->SetMarkerColor(kRed);
        N1MCm22SUB2[erangeNeg][cbin]->SetLineColor(kRed);
        N1MCm22SUB2[erangeNeg][cbin]->SetMarkerStyle(24);
        N1MCm22SUB2[erangeNeg][cbin]->SetMarkerSize(1.2);
        N1MCm22SUB2[erangeNeg][cbin]->Draw("same");
        N1MCp22SUB2[erangePos][cbin]->SetMarkerColor(kBlue);
        N1MCp22SUB2[erangePos][cbin]->SetLineColor(kBlue);
        N1MCp22SUB2[erangePos][cbin]->SetMarkerStyle(24);
        N1MCp22SUB2[erangePos][cbin]->SetMarkerSize(1.2);
        N1MCp22SUB2[erangePos][cbin]->Draw("same");
        TPaveText * txN1MC22SUB2_centscan_c = new TPaveText(0.75, 0.89, 0.89, 0.98, "NDC");
        SetTPaveTxt(txN1MC22SUB2_centscan_c, 16);
        txN1MC22SUB2_centscan_c->AddText(Form("%d-%d%%",cminCENT[cbin],cmaxCENT[cbin]));
        txN1MC22SUB2_centscan_c->Draw();
    }
    cN1MC22SUB2_centscan->cd(2);
    TLegend * legN1MC22SUB2_centscan = new TLegend(0.05, 0.79, 0.39, 0.97);
    SetLegend(legN1MC22SUB2_centscan, 16);
    legN1MC22SUB2_centscan->AddEntry(N1MCm22SUB2[erangeNeg][0],"#Psi_{1}^{trk}{2SUB}: -2.4<#eta<-2.0","p");
    legN1MC22SUB2_centscan->AddEntry(N1MCp22SUB2[erangePos][0],"#Psi_{1}^{trk}{2SUB}: 2.0<#eta<2.4","p");
    legN1MC22SUB2_centscan->Draw();
    cN1MC22SUB2_centscan->cd(1);
    TPaveText * txN1MC22SUB2_centscan = new TPaveText(0.24, 0.77, 0.73, 0.97, "NDC");
    SetTPaveTxt(txN1MC22SUB2_centscan, 16);
    txN1MC22SUB2_centscan->AddText("PbPb #sqrt{s_{NN}}=5.02 TeV");
    txN1MC22SUB2_centscan->AddText(Form("POI: %1.1f<|#eta|<%1.1f",eminETA[erangeNeg],emaxETA[erangeNeg]));
    txN1MC22SUB2_centscan->Draw();
    cN1MC22SUB2_centscan->Print("results/Pt_Distributions/N1MC22/centscan_N1MCpm22SUB2.png","png");


    TCanvas * cratN1MC22SUB2_centscan = new TCanvas("cratN1MC22SUB2_centscan","cratN1MC22SUB2_centscan",1100,600);
    cratN1MC22SUB2_centscan->Divide(4,2,0,0);
    for (int cbin = 0; cbin<8; cbin++) {
        TPad * padratN1MC22SUB2_centscan = (TPad *) cratN1MC22SUB2_centscan->cd(cbin+1);
        padratN1MC22SUB2_centscan->SetGrid();
        TH1D * hratN1MC22SUB2_centscan = (TH1D *) h0->Clone(Form("hratN1MC22SUB2_centscan_%d",cbin));
        hratN1MC22SUB2_centscan->SetYTitle("v_{1}^{even}{+#eta} / v_{1}^{even}{-#eta}");
        hratN1MC22SUB2_centscan->GetYaxis()->SetRangeUser(0.6, 1.4);
        hratN1MC22SUB2_centscan->Draw();
        TH1D * ratN1MC22SUB2_centscan = (TH1D *) N1MCp22SUB2[erangePos][cbin]->Clone(Form("ratN1MCpm22SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        ratN1MC22SUB2_centscan->Divide(N1MCm22SUB2[erangeNeg][cbin]);
        ratN1MC22SUB2_centscan->SetMarkerColor(kBlack);
        ratN1MC22SUB2_centscan->SetLineColor(kBlack);
        ratN1MC22SUB2_centscan->SetMarkerStyle(20);
        ratN1MC22SUB2_centscan->SetMarkerSize(1.1);
        ratN1MC22SUB2_centscan->Draw("same");
        TPaveText * txratN1MC22SUB2_centscan_c = new TPaveText(0.75, 0.89, 0.89, 0.98, "NDC");
        SetTPaveTxt(txratN1MC22SUB2_centscan_c, 16);
        txratN1MC22SUB2_centscan_c->AddText(Form("%d-%d%%",cminCENT[cbin],cmaxCENT[cbin]));
        txratN1MC22SUB2_centscan_c->Draw();
        TLine * lnratN1MC22SUB2_centscan = new TLine(0.0, 1.0, 12.0, 1.0);
        lnratN1MC22SUB2_centscan->SetLineWidth(2);
        lnratN1MC22SUB2_centscan->Draw();
    }
    cratN1MC22SUB2_centscan->cd(1);
    TPaveText * txratN1MC22SUB2_centscan = new TPaveText(0.24, 0.77, 0.73, 0.97, "NDC");
    SetTPaveTxt(txratN1MC22SUB2_centscan, 16);
    txratN1MC22SUB2_centscan->AddText("PbPb #sqrt{s_{NN}}=5.02 TeV");
    txratN1MC22SUB2_centscan->AddText(Form("POI: %1.1f<|#eta|<%1.1f",eminETA[erangeNeg],emaxETA[erangeNeg]));
    txratN1MC22SUB2_centscan->Draw();
    cN1MC22SUB2_centscan->Print("results/Pt_Distributions/N1MC22/centscan_N1MCpm22SUB2_ratio.png","png");


    TCanvas * cN1MC22SUB3_centscan = new TCanvas("cN1MC22SUB3_centscan","cN1MC22SUB3_centscan",1100,600);
    cN1MC22SUB3_centscan->Divide(4,2,0,0);
    for (int cbin = 0; cbin<8; cbin++) {
        TPad * padN1MC22SUB3_centscan = (TPad *) cN1MC22SUB3_centscan->cd(cbin+1);
        padN1MC22SUB3_centscan->SetGrid();
        TH1D * hN1MC22SUB3_centscan = (TH1D *) h0->Clone(Form("hN1MC22SUB3_centscan_%d",cbin));
        hN1MC22SUB3_centscan->SetYTitle("v_{1}^{even}");
        hN1MC22SUB3_centscan->GetYaxis()->SetRangeUser(-0.03, 0.16);
        hN1MC22SUB3_centscan->Draw();
        N1MCm22SUB3[erangeNeg][cbin]->SetMarkerColor(kRed);
        N1MCm22SUB3[erangeNeg][cbin]->SetLineColor(kRed);
        N1MCm22SUB3[erangeNeg][cbin]->SetMarkerStyle(24);
        N1MCm22SUB3[erangeNeg][cbin]->SetMarkerSize(1.2);
        N1MCm22SUB3[erangeNeg][cbin]->Draw("same");
        N1MCp22SUB3[erangePos][cbin]->SetMarkerColor(kBlue);
        N1MCp22SUB3[erangePos][cbin]->SetLineColor(kBlue);
        N1MCp22SUB3[erangePos][cbin]->SetMarkerStyle(24);
        N1MCp22SUB3[erangePos][cbin]->SetMarkerSize(1.2);
        N1MCp22SUB3[erangePos][cbin]->Draw("same");
        TPaveText * txN1MC22SUB3_centscan_c = new TPaveText(0.75, 0.89, 0.89, 0.98, "NDC");
        SetTPaveTxt(txN1MC22SUB3_centscan_c, 16);
        txN1MC22SUB3_centscan_c->AddText(Form("%d-%d%%",cminCENT[cbin],cmaxCENT[cbin]));
        txN1MC22SUB3_centscan_c->Draw();
    }
    cN1MC22SUB3_centscan->cd(2);
    TLegend * legN1MC22SUB3_centscan = new TLegend(0.05, 0.79, 0.39, 0.97);
    SetLegend(legN1MC22SUB3_centscan, 16);
    legN1MC22SUB3_centscan->AddEntry(N1MCm22SUB3[erangeNeg][0],"#Psi_{1}^{trk}{3SUB}: -2.4<#eta<-2.0","p");
    legN1MC22SUB3_centscan->AddEntry(N1MCp22SUB3[erangePos][0],"#Psi_{1}^{trk}{3SUB}: 2.0<#eta<2.4","p");
    legN1MC22SUB3_centscan->Draw();
    cN1MC22SUB3_centscan->cd(1);
    TPaveText * txN1MC22SUB3_centscan = new TPaveText(0.24, 0.77, 0.73, 0.97, "NDC");
    SetTPaveTxt(txN1MC22SUB3_centscan, 16);
    txN1MC22SUB3_centscan->AddText("PbPb #sqrt{s_{NN}}=5.02 TeV");
    txN1MC22SUB3_centscan->AddText(Form("POI: %1.1f<|#eta|<%1.1f",eminETA[erangeNeg],emaxETA[erangeNeg]));
    txN1MC22SUB3_centscan->Draw();
    cN1MC22SUB3_centscan->Print("results/Pt_Distributions/N1MC22/centscan_N1MCpm22SUB3.png","png");


    TCanvas * cratN1MC22SUB3_centscan = new TCanvas("cratN1MC22SUB3_centscan","cratN1MC22SUB3_centscan",1100,600);
    cratN1MC22SUB3_centscan->Divide(4,2,0,0);
    for (int cbin = 0; cbin<8; cbin++) {
        TPad * padratN1MC22SUB3_centscan = (TPad *) cratN1MC22SUB3_centscan->cd(cbin+1);
        padratN1MC22SUB3_centscan->SetGrid();
        TH1D * hratN1MC22SUB3_centscan = (TH1D *) h0->Clone(Form("hratN1MC22SUB3_centscan_%d",cbin));
        hratN1MC22SUB3_centscan->SetYTitle("v_{1}^{even}{+#eta} / v_{1}^{even}{-#eta}");
        hratN1MC22SUB3_centscan->GetYaxis()->SetRangeUser(0.6, 1.4);
        hratN1MC22SUB3_centscan->Draw();
        TH1D * ratN1MC22SUB3_centscan = (TH1D *) N1MCp22SUB3[erangePos][cbin]->Clone(Form("ratN1MCpm22SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        ratN1MC22SUB3_centscan->Divide(N1MCm22SUB3[erangeNeg][cbin]);
        ratN1MC22SUB3_centscan->SetMarkerColor(kBlack);
        ratN1MC22SUB3_centscan->SetLineColor(kBlack);
        ratN1MC22SUB3_centscan->SetMarkerStyle(20);
        ratN1MC22SUB3_centscan->SetMarkerSize(1.1);
        ratN1MC22SUB3_centscan->Draw("same");
        TPaveText * txratN1MC22SUB3_centscan_c = new TPaveText(0.75, 0.89, 0.89, 0.98, "NDC");
        SetTPaveTxt(txratN1MC22SUB3_centscan_c, 16);
        txratN1MC22SUB3_centscan_c->AddText(Form("%d-%d%%",cminCENT[cbin],cmaxCENT[cbin]));
        txratN1MC22SUB3_centscan_c->Draw();
        TLine * lnratN1MC22SUB3_centscan = new TLine(0.0, 1.0, 12.0, 1.0);
        lnratN1MC22SUB3_centscan->SetLineWidth(2);
        lnratN1MC22SUB3_centscan->Draw();
    }
    cratN1MC22SUB3_centscan->cd(1);
    TPaveText * txratN1MC22SUB3_centscan = new TPaveText(0.24, 0.77, 0.73, 0.97, "NDC");
    SetTPaveTxt(txratN1MC22SUB3_centscan, 16);
    txratN1MC22SUB3_centscan->AddText("PbPb #sqrt{s_{NN}}=5.02 TeV");
    txratN1MC22SUB3_centscan->AddText(Form("POI: %1.1f<|#eta|<%1.1f",eminETA[erangeNeg],emaxETA[erangeNeg]));
    txratN1MC22SUB3_centscan->Draw();
    cN1MC22SUB3_centscan->Print("results/Pt_Distributions/N1MC22/centscan_N1MCpm22SUB3_ratio.png","png");


    // cent scan (N1MC18)
    erangeNeg = 15;
    erangePos = 14;
    TCanvas * cN1MC18SUB2_centscan = new TCanvas("cN1MC18SUB2_centscan","cN1MC18SUB2_centscan",1100,600);
    cN1MC18SUB2_centscan->Divide(4,2,0,0);
    for (int cbin = 0; cbin<8; cbin++) {
        TPad * padN1MC18SUB2_centscan = (TPad *) cN1MC18SUB2_centscan->cd(cbin+1);
        padN1MC18SUB2_centscan->SetGrid();
        TH1D * hN1MC18SUB2_centscan = (TH1D *) h0->Clone(Form("hN1MC18SUB2_centscan_%d",cbin));
        hN1MC18SUB2_centscan->SetYTitle("v_{1}^{even}");
        hN1MC18SUB2_centscan->GetYaxis()->SetRangeUser(-0.03, 0.16);
        hN1MC18SUB2_centscan->Draw();
        N1MCm18SUB2[erangeNeg][cbin]->SetMarkerColor(kRed);
        N1MCm18SUB2[erangeNeg][cbin]->SetLineColor(kRed);
        N1MCm18SUB2[erangeNeg][cbin]->SetMarkerStyle(24);
        N1MCm18SUB2[erangeNeg][cbin]->SetMarkerSize(1.2);
        N1MCm18SUB2[erangeNeg][cbin]->Draw("same");
        N1MCp18SUB2[erangePos][cbin]->SetMarkerColor(kBlue);
        N1MCp18SUB2[erangePos][cbin]->SetLineColor(kBlue);
        N1MCp18SUB2[erangePos][cbin]->SetMarkerStyle(24);
        N1MCp18SUB2[erangePos][cbin]->SetMarkerSize(1.2);
        N1MCp18SUB2[erangePos][cbin]->Draw("same");
        TPaveText * txN1MC18SUB2_centscan_c = new TPaveText(0.75, 0.89, 0.89, 0.98, "NDC");
        SetTPaveTxt(txN1MC18SUB2_centscan_c, 16);
        txN1MC18SUB2_centscan_c->AddText(Form("%d-%d%%",cminCENT[cbin],cmaxCENT[cbin]));
        txN1MC18SUB2_centscan_c->Draw();
    }
    cN1MC18SUB2_centscan->cd(2);
    TLegend * legN1MC18SUB2_centscan = new TLegend(0.05, 0.79, 0.39, 0.97);
    SetLegend(legN1MC18SUB2_centscan, 16);
    legN1MC18SUB2_centscan->AddEntry(N1MCm18SUB2[erangeNeg][0],"#Psi_{1}^{trk}{2SUB}: -2.0<#eta<-1.6","p");
    legN1MC18SUB2_centscan->AddEntry(N1MCp18SUB2[erangePos][0],"#Psi_{1}^{trk}{2SUB}: 1.6<#eta<2.0","p");
    legN1MC18SUB2_centscan->Draw();
    cN1MC18SUB2_centscan->cd(1);
    TPaveText * txN1MC18SUB2_centscan = new TPaveText(0.24, 0.77, 0.73, 0.97, "NDC");
    SetTPaveTxt(txN1MC18SUB2_centscan, 16);
    txN1MC18SUB2_centscan->AddText("PbPb #sqrt{s_{NN}}=5.02 TeV");
    txN1MC18SUB2_centscan->AddText(Form("POI: %1.1f<|#eta|<%1.1f",eminETA[erangeNeg],emaxETA[erangeNeg]));
    txN1MC18SUB2_centscan->Draw();
    cN1MC18SUB2_centscan->Print("results/Pt_Distributions/N1MC18/centscan_N1MCpm18SUB2.png","png");


    TCanvas * cratN1MC18SUB2_centscan = new TCanvas("cratN1MC18SUB2_centscan","cratN1MC18SUB2_centscan",1100,600);
    cratN1MC18SUB2_centscan->Divide(4,2,0,0);
    for (int cbin = 0; cbin<8; cbin++) {
        TPad * padratN1MC18SUB2_centscan = (TPad *) cratN1MC18SUB2_centscan->cd(cbin+1);
        padratN1MC18SUB2_centscan->SetGrid();
        TH1D * hratN1MC18SUB2_centscan = (TH1D *) h0->Clone(Form("hratN1MC18SUB2_centscan_%d",cbin));
        hratN1MC18SUB2_centscan->SetYTitle("v_{1}^{even}{+#eta} / v_{1}^{even}{-#eta}");
        hratN1MC18SUB2_centscan->GetYaxis()->SetRangeUser(0.6, 1.4);
        hratN1MC18SUB2_centscan->Draw();
        TH1D * ratN1MC18SUB2_centscan = (TH1D *) N1MCp18SUB2[erangePos][cbin]->Clone(Form("ratN1MCpm18SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        ratN1MC18SUB2_centscan->Divide(N1MCm18SUB2[erangeNeg][cbin]);
        ratN1MC18SUB2_centscan->SetMarkerColor(kBlack);
        ratN1MC18SUB2_centscan->SetLineColor(kBlack);
        ratN1MC18SUB2_centscan->SetMarkerStyle(20);
        ratN1MC18SUB2_centscan->SetMarkerSize(1.1);
        ratN1MC18SUB2_centscan->Draw("same");
        TPaveText * txratN1MC18SUB2_centscan_c = new TPaveText(0.75, 0.89, 0.89, 0.98, "NDC");
        SetTPaveTxt(txratN1MC18SUB2_centscan_c, 16);
        txratN1MC18SUB2_centscan_c->AddText(Form("%d-%d%%",cminCENT[cbin],cmaxCENT[cbin]));
        txratN1MC18SUB2_centscan_c->Draw();
        TLine * lnratN1MC18SUB2_centscan = new TLine(0.0, 1.0, 12.0, 1.0);
        lnratN1MC18SUB2_centscan->SetLineWidth(2);
        lnratN1MC18SUB2_centscan->Draw();
    }
    cratN1MC18SUB2_centscan->cd(1);
    TPaveText * txratN1MC18SUB2_centscan = new TPaveText(0.24, 0.77, 0.73, 0.97, "NDC");
    SetTPaveTxt(txratN1MC18SUB2_centscan, 16);
    txratN1MC18SUB2_centscan->AddText("PbPb #sqrt{s_{NN}}=5.02 TeV");
    txratN1MC18SUB2_centscan->AddText(Form("POI: %1.1f<|#eta|<%1.1f",eminETA[erangeNeg],emaxETA[erangeNeg]));
    txratN1MC18SUB2_centscan->Draw();
    cN1MC18SUB2_centscan->Print("results/Pt_Distributions/N1MC18/centscan_N1MCpm18SUB2_ratio.png","png");


    TCanvas * cN1MC18SUB3_centscan = new TCanvas("cN1MC18SUB3_centscan","cN1MC18SUB3_centscan",1100,600);
    cN1MC18SUB3_centscan->Divide(4,2,0,0);
    for (int cbin = 0; cbin<8; cbin++) {
        TPad * padN1MC18SUB3_centscan = (TPad *) cN1MC18SUB3_centscan->cd(cbin+1);
        padN1MC18SUB3_centscan->SetGrid();
        TH1D * hN1MC18SUB3_centscan = (TH1D *) h0->Clone(Form("hN1MC18SUB3_centscan_%d",cbin));
        hN1MC18SUB3_centscan->SetYTitle("v_{1}^{even}");
        hN1MC18SUB3_centscan->GetYaxis()->SetRangeUser(-0.03, 0.16);
        hN1MC18SUB3_centscan->Draw();
        N1MCm18SUB3[erangeNeg][cbin]->SetMarkerColor(kRed);
        N1MCm18SUB3[erangeNeg][cbin]->SetLineColor(kRed);
        N1MCm18SUB3[erangeNeg][cbin]->SetMarkerStyle(24);
        N1MCm18SUB3[erangeNeg][cbin]->SetMarkerSize(1.2);
        N1MCm18SUB3[erangeNeg][cbin]->Draw("same");
        N1MCp18SUB3[erangePos][cbin]->SetMarkerColor(kBlue);
        N1MCp18SUB3[erangePos][cbin]->SetLineColor(kBlue);
        N1MCp18SUB3[erangePos][cbin]->SetMarkerStyle(24);
        N1MCp18SUB3[erangePos][cbin]->SetMarkerSize(1.2);
        N1MCp18SUB3[erangePos][cbin]->Draw("same");
        TPaveText * txN1MC18SUB3_centscan_c = new TPaveText(0.75, 0.89, 0.89, 0.98, "NDC");
        SetTPaveTxt(txN1MC18SUB3_centscan_c, 16);
        txN1MC18SUB3_centscan_c->AddText(Form("%d-%d%%",cminCENT[cbin],cmaxCENT[cbin]));
        txN1MC18SUB3_centscan_c->Draw();
    }
    cN1MC18SUB3_centscan->cd(2);
    TLegend * legN1MC18SUB3_centscan = new TLegend(0.05, 0.79, 0.39, 0.97);
    SetLegend(legN1MC18SUB3_centscan, 16);
    legN1MC18SUB3_centscan->AddEntry(N1MCm18SUB3[erangeNeg][0],"#Psi_{1}^{trk}{3SUB}: -2.0<#eta<-1.6","p");
    legN1MC18SUB3_centscan->AddEntry(N1MCp18SUB3[erangePos][0],"#Psi_{1}^{trk}{3SUB}: 1.6<#eta<2.0","p");
    legN1MC18SUB3_centscan->Draw();
    cN1MC18SUB3_centscan->cd(1);
    TPaveText * txN1MC18SUB3_centscan = new TPaveText(0.24, 0.77, 0.73, 0.97, "NDC");
    SetTPaveTxt(txN1MC18SUB3_centscan, 16);
    txN1MC18SUB3_centscan->AddText("PbPb #sqrt{s_{NN}}=5.02 TeV");
    txN1MC18SUB3_centscan->AddText(Form("POI: %1.1f<|#eta|<%1.1f",eminETA[erangeNeg],emaxETA[erangeNeg]));
    txN1MC18SUB3_centscan->Draw();
    cN1MC18SUB3_centscan->Print("results/Pt_Distributions/N1MC18/centscan_N1MCpm18SUB3.png","png");


    TCanvas * cratN1MC18SUB3_centscan = new TCanvas("cratN1MC18SUB3_centscan","cratN1MC18SUB3_centscan",1100,600);
    cratN1MC18SUB3_centscan->Divide(4,2,0,0);
    for (int cbin = 0; cbin<8; cbin++) {
        TPad * padratN1MC18SUB3_centscan = (TPad *) cratN1MC18SUB3_centscan->cd(cbin+1);
        padratN1MC18SUB3_centscan->SetGrid();
        TH1D * hratN1MC18SUB3_centscan = (TH1D *) h0->Clone(Form("hratN1MC18SUB3_centscan_%d",cbin));
        hratN1MC18SUB3_centscan->SetYTitle("v_{1}^{even}{+#eta} / v_{1}^{even}{-#eta}");
        hratN1MC18SUB3_centscan->GetYaxis()->SetRangeUser(0.6, 1.4);
        hratN1MC18SUB3_centscan->Draw();
        TH1D * ratN1MC18SUB3_centscan = (TH1D *) N1MCp18SUB3[erangePos][cbin]->Clone(Form("ratN1MCpm18SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]));
        ratN1MC18SUB3_centscan->Divide(N1MCm18SUB3[erangeNeg][cbin]);
        ratN1MC18SUB3_centscan->SetMarkerColor(kBlack);
        ratN1MC18SUB3_centscan->SetLineColor(kBlack);
        ratN1MC18SUB3_centscan->SetMarkerStyle(20);
        ratN1MC18SUB3_centscan->SetMarkerSize(1.1);
        ratN1MC18SUB3_centscan->Draw("same");
        TPaveText * txratN1MC18SUB3_centscan_c = new TPaveText(0.75, 0.89, 0.89, 0.98, "NDC");
        SetTPaveTxt(txratN1MC18SUB3_centscan_c, 16);
        txratN1MC18SUB3_centscan_c->AddText(Form("%d-%d%%",cminCENT[cbin],cmaxCENT[cbin]));
        txratN1MC18SUB3_centscan_c->Draw();
        TLine * lnratN1MC18SUB3_centscan = new TLine(0.0, 1.0, 12.0, 1.0);
        lnratN1MC18SUB3_centscan->SetLineWidth(2);
        lnratN1MC18SUB3_centscan->Draw();
    }
    cratN1MC18SUB3_centscan->cd(1);
    TPaveText * txratN1MC18SUB3_centscan = new TPaveText(0.24, 0.77, 0.73, 0.97, "NDC");
    SetTPaveTxt(txratN1MC18SUB3_centscan, 16);
    txratN1MC18SUB3_centscan->AddText("PbPb #sqrt{s_{NN}}=5.02 TeV");
    txratN1MC18SUB3_centscan->AddText(Form("POI: %1.1f<|#eta|<%1.1f",eminETA[erangeNeg],emaxETA[erangeNeg]));
    txratN1MC18SUB3_centscan->Draw();
    cN1MC18SUB3_centscan->Print("results/Pt_Distributions/N1MC18/centscan_N1MCpm18SUB3_ratio.png","png");

}
