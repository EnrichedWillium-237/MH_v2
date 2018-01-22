# include "TCanvas.h"
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

static const int cbinsCENT = 13;
static const int cminCENT[] = {0,  5, 10, 15, 20, 25, 30, 35, 40, 50, 60,  0, 20,  60};
static const int cmaxCENT[] = {5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 20, 60, 100};
static const int ndeltaEta = 11;
static const double deltaEta[] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4};
static const int ndeltaEtaBins = 12;
static const double deltaEtaBins[] = {-0.2, 0.2, 0.6, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0, 3.4, 3.8, 4.2, 4.6};
static const int netabins = 12;
static const double etabins[] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0,  2.4};
static const double etaMid[] = {-2.2, -1.8, -1.4, -1.0, -0.6, -0.2,  0.2,  0.6,  1.0,  1.4,  1.8,  2.2};
static int NANALS = 24;
TString ANAL[] = {"N1MCm22SUB2", "N1MCm18SUB2", "N1MCm14SUB2", "N1MCm10SUB2", "N1MCm06SUB2", "N1MCm02SUB2",
                 "N1MCp02SUB2", "N1MCp06SUB2", "N1MCp10SUB2", "N1MCp14SUB2", "N1MCp18SUB2", "N1MCp22SUB2",
                 "N1MCm22SUB3", "N1MCm18SUB3", "N1MCm14SUB3", "N1MCm10SUB3", "N1MCm06SUB3", "N1MCm02SUB3",
                 "N1MCp02SUB3", "N1MCp06SUB3", "N1MCp10SUB3", "N1MCp14SUB3", "N1MCp18SUB3", "N1MCp22SUB3"
};

TH1D * gN1MCm22SUB2[cbinsCENT];
TH1D * gN1MCm18SUB2[cbinsCENT];
TH1D * gN1MCm14SUB2[cbinsCENT];
TH1D * gN1MCm10SUB2[cbinsCENT];
TH1D * gN1MCm06SUB2[cbinsCENT];
TH1D * gN1MCm02SUB2[cbinsCENT];
TH1D * gN1MCp02SUB2[cbinsCENT];
TH1D * gN1MCp06SUB2[cbinsCENT];
TH1D * gN1MCp10SUB2[cbinsCENT];
TH1D * gN1MCp14SUB2[cbinsCENT];
TH1D * gN1MCp18SUB2[cbinsCENT];
TH1D * gN1MCp22SUB2[cbinsCENT];

TH1D * gN1MCm22SUB3[cbinsCENT];
TH1D * gN1MCm18SUB3[cbinsCENT];
TH1D * gN1MCm14SUB3[cbinsCENT];
TH1D * gN1MCm10SUB3[cbinsCENT];
TH1D * gN1MCm06SUB3[cbinsCENT];
TH1D * gN1MCm02SUB3[cbinsCENT];
TH1D * gN1MCp02SUB3[cbinsCENT];
TH1D * gN1MCp06SUB3[cbinsCENT];
TH1D * gN1MCp10SUB3[cbinsCENT];
TH1D * gN1MCp14SUB3[cbinsCENT];
TH1D * gN1MCp18SUB3[cbinsCENT];
TH1D * gN1MCp22SUB3[cbinsCENT];

void EtaGap()
{

    TH1::SetDefaultSumw2();

    for (int cbin = 0; cbin<cbinsCENT; cbin++) {
        gN1MCm22SUB2[cbin] = new TH1D(Form("N1MCm22SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins, deltaEtaBins);
        gN1MCm18SUB2[cbin] = new TH1D(Form("N1MCm18SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins-1, deltaEtaBins);
        gN1MCm14SUB2[cbin] = new TH1D(Form("N1MCm14SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins-2, deltaEtaBins);
        gN1MCm10SUB2[cbin] = new TH1D(Form("N1MCm10SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins-3, deltaEtaBins);
        gN1MCm06SUB2[cbin] = new TH1D(Form("N1MCm06SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins-4, deltaEtaBins);
        gN1MCm02SUB2[cbin] = new TH1D(Form("N1MCm02SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins-5, deltaEtaBins);
        gN1MCp02SUB2[cbin] = new TH1D(Form("N1MCp02SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins-5, deltaEtaBins);
        gN1MCp06SUB2[cbin] = new TH1D(Form("N1MCp06SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins-4, deltaEtaBins);
        gN1MCp10SUB2[cbin] = new TH1D(Form("N1MCp10SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins-3, deltaEtaBins);
        gN1MCp14SUB2[cbin] = new TH1D(Form("N1MCp14SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins-2, deltaEtaBins);
        gN1MCp18SUB2[cbin] = new TH1D(Form("N1MCp18SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins-1, deltaEtaBins);
        gN1MCp22SUB2[cbin] = new TH1D(Form("N1MCp22SUB2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins, deltaEtaBins);

        gN1MCm22SUB3[cbin] = new TH1D(Form("N1MCm22SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins, deltaEtaBins);
        gN1MCm18SUB3[cbin] = new TH1D(Form("N1MCm18SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins-1, deltaEtaBins);
        gN1MCm14SUB3[cbin] = new TH1D(Form("N1MCm14SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins-2, deltaEtaBins);
        gN1MCm10SUB3[cbin] = new TH1D(Form("N1MCm10SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins-3, deltaEtaBins);
        gN1MCm06SUB3[cbin] = new TH1D(Form("N1MCm06SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins-4, deltaEtaBins);
        gN1MCm02SUB3[cbin] = new TH1D(Form("N1MCm02SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins-5, deltaEtaBins);
        gN1MCp02SUB3[cbin] = new TH1D(Form("N1MCp02SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins-5, deltaEtaBins);
        gN1MCp06SUB3[cbin] = new TH1D(Form("N1MCp06SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins-4, deltaEtaBins);
        gN1MCp10SUB3[cbin] = new TH1D(Form("N1MCp10SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins-3, deltaEtaBins);
        gN1MCp14SUB3[cbin] = new TH1D(Form("N1MCp14SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins-2, deltaEtaBins);
        gN1MCp18SUB3[cbin] = new TH1D(Form("N1MCp18SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins-1, deltaEtaBins);
        gN1MCp22SUB3[cbin] = new TH1D(Form("N1MCp22SUB3_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", ndeltaEtaBins, deltaEtaBins);
    }

    for (int i = 0; i<NANALS; i++) {

        for (int cbin = 0; cbin<cbinsCENT; cbin++) {
            TString tag = Form("figures_MH/%s/EtaDistributions/data/EtaInt_%s_%d_%d.dat",ANAL[i].Data(),ANAL[i].Data(),cminCENT[cbin],cmaxCENT[cbin]);
            // cout<<tag.Data()<<endl;
            ifstream fin(tag.Data());

            double etaval[20] = {0};
            double v1val[20] = {0};
            double v1vale[20] = {0};
            double v1Aval[20] = {0};
            double v1Avale[20] = {0};
            double v1Bval[20] = {0};
            double v1Bvale[20] = {0};
            double eta, v1, v1e, v1A, v1Ae, v1B, v1Be;
            int neta = 0;
            while (fin >> eta >> v1 >> v1e >> v1A >> v1Ae >> v1B >> v1Be) {
                // if (i == 0 && cbin == 0) cout<<"neta = "<<neta<<"\tdeltaEta = "<<deltaEta[neta]<<"\tv1A = "<<v1A<<"\tv1Ae = "<<v1Ae<<endl;
                etaval[neta] = eta;
                v1val[neta] = fabs(v1A);
                v1vale[neta] = v1Ae;
                neta++;
            }
            for (int ebin = 0; ebin<neta; ebin++) {
                if (i == 0) { // neg side
                    gN1MCm22SUB2[cbin]->SetBinContent(ebin+1, v1val[ebin+i]);
                    gN1MCm22SUB2[cbin]->SetBinError(ebin+1, v1vale[ebin+i]);
                } else if (i == 1) {
                    gN1MCm18SUB2[cbin]->SetBinContent(ebin+1, v1val[ebin+i]);
                    gN1MCm18SUB2[cbin]->SetBinError(ebin+1, v1vale[ebin+i]);
                } else if (i == 2) {
                    gN1MCm14SUB2[cbin]->SetBinContent(ebin+1, v1val[ebin+i]);
                    gN1MCm14SUB2[cbin]->SetBinError(ebin+1, v1vale[ebin+i]);
                } else if (i == 3) {
                    gN1MCm10SUB2[cbin]->SetBinContent(ebin+1, v1val[ebin+i]);
                    gN1MCm10SUB2[cbin]->SetBinError(ebin+1, v1vale[ebin+i]);
                } else if (i == 4) {
                    gN1MCm06SUB2[cbin]->SetBinContent(ebin+1, v1val[ebin+i]);
                    gN1MCm06SUB2[cbin]->SetBinError(ebin+1, v1vale[ebin+i]);
                } else if (i == 5) {
                    gN1MCm02SUB2[cbin]->SetBinContent(ebin+1, v1val[ebin+i]);
                    gN1MCm02SUB2[cbin]->SetBinError(ebin+1, v1vale[ebin+i]);
                } else if (i == 6) { // pos side
                    gN1MCp02SUB2[cbin]->SetBinContent(ebin+1, v1val[neta-ebin-6]);
                    gN1MCp02SUB2[cbin]->SetBinError(ebin+1, v1vale[neta-ebin-6]);
                } else if (i == 7) {
                    gN1MCp06SUB2[cbin]->SetBinContent(ebin+1, v1val[neta-ebin-5]);
                    gN1MCp06SUB2[cbin]->SetBinError(ebin+1, v1vale[neta-ebin-5]);
                } else if (i == 8) {
                    gN1MCp10SUB2[cbin]->SetBinContent(ebin+1, v1val[neta-ebin-4]);
                    gN1MCp10SUB2[cbin]->SetBinError(ebin+1, v1vale[neta-ebin-4]);
                } else if (i == 9) {
                    gN1MCp14SUB2[cbin]->SetBinContent(ebin+1, v1val[neta-ebin-3]);
                    gN1MCp14SUB2[cbin]->SetBinError(ebin+1, v1vale[neta-ebin-3]);
                } else if (i == 10) {
                    gN1MCp18SUB2[cbin]->SetBinContent(ebin+1, v1val[neta-ebin-2]);
                    gN1MCp18SUB2[cbin]->SetBinError(ebin+1, v1vale[neta-ebin-2]);
                } else if (i == 11) {
                    gN1MCp22SUB2[cbin]->SetBinContent(ebin+1, v1val[neta-ebin-1]);
                    gN1MCp22SUB2[cbin]->SetBinError(ebin+1, v1vale[neta-ebin-1]);

                } else if (i == 12) {
                    gN1MCm22SUB3[cbin]->SetBinContent(ebin+1, v1val[ebin+i-12]);
                    gN1MCm22SUB3[cbin]->SetBinError(ebin+1, v1vale[ebin+i-12]);
                } else if (i == 13) {
                    gN1MCm18SUB3[cbin]->SetBinContent(ebin+1, v1val[ebin+i-12]);
                    gN1MCm18SUB3[cbin]->SetBinError(ebin+1, v1vale[ebin+i-12]);
                } else if (i == 14) {
                    gN1MCm14SUB3[cbin]->SetBinContent(ebin+1, v1val[ebin+i-12]);
                    gN1MCm14SUB3[cbin]->SetBinError(ebin+1, v1vale[ebin+i-12]);
                } else if (i == 15) {
                    gN1MCm10SUB3[cbin]->SetBinContent(ebin+1, v1val[ebin+i-12]);
                    gN1MCm10SUB3[cbin]->SetBinError(ebin+1, v1vale[ebin+i-12]);
                } else if (i == 16) {
                    gN1MCm06SUB3[cbin]->SetBinContent(ebin+1, v1val[ebin+i-12]);
                    gN1MCm06SUB3[cbin]->SetBinError(ebin+1, v1vale[ebin+i-12]);
                } else if (i == 17) {
                    gN1MCm02SUB3[cbin]->SetBinContent(ebin+1, v1val[ebin+i-12]);
                    gN1MCm02SUB3[cbin]->SetBinError(ebin+1, v1vale[ebin+i-12]);
                } else if (i == 18) { // pos side
                    gN1MCp02SUB3[cbin]->SetBinContent(ebin+1, v1val[neta-ebin-6]);
                    gN1MCp02SUB3[cbin]->SetBinError(ebin+1, v1vale[neta-ebin-6]);
                } else if (i == 19) {
                    gN1MCp06SUB3[cbin]->SetBinContent(ebin+1, v1val[neta-ebin-5]);
                    gN1MCp06SUB3[cbin]->SetBinError(ebin+1, v1vale[neta-ebin-5]);
                } else if (i == 20) {
                    gN1MCp10SUB3[cbin]->SetBinContent(ebin+1, v1val[neta-ebin-4]);
                    gN1MCp10SUB3[cbin]->SetBinError(ebin+1, v1vale[neta-ebin-4]);
                } else if (i == 21) {
                    gN1MCp14SUB3[cbin]->SetBinContent(ebin+1, v1val[neta-ebin-3]);
                    gN1MCp14SUB3[cbin]->SetBinError(ebin+1, v1vale[neta-ebin-3]);
                } else if (i == 22) {
                    gN1MCp18SUB3[cbin]->SetBinContent(ebin+1, v1val[neta-ebin-2]);
                    gN1MCp18SUB3[cbin]->SetBinError(ebin+1, v1vale[neta-ebin-2]);
                } else if (i == 23) {
                    gN1MCp22SUB3[cbin]->SetBinContent(ebin+1, v1val[neta-ebin-1]);
                    gN1MCp22SUB3[cbin]->SetBinError(ebin+1, v1vale[neta-ebin-1]);
                    // if (cbin == 0) cout<<"i: "<<i<<"\tebin: "<<ebin<<"\tdeltaEta: "<<deltaEta[ebin]<<"\tneta-ebin-1: "<<neta-ebin-1<<"\tv1val: "<<v1val[neta-ebin-1]<<endl;
                } else {continue;}
            }

        }
    }

    // draw options
    for (int cbin = 0; cbin<cbinsCENT; cbin++) {
        gN1MCm22SUB2[cbin]->SetMarkerColor(kRed);
        gN1MCm18SUB2[cbin]->SetMarkerColor(kBlue);
        gN1MCm14SUB2[cbin]->SetMarkerColor(kGreen+2);
        gN1MCm10SUB2[cbin]->SetMarkerColor(kMagenta);
        gN1MCm06SUB2[cbin]->SetMarkerColor(kOrange+7);
        gN1MCm02SUB2[cbin]->SetMarkerColor(kBlack);
        gN1MCp02SUB2[cbin]->SetMarkerColor(kBlack);
        gN1MCp06SUB2[cbin]->SetMarkerColor(kOrange+7);
        gN1MCp10SUB2[cbin]->SetMarkerColor(kMagenta);
        gN1MCp14SUB2[cbin]->SetMarkerColor(kGreen+2);
        gN1MCp18SUB2[cbin]->SetMarkerColor(kBlue);
        gN1MCp22SUB2[cbin]->SetMarkerColor(kRed);

        gN1MCm22SUB2[cbin]->SetLineColor(kRed);
        gN1MCm18SUB2[cbin]->SetLineColor(kBlue);
        gN1MCm14SUB2[cbin]->SetLineColor(kGreen+2);
        gN1MCm10SUB2[cbin]->SetLineColor(kMagenta);
        gN1MCm06SUB2[cbin]->SetLineColor(kOrange+7);
        gN1MCm02SUB2[cbin]->SetLineColor(kBlack);
        gN1MCp02SUB2[cbin]->SetLineColor(kBlack);
        gN1MCp06SUB2[cbin]->SetLineColor(kOrange+7);
        gN1MCp10SUB2[cbin]->SetLineColor(kMagenta);
        gN1MCp14SUB2[cbin]->SetLineColor(kGreen+2);
        gN1MCp18SUB2[cbin]->SetLineColor(kBlue);
        gN1MCp22SUB2[cbin]->SetLineColor(kRed);

        gN1MCm22SUB2[cbin]->SetMarkerStyle(20);
        gN1MCm18SUB2[cbin]->SetMarkerStyle(21);
        gN1MCm14SUB2[cbin]->SetMarkerStyle(34);
        gN1MCm10SUB2[cbin]->SetMarkerStyle(33);
        gN1MCm06SUB2[cbin]->SetMarkerStyle(21);
        gN1MCm02SUB2[cbin]->SetMarkerStyle(34);
        gN1MCp02SUB2[cbin]->SetMarkerStyle(28);
        gN1MCp06SUB2[cbin]->SetMarkerStyle(25);
        gN1MCp10SUB2[cbin]->SetMarkerStyle(27);
        gN1MCp14SUB2[cbin]->SetMarkerStyle(28);
        gN1MCp18SUB2[cbin]->SetMarkerStyle(25);
        gN1MCp22SUB2[cbin]->SetMarkerStyle(24);

        gN1MCm22SUB2[cbin]->SetMarkerSize(1.3);
        gN1MCm18SUB2[cbin]->SetMarkerSize(1.2);
        gN1MCm14SUB2[cbin]->SetMarkerSize(1.7);
        gN1MCm10SUB2[cbin]->SetMarkerSize(1.8);
        gN1MCm06SUB2[cbin]->SetMarkerSize(1.2);
        gN1MCm02SUB2[cbin]->SetMarkerSize(1.7);
        gN1MCp02SUB2[cbin]->SetMarkerSize(1.7);
        gN1MCp06SUB2[cbin]->SetMarkerSize(1.2);
        gN1MCp10SUB2[cbin]->SetMarkerSize(1.8);
        gN1MCp14SUB2[cbin]->SetMarkerSize(1.7);
        gN1MCp18SUB2[cbin]->SetMarkerSize(1.2);
        gN1MCp22SUB2[cbin]->SetMarkerSize(1.3);
    }


    if (!fopen("figures_MH/EtaGapStudy/","r")) system("mkdir figures_MH/EtaGapStudy");

    int setcent = 5;

    TCanvas * cMCmScan = new TCanvas("cMCmScan","cMCmScan",650,600);
    TPad * padMCmScan = (TPad *) cMCmScan->cd();
    padMCmScan->SetGrid();
    TH1D * hMCmScan = new TH1D("hMCmScan", "", 100, 0, 4.8);
    hMCmScan->SetTitle("");
    hMCmScan->SetStats(0);
    hMCmScan->GetYaxis()->SetRangeUser(0, 0.04);
    hMCmScan->SetXTitle("#Delta#eta");
    hMCmScan->SetYTitle("|v_{1}^{even}|");
    hMCmScan->Draw();
    gN1MCm22SUB3[setcent-1]->Draw("same");
    gN1MCm18SUB3[setcent-1]->Draw("same");
    gN1MCm14SUB3[setcent-1]->Draw("same");
    gN1MCm10SUB3[setcent-1]->Draw("same");
    gN1MCm06SUB3[setcent-1]->Draw("same");
    gN1MCm02SUB3[setcent-1]->Draw("same");
    TLegend * legMCmScan = new TLegend(0.66, 0.67, 0.93, 0.93);
    SetLegend(legMCmScan,18);
    legMCmScan->AddEntry(gN1MCm22SUB3[setcent-1],"N1MCm22SUB3","p");
    legMCmScan->AddEntry(gN1MCm18SUB3[setcent-1],"N1MCm18SUB3","p");
    legMCmScan->AddEntry(gN1MCm14SUB3[setcent-1],"N1MCm14SUB3","p");
    legMCmScan->AddEntry(gN1MCm10SUB3[setcent-1],"N1MCm10SUB3","p");
    legMCmScan->AddEntry(gN1MCm06SUB3[setcent-1],"N1MCm06SUB3","p");
    legMCmScan->AddEntry(gN1MCm02SUB3[setcent-1],"N1MCm02SUB3","p");
    legMCmScan->Draw();
    TPaveText * txMCmScan = new TPaveText(0.2, 0.81, 0.47, 0.92, "NDC");
    SetTPaveTxt(txMCmScan, 20);
    txMCmScan->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txMCmScan->AddText(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    txMCmScan->Draw();
    cMCmScan->Print(Form("figures_MH/EtaGapStudy/MCm_deltaEta_scan_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");


    TCanvas * cMCpScan = new TCanvas("cMCpScan","cMCpScan",650,600);
    TPad * padMCpScan = (TPad *) cMCpScan->cd();
    padMCpScan->SetGrid();
    TH1D * hMCpScan = new TH1D("hMCpScan", "", 100, 0, 4.8);
    hMCpScan->SetTitle("");
    hMCpScan->SetStats(0);
    hMCpScan->GetYaxis()->SetRangeUser(0, 0.04);
    hMCpScan->SetXTitle("#Delta#eta");
    hMCpScan->SetYTitle("|v_{1}^{even}|");
    hMCpScan->Draw();
    gN1MCp22SUB3[setcent-1]->Draw("same");
    gN1MCp18SUB3[setcent-1]->Draw("same");
    gN1MCp14SUB3[setcent-1]->Draw("same");
    gN1MCp10SUB3[setcent-1]->Draw("same");
    gN1MCp06SUB3[setcent-1]->Draw("same");
    gN1MCp02SUB3[setcent-1]->Draw("same");
    TLegend * legMCpScan = new TLegend(0.66, 0.67, 0.93, 0.93);
    SetLegend(legMCpScan,18);
    legMCpScan->AddEntry(gN1MCp22SUB3[setcent-1],"N1MCp22SUB3","p");
    legMCpScan->AddEntry(gN1MCp18SUB3[setcent-1],"N1MCp18SUB3","p");
    legMCpScan->AddEntry(gN1MCp14SUB3[setcent-1],"N1MCp14SUB3","p");
    legMCpScan->AddEntry(gN1MCp10SUB3[setcent-1],"N1MCp10SUB3","p");
    legMCpScan->AddEntry(gN1MCp06SUB3[setcent-1],"N1MCp06SUB3","p");
    legMCpScan->AddEntry(gN1MCp02SUB3[setcent-1],"N1MCp02SUB3","p");
    legMCpScan->Draw();
    TPaveText * txMCpScan = new TPaveText(0.2, 0.81, 0.47, 0.92, "NDC");
    SetTPaveTxt(txMCpScan, 20);
    txMCpScan->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txMCpScan->AddText(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    txMCpScan->Draw();
    cMCpScan->Print(Form("figures_MH/EtaGapStudy/MCp_deltaEta_scan_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");


    TCanvas * cMCmCent22 = new TCanvas("cMCmCent22","cMCmCent22",650,600);
    TPad * padMCmCent22 = (TPad *) cMCmCent22->cd();
    padMCmCent22->SetGrid();
    TH1D * hMCmCent22 = new TH1D("hMCmCent22", "", 100, 0, 4.8);
    hMCmCent22->SetTitle("");
    hMCmCent22->SetStats(0);
    hMCmCent22->GetYaxis()->SetRangeUser(0, 0.04);
    hMCmCent22->SetXTitle("#Delta#eta");
    hMCmCent22->SetYTitle("|v_{1}^{even}|");
    hMCmCent22->Draw();
    TH1D * gN1MCm22SUB3_tmp_0 = (TH1D *) gN1MCm22SUB3[0]->Clone("gN1MCm22SUB3_tmp_0");
    TH1D * gN1MCm22SUB3_tmp_7 = (TH1D *) gN1MCm22SUB3[7]->Clone("gN1MCm22SUB3_tmp_7");
    TH1D * gN1MCm22SUB3_tmp_8 = (TH1D *) gN1MCm22SUB3[8]->Clone("gN1MCm22SUB3_tmp_8");
    TH1D * gN1MCm22SUB3_tmp_9 = (TH1D *) gN1MCm22SUB3[9]->Clone("gN1MCm22SUB3_tmp_9");
    gN1MCm22SUB3_tmp_0->SetMarkerColor(kRed);
    gN1MCm22SUB3_tmp_7->SetMarkerColor(kBlue);
    gN1MCm22SUB3_tmp_8->SetMarkerColor(kGreen+2);
    gN1MCm22SUB3_tmp_9->SetMarkerColor(kMagenta);
    gN1MCm22SUB3_tmp_0->SetLineColor(kRed);
    gN1MCm22SUB3_tmp_7->SetLineColor(kBlue);
    gN1MCm22SUB3_tmp_8->SetLineColor(kGreen+2);
    gN1MCm22SUB3_tmp_9->SetLineColor(kMagenta);
    gN1MCm22SUB3_tmp_0->SetMarkerStyle(21);
    gN1MCm22SUB3_tmp_7->SetMarkerStyle(27);
    gN1MCm22SUB3_tmp_8->SetMarkerStyle(34);
    gN1MCm22SUB3_tmp_9->SetMarkerStyle(24);
    gN1MCm22SUB3_tmp_0->SetMarkerSize(1.2);
    gN1MCm22SUB3_tmp_7->SetMarkerSize(1.8);
    gN1MCm22SUB3_tmp_8->SetMarkerSize(1.7);
    gN1MCm22SUB3_tmp_9->SetMarkerSize(1.3);
    gN1MCm22SUB3_tmp_0->Draw("same");
    gN1MCm22SUB3_tmp_7->Draw("same");
    gN1MCm22SUB3_tmp_8->Draw("same");
    gN1MCm22SUB3_tmp_9->Draw("same");
    TLegend * legMCmCent22 = new TLegend(0.74, 0.74, 0.93, 0.93);
    SetLegend(legMCmCent22, 18);
    legMCmCent22->AddEntry(gN1MCm22SUB3_tmp_0,Form(" %d-%d%%",cminCENT[0],cmaxCENT[0]),"p");
    legMCmCent22->AddEntry(gN1MCm22SUB3_tmp_7,Form(" %d-%d%%",cminCENT[7],cmaxCENT[7]),"p");
    legMCmCent22->AddEntry(gN1MCm22SUB3_tmp_8,Form(" %d-%d%%",cminCENT[8],cmaxCENT[8]),"p");
    legMCmCent22->AddEntry(gN1MCm22SUB3_tmp_9,Form(" %d-%d%%",cminCENT[9],cmaxCENT[9]),"p");
    legMCmCent22->Draw();
    TPaveText * txMCmCent22 = new TPaveText(0.2, 0.81, 0.65, 0.92, "NDC");
    SetTPaveTxt(txMCmCent22, 20);
    txMCmCent22->AddText("N1MCm22SUB3 (-2.4 < #eta < -2.0)");
    txMCmCent22->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txMCmCent22->Draw();
    cMCmCent22->Print("figures_MH/EtaGapStudy/MCm22_Cent.png","png");


    TCanvas * cMCpCent22 = new TCanvas("cMCpCent22","cMCpCent22",650,600);
    TPad * padMCpCent22 = (TPad *) cMCpCent22->cd();
    padMCpCent22->SetGrid();
    TH1D * hMCpCent22 = new TH1D("hMCpCent22", "", 100, 0, 4.8);
    hMCpCent22->SetTitle("");
    hMCpCent22->SetStats(0);
    hMCpCent22->GetYaxis()->SetRangeUser(0, 0.04);
    hMCpCent22->SetXTitle("#Delta#eta");
    hMCpCent22->SetYTitle("|v_{1}^{even}|");
    hMCpCent22->Draw();
    TH1D * gN1MCp22SUB3_tmp_0 = (TH1D *) gN1MCp22SUB3[0]->Clone("gN1MCp22SUB3_tmp_0");
    TH1D * gN1MCp22SUB3_tmp_7 = (TH1D *) gN1MCp22SUB3[4]->Clone("gN1MCp22SUB3_tmp_7");
    TH1D * gN1MCp22SUB3_tmp_8 = (TH1D *) gN1MCp22SUB3[7]->Clone("gN1MCp22SUB3_tmp_8");
    TH1D * gN1MCp22SUB3_tmp_9 = (TH1D *) gN1MCp22SUB3[9]->Clone("gN1MCp22SUB3_tmp_9");
    gN1MCp22SUB3_tmp_0->SetMarkerColor(kRed);
    gN1MCp22SUB3_tmp_7->SetMarkerColor(kBlue);
    gN1MCp22SUB3_tmp_8->SetMarkerColor(kGreen+2);
    gN1MCp22SUB3_tmp_9->SetMarkerColor(kMagenta);
    gN1MCp22SUB3_tmp_0->SetLineColor(kRed);
    gN1MCp22SUB3_tmp_7->SetLineColor(kBlue);
    gN1MCp22SUB3_tmp_8->SetLineColor(kGreen+2);
    gN1MCp22SUB3_tmp_9->SetLineColor(kMagenta);
    gN1MCp22SUB3_tmp_0->SetMarkerStyle(21);
    gN1MCp22SUB3_tmp_7->SetMarkerStyle(27);
    gN1MCp22SUB3_tmp_8->SetMarkerStyle(34);
    gN1MCp22SUB3_tmp_9->SetMarkerStyle(24);
    gN1MCp22SUB3_tmp_0->SetMarkerSize(1.2);
    gN1MCp22SUB3_tmp_7->SetMarkerSize(1.8);
    gN1MCp22SUB3_tmp_8->SetMarkerSize(1.7);
    gN1MCp22SUB3_tmp_9->SetMarkerSize(1.3);
    gN1MCp22SUB3_tmp_0->Draw("same");
    gN1MCp22SUB3_tmp_7->Draw("same");
    gN1MCp22SUB3_tmp_8->Draw("same");
    gN1MCp22SUB3_tmp_9->Draw("same");
    TLegend * legMCpCent22 = new TLegend(0.74, 0.74, 0.93, 0.93);
    SetLegend(legMCpCent22, 18);
    legMCpCent22->AddEntry(gN1MCp22SUB3_tmp_0,Form(" %d-%d%%",cminCENT[0],cmaxCENT[0]),"p");
    legMCpCent22->AddEntry(gN1MCp22SUB3_tmp_7,Form(" %d-%d%%",cminCENT[7],cmaxCENT[7]),"p");
    legMCpCent22->AddEntry(gN1MCp22SUB3_tmp_8,Form(" %d-%d%%",cminCENT[8],cmaxCENT[8]),"p");
    legMCpCent22->AddEntry(gN1MCp22SUB3_tmp_9,Form(" %d-%d%%",cminCENT[9],cmaxCENT[9]),"p");
    legMCpCent22->Draw();
    TPaveText * txMCpCent22 = new TPaveText(0.2, 0.81, 0.65, 0.92, "NDC");
    SetTPaveTxt(txMCpCent22, 20);
    txMCpCent22->AddText("N1MCp22SUB3 (2.0 < #eta < 2.4)");
    txMCpCent22->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txMCpCent22->Draw();
    cMCpCent22->Print("figures_MH/EtaGapStudy/MCp22_Cent.png","png");


    TCanvas * cMCmCent18 = new TCanvas("cMCmCent18","cMCmCent18",650,600);
    TPad * padMCmCent18 = (TPad *) cMCmCent18->cd();
    padMCmCent18->SetGrid();
    TH1D * hMCmCent18 = new TH1D("hMCmCent18", "", 100, 0, 4.8);
    hMCmCent18->SetTitle("");
    hMCmCent18->SetStats(0);
    hMCmCent18->GetYaxis()->SetRangeUser(0, 0.04);
    hMCmCent18->SetXTitle("#Delta#eta");
    hMCmCent18->SetYTitle("|v_{1}^{even}|");
    hMCmCent18->Draw();
    TH1D * gN1MCm18SUB3_tmp_0 = (TH1D *) gN1MCm18SUB3[0]->Clone("gN1MCm18SUB3_tmp_0");
    TH1D * gN1MCm18SUB3_tmp_7 = (TH1D *) gN1MCm18SUB3[7]->Clone("gN1MCm18SUB3_tmp_7");
    TH1D * gN1MCm18SUB3_tmp_8 = (TH1D *) gN1MCm18SUB3[8]->Clone("gN1MCm18SUB3_tmp_8");
    TH1D * gN1MCm18SUB3_tmp_9 = (TH1D *) gN1MCm18SUB3[9]->Clone("gN1MCm18SUB3_tmp_9");
    gN1MCm18SUB3_tmp_0->SetMarkerColor(kRed);
    gN1MCm18SUB3_tmp_7->SetMarkerColor(kBlue);
    gN1MCm18SUB3_tmp_8->SetMarkerColor(kGreen+2);
    gN1MCm18SUB3_tmp_9->SetMarkerColor(kMagenta);
    gN1MCm18SUB3_tmp_0->SetLineColor(kRed);
    gN1MCm18SUB3_tmp_7->SetLineColor(kBlue);
    gN1MCm18SUB3_tmp_8->SetLineColor(kGreen+2);
    gN1MCm18SUB3_tmp_9->SetLineColor(kMagenta);
    gN1MCm18SUB3_tmp_0->SetMarkerStyle(21);
    gN1MCm18SUB3_tmp_7->SetMarkerStyle(27);
    gN1MCm18SUB3_tmp_8->SetMarkerStyle(34);
    gN1MCm18SUB3_tmp_9->SetMarkerStyle(24);
    gN1MCm18SUB3_tmp_0->SetMarkerSize(1.2);
    gN1MCm18SUB3_tmp_7->SetMarkerSize(1.8);
    gN1MCm18SUB3_tmp_8->SetMarkerSize(1.7);
    gN1MCm18SUB3_tmp_9->SetMarkerSize(1.3);
    gN1MCm18SUB3_tmp_0->Draw("same");
    gN1MCm18SUB3_tmp_7->Draw("same");
    gN1MCm18SUB3_tmp_8->Draw("same");
    gN1MCm18SUB3_tmp_9->Draw("same");
    TLegend * legMCmCent18 = new TLegend(0.74, 0.74, 0.93, 0.93);
    SetLegend(legMCmCent18, 18);
    legMCmCent18->AddEntry(gN1MCm18SUB3_tmp_0,Form(" %d-%d%%",cminCENT[0],cmaxCENT[0]),"p");
    legMCmCent18->AddEntry(gN1MCm18SUB3_tmp_7,Form(" %d-%d%%",cminCENT[7],cmaxCENT[7]),"p");
    legMCmCent18->AddEntry(gN1MCm18SUB3_tmp_8,Form(" %d-%d%%",cminCENT[8],cmaxCENT[8]),"p");
    legMCmCent18->AddEntry(gN1MCm18SUB3_tmp_9,Form(" %d-%d%%",cminCENT[9],cmaxCENT[9]),"p");
    legMCmCent18->Draw();
    TPaveText * txMCmCent18 = new TPaveText(0.2, 0.81, 0.65, 0.92, "NDC");
    SetTPaveTxt(txMCmCent18, 20);
    txMCmCent18->AddText("N1MCm18SUB3 (-2.0 < #eta < -1.6)");
    txMCmCent18->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txMCmCent18->Draw();
    cMCmCent18->Print("figures_MH/EtaGapStudy/MCm18_Cent.png","png");


    TCanvas * cMCpCent18 = new TCanvas("cMCpCent18","cMCpCent18",650,600);
    TPad * padMCpCent18 = (TPad *) cMCpCent18->cd();
    padMCpCent18->SetGrid();
    TH1D * hMCpCent18 = new TH1D("hMCpCent18", "", 100, 0, 4.8);
    hMCpCent18->SetTitle("");
    hMCpCent18->SetStats(0);
    hMCpCent18->GetYaxis()->SetRangeUser(0, 0.04);
    hMCpCent18->SetXTitle("#Delta#eta");
    hMCpCent18->SetYTitle("|v_{1}^{even}|");
    hMCpCent18->Draw();
    TH1D * gN1MCp18SUB3_tmp_0 = (TH1D *) gN1MCp18SUB3[0]->Clone("gN1MCp18SUB3_tmp_0");
    TH1D * gN1MCp18SUB3_tmp_7 = (TH1D *) gN1MCp18SUB3[7]->Clone("gN1MCp18SUB3_tmp_7");
    TH1D * gN1MCp18SUB3_tmp_8 = (TH1D *) gN1MCp18SUB3[8]->Clone("gN1MCp18SUB3_tmp_8");
    TH1D * gN1MCp18SUB3_tmp_9 = (TH1D *) gN1MCp18SUB3[9]->Clone("gN1MCp18SUB3_tmp_9");
    gN1MCp18SUB3_tmp_0->SetMarkerColor(kBlue);
    gN1MCp18SUB3_tmp_7->SetMarkerColor(kRed);
    gN1MCp18SUB3_tmp_8->SetMarkerColor(kGreen+2);
    gN1MCp18SUB3_tmp_9->SetMarkerColor(kMagenta);
    gN1MCp18SUB3_tmp_0->SetLineColor(kBlue);
    gN1MCp18SUB3_tmp_7->SetLineColor(kRed);
    gN1MCp18SUB3_tmp_8->SetLineColor(kGreen+2);
    gN1MCp18SUB3_tmp_9->SetLineColor(kMagenta);
    gN1MCp18SUB3_tmp_0->SetMarkerStyle(21);
    gN1MCp18SUB3_tmp_7->SetMarkerStyle(27);
    gN1MCp18SUB3_tmp_8->SetMarkerStyle(34);
    gN1MCp18SUB3_tmp_9->SetMarkerStyle(24);
    gN1MCp18SUB3_tmp_0->SetMarkerSize(1.2);
    gN1MCp18SUB3_tmp_7->SetMarkerSize(1.8);
    gN1MCp18SUB3_tmp_8->SetMarkerSize(1.7);
    gN1MCp18SUB3_tmp_9->SetMarkerSize(1.3);
    gN1MCp18SUB3_tmp_0->Draw("same");
    gN1MCp18SUB3_tmp_7->Draw("same");
    gN1MCp18SUB3_tmp_8->Draw("same");
    gN1MCp18SUB3_tmp_9->Draw("same");
    TLegend * legMCpCent18 = new TLegend(0.74, 0.74, 0.93, 0.93);
    SetLegend(legMCpCent18, 18);
    legMCpCent18->AddEntry(gN1MCp18SUB3_tmp_0,Form(" %d-%d%%",cminCENT[0],cmaxCENT[0]),"p");
    legMCpCent18->AddEntry(gN1MCp18SUB3_tmp_7,Form(" %d-%d%%",cminCENT[7],cmaxCENT[7]),"p");
    legMCpCent18->AddEntry(gN1MCp18SUB3_tmp_8,Form(" %d-%d%%",cminCENT[8],cmaxCENT[8]),"p");
    legMCpCent18->AddEntry(gN1MCp18SUB3_tmp_9,Form(" %d-%d%%",cminCENT[9],cmaxCENT[9]),"p");
    legMCpCent18->Draw();
    TPaveText * txMCpCent18 = new TPaveText(0.2, 0.81, 0.65, 0.92, "NDC");
    SetTPaveTxt(txMCpCent18, 20);
    txMCpCent18->AddText("N1MCp18SUB3 (1.6 < #eta < 2.0)");
    txMCpCent18->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txMCpCent18->Draw();
    cMCpCent18->Print("figures_MH/EtaGapStudy/MCp18_Cent.png","png");


    setcent = 5;
    TCanvas * cComparepmMC22 = new TCanvas("cComparepmMC22","cComparepmMC22",650,600);
    TPad * padComparepmMC22_1 = (TPad *) cComparepmMC22->cd();
    padComparepmMC22_1->SetGrid();
    TH1D * hComparepmMC22 = new TH1D("hComparepmMC22", "", 100, 0, 4.8);
    hComparepmMC22->SetTitle("");
    hComparepmMC22->SetStats(0);
    hComparepmMC22->GetYaxis()->SetRangeUser(0, 0.04);
    hComparepmMC22->SetXTitle("#Delta#eta");
    hComparepmMC22->SetYTitle("|v_{1}^{even}|");
    hComparepmMC22->Draw();
    gN1MCm22SUB3[setcent]->SetMarkerColor(kRed);
    gN1MCm22SUB3[setcent]->SetLineColor(kRed);
    gN1MCm22SUB3[setcent]->SetMarkerStyle(20);
    gN1MCm22SUB3[setcent]->SetMarkerSize(1.3);
    gN1MCm22SUB3[setcent]->Draw("same");
    gN1MCp22SUB3[setcent]->SetMarkerColor(kRed);
    gN1MCp22SUB3[setcent]->SetLineColor(kRed);
    gN1MCp22SUB3[setcent]->SetMarkerStyle(24);
    gN1MCp22SUB3[setcent]->SetMarkerSize(1.3);
    gN1MCp22SUB3[setcent]->Draw("same");
    TLegend * legComparepmMC22 = new TLegend(0.48, 0.82, 0.78, 0.93);
    SetLegend(legComparepmMC22, 18);
    legComparepmMC22->AddEntry(gN1MCm22SUB3[setcent],"N1MCm22SUB3 (-2.4 < #eta < -2.0)","p");
    legComparepmMC22->AddEntry(gN1MCp22SUB3[setcent],"N1MCp22SUB3 (2.0 < #eta < 2.4)","p");
    legComparepmMC22->Draw();
    TPaveText * txComparepmMC22 = new TPaveText(0.2, 0.81, 0.49, 0.92, "NDC");
    SetTPaveTxt(txComparepmMC22, 20);
    txComparepmMC22->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txComparepmMC22->AddText(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    txComparepmMC22->Draw();
    cComparepmMC22->Print(Form("figures_MH/EtaGapStudy/Compare_pm_MC22_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");


    TCanvas * cComparepmMC18 = new TCanvas("cComparepmMC18","cComparepmMC18",650,600);
    TPad * padComparepmMC18_1 = (TPad *) cComparepmMC18->cd();
    padComparepmMC18_1->SetGrid();
    TH1D * hComparepmMC18 = new TH1D("hComparepmMC18", "", 100, 0, 4.8);
    hComparepmMC18->SetTitle("");
    hComparepmMC18->SetStats(0);
    hComparepmMC18->GetYaxis()->SetRangeUser(0, 0.04);
    hComparepmMC18->SetXTitle("#Delta#eta");
    hComparepmMC18->SetYTitle("|v_{1}^{even}|");
    hComparepmMC18->Draw();
    gN1MCm18SUB3[setcent]->SetMarkerColor(kBlue);
    gN1MCm18SUB3[setcent]->SetLineColor(kBlue);
    gN1MCm18SUB3[setcent]->SetMarkerStyle(21);
    gN1MCm18SUB3[setcent]->SetMarkerSize(1.2);
    gN1MCm18SUB3[setcent]->Draw("same");
    gN1MCp18SUB3[setcent]->SetMarkerColor(kBlue);
    gN1MCp18SUB3[setcent]->SetLineColor(kBlue);
    gN1MCp18SUB3[setcent]->SetMarkerStyle(25);
    gN1MCp18SUB3[setcent]->SetMarkerSize(1.2);
    gN1MCp18SUB3[setcent]->Draw("same");
    TLegend * legComparepmMC18 = new TLegend(0.48, 0.82, 0.78, 0.93);
    SetLegend(legComparepmMC18, 18);
    legComparepmMC18->AddEntry(gN1MCm18SUB3[setcent],"N1MCm18SUB3 (-2.0 < #eta < -1.6)","p");
    legComparepmMC18->AddEntry(gN1MCp18SUB3[setcent],"N1MCp18SUB3 (1.6 < #eta < 2.0)","p");
    legComparepmMC18->Draw();
    TPaveText * txComparepmMC18 = new TPaveText(0.2, 0.81, 0.49, 0.92, "NDC");
    SetTPaveTxt(txComparepmMC18, 20);
    txComparepmMC18->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txComparepmMC18->AddText(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    txComparepmMC18->Draw();
    cComparepmMC18->Print(Form("figures_MH/EtaGapStudy/Compare_pm_MC18_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");


    TCanvas * cComparepmMC14 = new TCanvas("cComparepmMC14","cComparepmMC14",650,600);
    TPad * padComparepmMC14_1 = (TPad *) cComparepmMC14->cd();
    padComparepmMC14_1->SetGrid();
    TH1D * hComparepmMC14 = new TH1D("hComparepmMC14", "", 100, 0, 4.8);
    hComparepmMC14->SetTitle("");
    hComparepmMC14->SetStats(0);
    hComparepmMC14->GetYaxis()->SetRangeUser(0, 0.04);
    hComparepmMC14->SetXTitle("#Delta#eta");
    hComparepmMC14->SetYTitle("|v_{1}^{even}|");
    hComparepmMC14->Draw();
    gN1MCm14SUB3[setcent]->SetMarkerColor(kGreen+2);
    gN1MCm14SUB3[setcent]->SetLineColor(kGreen+2);
    gN1MCm14SUB3[setcent]->SetMarkerStyle(34);
    gN1MCm14SUB3[setcent]->SetMarkerSize(1.7);
    gN1MCm14SUB3[setcent]->Draw("same");
    gN1MCp14SUB3[setcent]->SetMarkerColor(kGreen+2);
    gN1MCp14SUB3[setcent]->SetLineColor(kGreen+2);
    gN1MCp14SUB3[setcent]->SetMarkerStyle(28);
    gN1MCp14SUB3[setcent]->SetMarkerSize(1.7);
    gN1MCp14SUB3[setcent]->Draw("same");
    TLegend * legComparepmMC14 = new TLegend(0.48, 0.82, 0.78, 0.93);
    SetLegend(legComparepmMC14, 18);
    legComparepmMC14->AddEntry(gN1MCm14SUB3[setcent],"N1MCm14SUB3 (-1.6 < #eta < -1.2)","p");
    legComparepmMC14->AddEntry(gN1MCp14SUB3[setcent],"N1MCp14SUB3 (1.2 < #eta < 1.6)","p");
    legComparepmMC14->Draw();
    TPaveText * txComparepmMC14 = new TPaveText(0.2, 0.81, 0.49, 0.92, "NDC");
    SetTPaveTxt(txComparepmMC14, 20);
    txComparepmMC14->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txComparepmMC14->AddText(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    txComparepmMC14->Draw();
    cComparepmMC14->Print(Form("figures_MH/EtaGapStudy/Compare_pm_MC14_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");


    TCanvas * cComparepmMC10 = new TCanvas("cComparepmMC10","cComparepmMC10",650,600);
    TPad * padComparepmMC10_1 = (TPad *) cComparepmMC10->cd();
    padComparepmMC10_1->SetGrid();
    TH1D * hComparepmMC10 = new TH1D("hComparepmMC10", "", 100, 0, 4.8);
    hComparepmMC10->SetTitle("");
    hComparepmMC10->SetStats(0);
    hComparepmMC10->GetYaxis()->SetRangeUser(0, 0.04);
    hComparepmMC10->SetXTitle("#Delta#eta");
    hComparepmMC10->SetYTitle("|v_{1}^{even}|");
    hComparepmMC10->Draw();
    gN1MCm10SUB3[setcent]->Draw("same");
    gN1MCp10SUB3[setcent]->Draw("same");
    TLegend * legComparepmMC10 = new TLegend(0.48, 0.82, 0.78, 0.93);
    SetLegend(legComparepmMC10, 18);
    legComparepmMC10->AddEntry(gN1MCm10SUB3[setcent],"N1MCm10SUB3 (-1.2 < #eta < -0.8)","p");
    legComparepmMC10->AddEntry(gN1MCp10SUB3[setcent],"N1MCp10SUB3 (0.8 < #eta < 1.2)","p");
    legComparepmMC10->Draw();
    TPaveText * txComparepmMC10 = new TPaveText(0.2, 0.81, 0.49, 0.92, "NDC");
    SetTPaveTxt(txComparepmMC10, 20);
    txComparepmMC10->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txComparepmMC10->AddText(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    txComparepmMC10->Draw();
    cComparepmMC10->Print(Form("figures_MH/EtaGapStudy/Compare_pm_MC10_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");


    TCanvas * cComparepmMC06 = new TCanvas("cComparepmMC06","cComparepmMC06",650,600);
    TPad * padComparepmMC06_1 = (TPad *) cComparepmMC06->cd();
    padComparepmMC06_1->SetGrid();
    TH1D * hComparepmMC06 = new TH1D("hComparepmMC06", "", 100, 0, 4.8);
    hComparepmMC06->SetTitle("");
    hComparepmMC06->SetStats(0);
    hComparepmMC06->GetYaxis()->SetRangeUser(0, 0.04);
    hComparepmMC06->SetXTitle("#Delta#eta");
    hComparepmMC06->SetYTitle("|v_{1}^{even}|");
    hComparepmMC06->Draw();
    gN1MCm06SUB3[setcent]->Draw("same");
    gN1MCp06SUB3[setcent]->Draw("same");
    TLegend * legComparepmMC06 = new TLegend(0.48, 0.82, 0.78, 0.93);
    SetLegend(legComparepmMC06, 18);
    legComparepmMC06->AddEntry(gN1MCm06SUB3[setcent],"N1MCm06SUB3 (-0.8 < #eta < -0.4)","p");
    legComparepmMC06->AddEntry(gN1MCp06SUB3[setcent],"N1MCp06SUB3 (0.4 < #eta < 0.8)","p");
    legComparepmMC06->Draw();
    TPaveText * txComparepmMC06 = new TPaveText(0.2, 0.81, 0.49, 0.92, "NDC");
    SetTPaveTxt(txComparepmMC06, 20);
    txComparepmMC06->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txComparepmMC06->AddText(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    txComparepmMC06->Draw();
    cComparepmMC06->Print(Form("figures_MH/EtaGapStudy/Compare_pm_MC06_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");


    TCanvas * cComparepmMC02 = new TCanvas("cComparepmMC02","cComparepmMC02",650,600);
    TPad * padComparepmMC02_1 = (TPad *) cComparepmMC02->cd();
    padComparepmMC02_1->SetGrid();
    TH1D * hComparepmMC02 = new TH1D("hComparepmMC02", "", 100, 0, 4.8);
    hComparepmMC02->SetTitle("");
    hComparepmMC02->SetStats(0);
    hComparepmMC02->GetYaxis()->SetRangeUser(0, 0.04);
    hComparepmMC02->SetXTitle("#Delta#eta");
    hComparepmMC02->SetYTitle("|v_{1}^{even}|");
    hComparepmMC02->Draw();
    gN1MCm02SUB3[setcent]->Draw("same");
    gN1MCp02SUB3[setcent]->Draw("same");
    TLegend * legComparepmMC02 = new TLegend(0.48, 0.82, 0.78, 0.93);
    SetLegend(legComparepmMC02, 18);
    legComparepmMC02->AddEntry(gN1MCm02SUB3[setcent],"N1MCm02SUB3 (-0.4 < #eta < 0.0)","p");
    legComparepmMC02->AddEntry(gN1MCp02SUB3[setcent],"N1MCp02SUB3 (0.0 < #eta < 0.4)","p");
    legComparepmMC02->Draw();
    TPaveText * txComparepmMC02 = new TPaveText(0.2, 0.81, 0.49, 0.92, "NDC");
    SetTPaveTxt(txComparepmMC02, 20);
    txComparepmMC02->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txComparepmMC02->AddText(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    txComparepmMC02->Draw();
    cComparepmMC02->Print(Form("figures_MH/EtaGapStudy/Compare_pm_MC02_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");


    TCanvas * cCompareSUBMCm22 = new TCanvas("cCompareSUBMCm22","cCompareSUBMCm22",650,600);
    TPad * padCompareSUBMCm22 = (TPad *) cCompareSUBMCm22->cd();
    padCompareSUBMCm22->SetGrid();
    TH1D * hCompareSUBMCm22 = new TH1D("hCompareSUBMCm22", "", 100, 0, 4.8);
    hCompareSUBMCm22->SetTitle("");
    hCompareSUBMCm22->SetStats(0);
    hCompareSUBMCm22->GetYaxis()->SetRangeUser(0.0, 0.04);
    hCompareSUBMCm22->SetXTitle("#Delta#eta");
    hCompareSUBMCm22->SetYTitle("|v_{1}^{even}|");
    hCompareSUBMCm22->Draw();
    gN1MCm22SUB2[setcent]->Draw("same");
    gN1MCm22SUB3[setcent]->SetMarkerColor(kBlack);
    gN1MCm22SUB3[setcent]->SetLineColor(kBlack);
    gN1MCm22SUB3[setcent]->SetMarkerStyle(25);
    gN1MCm22SUB3[setcent]->SetMarkerSize(1.3);
    gN1MCm22SUB3[setcent]->Draw("same");
    TLegend * legCompareSUBMCm22 = new TLegend(0.63, 0.82, 0.93, 0.93);
    SetLegend(legCompareSUBMCm22, 18);
    legCompareSUBMCm22->AddEntry(gN1MCm22SUB2[setcent],"N1MCm22SUB2","p");
    legCompareSUBMCm22->AddEntry(gN1MCm22SUB3[setcent],"N1MCm22SUB3","p");
    legCompareSUBMCm22->Draw();
    TPaveText * txCompareSUBMCm22 = new TPaveText(0.2, 0.81, 0.49, 0.92, "NDC");
    SetTPaveTxt(txCompareSUBMCm22, 20);
    txCompareSUBMCm22->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txCompareSUBMCm22->AddText(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    txCompareSUBMCm22->Draw();
    cCompareSUBMCm22->Print(Form("figures_MH/EtaGapStudy/Compare_sub_MCm22_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");


    TCanvas * cCompareSUBMCm18 = new TCanvas("cCompareSUBMCm18","cCompareSUBMCm18",650,600);
    TPad * padCompareSUBMCm18 = (TPad *) cCompareSUBMCm18->cd();
    padCompareSUBMCm18->SetGrid();
    TH1D * hCompareSUBMCm18 = new TH1D("hCompareSUBMCm18", "", 100, 0, 4.8);
    hCompareSUBMCm18->SetTitle("");
    hCompareSUBMCm18->SetStats(0);
    hCompareSUBMCm18->GetYaxis()->SetRangeUser(0.0, 0.04);
    hCompareSUBMCm18->SetXTitle("#Delta#eta");
    hCompareSUBMCm18->SetYTitle("|v_{1}^{even}|");
    hCompareSUBMCm18->Draw();
    gN1MCm18SUB2[setcent]->Draw("same");
    gN1MCm18SUB3[setcent]->SetMarkerColor(kBlack);
    gN1MCm18SUB3[setcent]->SetLineColor(kBlack);
    gN1MCm18SUB3[setcent]->SetMarkerStyle(25);
    gN1MCm18SUB3[setcent]->SetMarkerSize(1.3);
    gN1MCm18SUB3[setcent]->Draw("same");
    TLegend * legCompareSUBMCm18 = new TLegend(0.63, 0.82, 0.93, 0.93);
    SetLegend(legCompareSUBMCm18, 18);
    legCompareSUBMCm18->AddEntry(gN1MCm18SUB2[setcent],"N1MCm18SUB2","p");
    legCompareSUBMCm18->AddEntry(gN1MCm18SUB3[setcent],"N1MCm18SUB3","p");
    legCompareSUBMCm18->Draw();
    TPaveText * txCompareSUBMCm18 = new TPaveText(0.2, 0.81, 0.49, 0.92, "NDC");
    SetTPaveTxt(txCompareSUBMCm18, 20);
    txCompareSUBMCm18->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txCompareSUBMCm18->AddText(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    txCompareSUBMCm18->Draw();
    cCompareSUBMCm18->Print(Form("figures_MH/EtaGapStudy/Compare_sub_MCm18_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");


    TCanvas * cCompareSUBMCm14 = new TCanvas("cCompareSUBMCm14","cCompareSUBMCm14",650,600);
    TPad * padCompareSUBMCm14 = (TPad *) cCompareSUBMCm14->cd();
    padCompareSUBMCm14->SetGrid();
    TH1D * hCompareSUBMCm14 = new TH1D("hCompareSUBMCm14", "", 100, 0, 4.8);
    hCompareSUBMCm14->SetTitle("");
    hCompareSUBMCm14->SetStats(0);
    hCompareSUBMCm14->GetYaxis()->SetRangeUser(0.0, 0.04);
    hCompareSUBMCm14->SetXTitle("#Delta#eta");
    hCompareSUBMCm14->SetYTitle("|v_{1}^{even}|");
    hCompareSUBMCm14->Draw();
    gN1MCm14SUB2[setcent]->Draw("same");
    gN1MCm14SUB3[setcent]->SetMarkerColor(kBlack);
    gN1MCm14SUB3[setcent]->SetLineColor(kBlack);
    gN1MCm14SUB3[setcent]->SetMarkerStyle(25);
    gN1MCm14SUB3[setcent]->SetMarkerSize(1.3);
    gN1MCm14SUB3[setcent]->Draw("same");
    TLegend * legCompareSUBMCm14 = new TLegend(0.63, 0.82, 0.93, 0.93);
    SetLegend(legCompareSUBMCm14, 18);
    legCompareSUBMCm14->AddEntry(gN1MCm14SUB2[setcent],"N1MCm14SUB2","p");
    legCompareSUBMCm14->AddEntry(gN1MCm14SUB3[setcent],"N1MCm14SUB3","p");
    legCompareSUBMCm14->Draw();
    TPaveText * txCompareSUBMCm14 = new TPaveText(0.2, 0.81, 0.49, 0.92, "NDC");
    SetTPaveTxt(txCompareSUBMCm14, 20);
    txCompareSUBMCm14->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txCompareSUBMCm14->AddText(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    txCompareSUBMCm14->Draw();
    cCompareSUBMCm14->Print(Form("figures_MH/EtaGapStudy/Compare_sub_MCm14_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");


    TCanvas * cCompareSUBMCp22 = new TCanvas("cCompareSUBMCp22","cCompareSUBMCp22",650,600);
    TPad * padCompareSUBMCp22 = (TPad *) cCompareSUBMCp22->cd();
    padCompareSUBMCp22->SetGrid();
    TH1D * hCompareSUBMCp22 = new TH1D("hCompareSUBMCp22", "", 100, 0, 4.8);
    hCompareSUBMCp22->SetTitle("");
    hCompareSUBMCp22->SetStats(0);
    hCompareSUBMCp22->GetYaxis()->SetRangeUser(0.0, 0.04);
    hCompareSUBMCp22->SetXTitle("#Delta#eta");
    hCompareSUBMCp22->SetYTitle("|v_{1}^{even}|");
    hCompareSUBMCp22->Draw();
    gN1MCp22SUB2[setcent]->SetMarkerStyle(20);
    gN1MCp22SUB2[setcent]->Draw("same");
    gN1MCp22SUB3[setcent]->SetMarkerColor(kBlack);
    gN1MCp22SUB3[setcent]->SetLineColor(kBlack);
    gN1MCp22SUB3[setcent]->SetMarkerStyle(25);
    gN1MCp22SUB3[setcent]->SetMarkerSize(1.3);
    gN1MCp22SUB3[setcent]->Draw("same");
    TLegend * legCompareSUBMCp22 = new TLegend(0.63, 0.82, 0.93, 0.93);
    SetLegend(legCompareSUBMCp22, 18);
    legCompareSUBMCp22->AddEntry(gN1MCp22SUB2[setcent],"N1MCp22SUB2","p");
    legCompareSUBMCp22->AddEntry(gN1MCp22SUB3[setcent],"N1MCp22SUB3","p");
    legCompareSUBMCp22->Draw();
    TPaveText * txCompareSUBMCp22 = new TPaveText(0.2, 0.81, 0.49, 0.92, "NDC");
    SetTPaveTxt(txCompareSUBMCp22, 20);
    txCompareSUBMCp22->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txCompareSUBMCp22->AddText(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    txCompareSUBMCp22->Draw();
    cCompareSUBMCp22->Print(Form("figures_MH/EtaGapStudy/Compare_sub_MCp22_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");


    TCanvas * cCompareSUBMCp18 = new TCanvas("cCompareSUBMCp18","cCompareSUBMCp18",650,600);
    TPad * padCompareSUBMCp18 = (TPad *) cCompareSUBMCp18->cd();
    padCompareSUBMCp18->SetGrid();
    TH1D * hCompareSUBMCp18 = new TH1D("hCompareSUBMCp18", "", 100, 0, 4.8);
    hCompareSUBMCp18->SetTitle("");
    hCompareSUBMCp18->SetStats(0);
    hCompareSUBMCp18->GetYaxis()->SetRangeUser(0.0, 0.04);
    hCompareSUBMCp18->SetXTitle("#Delta#eta");
    hCompareSUBMCp18->SetYTitle("|v_{1}^{even}|");
    hCompareSUBMCp18->Draw();
    gN1MCp18SUB2[setcent]->Draw("same");
    gN1MCp18SUB3[setcent]->SetMarkerColor(kBlack);
    gN1MCp18SUB3[setcent]->SetLineColor(kBlack);
    gN1MCp18SUB3[setcent]->SetMarkerStyle(25);
    gN1MCp18SUB3[setcent]->SetMarkerSize(1.3);
    gN1MCp18SUB3[setcent]->Draw("same");
    TLegend * legCompareSUBMCp18 = new TLegend(0.63, 0.82, 0.93, 0.93);
    SetLegend(legCompareSUBMCp18, 18);
    legCompareSUBMCp18->AddEntry(gN1MCp18SUB2[setcent],"N1MCp18SUB2","p");
    legCompareSUBMCp18->AddEntry(gN1MCp18SUB3[setcent],"N1MCp18SUB3","p");
    legCompareSUBMCp18->Draw();
    TPaveText * txCompareSUBMCp18 = new TPaveText(0.2, 0.81, 0.49, 0.92, "NDC");
    SetTPaveTxt(txCompareSUBMCp18, 20);
    txCompareSUBMCp18->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txCompareSUBMCp18->AddText(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    txCompareSUBMCp18->Draw();
    cCompareSUBMCp18->Print(Form("figures_MH/EtaGapStudy/Compare_sub_MCp18_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");


    TCanvas * cCompareSUBMCp14 = new TCanvas("cCompareSUBMCp14","cCompareSUBMCp14",650,600);
    TPad * padCompareSUBMCp14 = (TPad *) cCompareSUBMCp14->cd();
    padCompareSUBMCp14->SetGrid();
    TH1D * hCompareSUBMCp14 = new TH1D("hCompareSUBMCp14", "", 100, 0, 4.8);
    hCompareSUBMCp14->SetTitle("");
    hCompareSUBMCp14->SetStats(0);
    hCompareSUBMCp14->GetYaxis()->SetRangeUser(0.0, 0.04);
    hCompareSUBMCp14->SetXTitle("#Delta#eta");
    hCompareSUBMCp14->SetYTitle("|v_{1}^{even}|");
    hCompareSUBMCp14->Draw();
    gN1MCp14SUB2[setcent]->Draw("same");
    gN1MCp14SUB3[setcent]->SetMarkerColor(kBlack);
    gN1MCp14SUB3[setcent]->SetLineColor(kBlack);
    gN1MCp14SUB3[setcent]->SetMarkerStyle(25);
    gN1MCp14SUB3[setcent]->SetMarkerSize(1.3);
    gN1MCp14SUB3[setcent]->Draw("same");
    TLegend * legCompareSUBMCp14 = new TLegend(0.63, 0.82, 0.93, 0.93);
    SetLegend(legCompareSUBMCp14, 18);
    legCompareSUBMCp14->AddEntry(gN1MCp14SUB2[setcent],"N1MCp14SUB2","p");
    legCompareSUBMCp14->AddEntry(gN1MCp14SUB3[setcent],"N1MCp14SUB3","p");
    legCompareSUBMCp14->Draw();
    TPaveText * txCompareSUBMCp14 = new TPaveText(0.2, 0.81, 0.49, 0.92, "NDC");
    SetTPaveTxt(txCompareSUBMCp14, 20);
    txCompareSUBMCp14->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txCompareSUBMCp14->AddText(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    txCompareSUBMCp14->Draw();
    cCompareSUBMCp14->Print(Form("figures_MH/EtaGapStudy/Compare_sub_MCp14_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");


    setcent = 5;
    TCanvas * cComparepmMC22_SUB3 = new TCanvas("cComparepmMC22_SUB3","cComparepmMC22_SUB3",650,600);
    TPad * padComparepmMC22_SUB3_1 = (TPad *) cComparepmMC22_SUB3->cd();
    padComparepmMC22_SUB3_1->SetGrid();
    TH1D * hComparepmMC22_SUB3 = new TH1D("hComparepmMC22_SUB3", "", 100, 0, 4.8);
    hComparepmMC22_SUB3->SetTitle("");
    hComparepmMC22_SUB3->SetStats(0);
    hComparepmMC22_SUB3->GetYaxis()->SetRangeUser(0, 0.04);
    hComparepmMC22_SUB3->SetXTitle("#Delta#eta");
    hComparepmMC22_SUB3->SetYTitle("|v_{1}^{even}|");
    hComparepmMC22_SUB3->Draw();
    gN1MCm22SUB3[setcent]->SetMarkerColor(kRed);
    gN1MCp22SUB3[setcent]->SetMarkerColor(kBlack);
    gN1MCm22SUB3[setcent]->SetLineColor(kRed);
    gN1MCp22SUB3[setcent]->SetLineColor(kBlack);
    gN1MCm22SUB3[setcent]->SetMarkerStyle(20);
    gN1MCp22SUB3[setcent]->SetMarkerStyle(24);
    gN1MCm22SUB3[setcent]->SetMarkerSize(1.3);
    gN1MCp22SUB3[setcent]->SetMarkerSize(1.3);
    gN1MCm22SUB3[setcent]->Draw("same");
    gN1MCp22SUB3[setcent]->Draw("same");
    TLegend * legComparepmMC22_SUB3 = new TLegend(0.48, 0.82, 0.78, 0.93);
    SetLegend(legComparepmMC22_SUB3, 18);
    legComparepmMC22_SUB3->AddEntry(gN1MCm22SUB3[setcent],"N1MCm22SUB3 (-2.4 < #eta < -2.0)","p");
    legComparepmMC22_SUB3->AddEntry(gN1MCp22SUB3[setcent],"N1MCp22SUB3 (2.0 < #eta < 2.4)","p");
    legComparepmMC22_SUB3->Draw();
    TPaveText * txComparepmMC22_SUB3 = new TPaveText(0.2, 0.81, 0.49, 0.92, "NDC");
    SetTPaveTxt(txComparepmMC22_SUB3, 20);
    txComparepmMC22_SUB3->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txComparepmMC22_SUB3->AddText(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    txComparepmMC22_SUB3->Draw();
    cComparepmMC22_SUB3->Print(Form("figures_MH/EtaGapStudy/Compare_pm_MC22_sub3_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");


    TCanvas * cComparepmMC18_SUB3 = new TCanvas("cComparepmMC18_SUB3","cComparepmMC18_SUB3",650,600);
    TPad * padComparepmMC18_SUB3_1 = (TPad *) cComparepmMC18_SUB3->cd();
    padComparepmMC18_SUB3_1->SetGrid();
    TH1D * hComparepmMC18_SUB3 = new TH1D("hComparepmMC18_SUB3", "", 100, 0, 4.8);
    hComparepmMC18_SUB3->SetTitle("");
    hComparepmMC18_SUB3->SetStats(0);
    hComparepmMC18_SUB3->GetYaxis()->SetRangeUser(0, 0.04);
    hComparepmMC18_SUB3->SetXTitle("#Delta#eta");
    hComparepmMC18_SUB3->SetYTitle("|v_{1}^{even}|");
    hComparepmMC18_SUB3->Draw();
    gN1MCm18SUB3[setcent]->SetMarkerColor(kBlue);
    gN1MCp18SUB3[setcent]->SetMarkerColor(kBlack);
    gN1MCm18SUB3[setcent]->SetLineColor(kRed);
    gN1MCp18SUB3[setcent]->SetLineColor(kBlack);
    gN1MCm18SUB3[setcent]->SetMarkerStyle(21);
    gN1MCp18SUB3[setcent]->SetMarkerStyle(25);
    gN1MCm18SUB3[setcent]->SetMarkerSize(1.2);
    gN1MCp18SUB3[setcent]->SetMarkerSize(1.2);
    gN1MCm18SUB3[setcent]->Draw("same");
    gN1MCp18SUB3[setcent]->Draw("same");
    TLegend * legComparepmMC18_SUB3 = new TLegend(0.48, 0.82, 0.78, 0.93);
    SetLegend(legComparepmMC18_SUB3, 18);
    legComparepmMC18_SUB3->AddEntry(gN1MCm18SUB3[setcent],"N1MCm18SUB3 (-2.0 < #eta < -1.6)","p");
    legComparepmMC18_SUB3->AddEntry(gN1MCp18SUB3[setcent],"N1MCp18SUB3 (1.6 < #eta < 2.0)","p");
    legComparepmMC18_SUB3->Draw();
    TPaveText * txComparepmMC18_SUB3 = new TPaveText(0.2, 0.81, 0.49, 0.92, "NDC");
    SetTPaveTxt(txComparepmMC18_SUB3, 20);
    txComparepmMC18_SUB3->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txComparepmMC18_SUB3->AddText(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    txComparepmMC18_SUB3->Draw();
    cComparepmMC18_SUB3->Print(Form("figures_MH/EtaGapStudy/Compare_pm_MC18_sub3_%d_%d.png",cminCENT[setcent],cmaxCENT[setcent]),"png");

}
