# include "TArrow.h"
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
                if (i == 0 && cbin == 0) cout<<"neta = "<<neta<<"\tdeltaEta = "<<deltaEta[neta]<<"\tv1A = "<<v1A<<"\tv1Ae = "<<v1Ae<<endl;
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
                    gN1MCm02SUB2[cbin]->SetBinContent(ebin+1, v1val[ebin+i]);
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
                    if (cbin == 0) cout<<"i: "<<i<<"\tebin: "<<ebin<<"\tdeltaEta: "<<deltaEta[ebin]<<"\tv1val: "<<v1val[neta-ebin-2]<<endl;
                } else if (i == 11) {
                    gN1MCp22SUB2[cbin]->SetBinContent(ebin+1, v1val[neta-ebin-1]);
                    gN1MCp22SUB2[cbin]->SetBinError(ebin+1, v1vale[neta-ebin-1]);
                } else {continue;}
            }


            /*for (int j = 0; j<neta; j++) {
                cout<<"deltaEta = "<<deltaEta[j]<<"\tv1val = "<<v1val[j]<<"\tv1vale = "<<v1vale[j]<<endl;
            }*/
            // for (int j = 0; j<neta; j++) {
            //     if (i == 11 && cbin == 0) {
            //         cout<<"deltaEta: "<<deltaEta[j]<<"\tv1val: "<<v1val[j]<<"\tv1vale: "<<v1vale[j]<<endl;
            //     }
            // }
        }
    }

    // draw options
    for (int cbin = 0; cbin<cbinsCENT; cbin++) {
        gN1MCm22SUB2[cbin]->SetMarkerColor(kBlue);
        gN1MCm18SUB2[cbin]->SetMarkerColor(kRed);
        gN1MCm14SUB2[cbin]->SetMarkerColor(kGreen+2);
        gN1MCm10SUB2[cbin]->SetMarkerColor(kMagenta);
        gN1MCm06SUB2[cbin]->SetMarkerColor(kOrange+7);
        gN1MCm02SUB2[cbin]->SetMarkerColor(kBlack);
        gN1MCp02SUB2[cbin]->SetMarkerColor(kBlack);
        gN1MCp06SUB2[cbin]->SetMarkerColor(kOrange+7);
        gN1MCp10SUB2[cbin]->SetMarkerColor(kMagenta);
        gN1MCp14SUB2[cbin]->SetMarkerColor(kGreen+2);
        gN1MCp18SUB2[cbin]->SetMarkerColor(kRed);
        gN1MCp22SUB2[cbin]->SetMarkerColor(kBlue);

        gN1MCm22SUB2[cbin]->SetLineColor(kBlue);
        gN1MCm18SUB2[cbin]->SetLineColor(kRed);
        gN1MCm14SUB2[cbin]->SetLineColor(kGreen+2);
        gN1MCm10SUB2[cbin]->SetLineColor(kMagenta);
        gN1MCm06SUB2[cbin]->SetLineColor(kOrange+7);
        gN1MCm02SUB2[cbin]->SetLineColor(kBlack);
        gN1MCp02SUB2[cbin]->SetLineColor(kBlack);
        gN1MCp06SUB2[cbin]->SetLineColor(kOrange+7);
        gN1MCp10SUB2[cbin]->SetLineColor(kMagenta);
        gN1MCp14SUB2[cbin]->SetLineColor(kGreen+2);
        gN1MCp18SUB2[cbin]->SetLineColor(kRed);
        gN1MCp22SUB2[cbin]->SetLineColor(kBlue);

        gN1MCm22SUB2[cbin]->SetMarkerStyle(21);
        gN1MCm18SUB2[cbin]->SetMarkerStyle(20);
        gN1MCm14SUB2[cbin]->SetMarkerStyle(34);
        gN1MCm10SUB2[cbin]->SetMarkerStyle(33);
        gN1MCm06SUB2[cbin]->SetMarkerStyle(21);
        gN1MCm02SUB2[cbin]->SetMarkerStyle(34);
        gN1MCp02SUB2[cbin]->SetMarkerStyle(28);
        gN1MCp06SUB2[cbin]->SetMarkerStyle(25);
        gN1MCp10SUB2[cbin]->SetMarkerStyle(27);
        gN1MCp14SUB2[cbin]->SetMarkerStyle(28);
        gN1MCp18SUB2[cbin]->SetMarkerStyle(24);
        gN1MCp22SUB2[cbin]->SetMarkerStyle(25);

        gN1MCm22SUB2[cbin]->SetMarkerSize(1.2);
        gN1MCm18SUB2[cbin]->SetMarkerSize(1.3);
        gN1MCm14SUB2[cbin]->SetMarkerSize(1.7);
        gN1MCm10SUB2[cbin]->SetMarkerSize(1.8);
        gN1MCm06SUB2[cbin]->SetMarkerSize(1.2);
        gN1MCm02SUB2[cbin]->SetMarkerSize(1.7);
        gN1MCp02SUB2[cbin]->SetMarkerSize(1.7);
        gN1MCp06SUB2[cbin]->SetMarkerSize(1.2);
        gN1MCp10SUB2[cbin]->SetMarkerSize(1.8);
        gN1MCp14SUB2[cbin]->SetMarkerSize(1.7);
        gN1MCp18SUB2[cbin]->SetMarkerSize(1.3);
        gN1MCp22SUB2[cbin]->SetMarkerSize(1.2);
    }


    if (!fopen("figures_MH/EtaGapStudy/","r")) system("mkdir figures_MH/EtaGapStudy");

    int setcent = 12;

    TCanvas * cMCmScan = new TCanvas("cMCmScan","cMCmScan",650,600);
    TPad * padMCmScan = (TPad *) cMCmScan->cd();
    padMCmScan->SetGrid();
    TH1D * hMCmScan = new TH1D("hMCmScan", "hMCmScan", 100, 0, 4.8);
    hMCmScan->SetTitle("");
    hMCmScan->SetStats(0);
    hMCmScan->GetYaxis()->SetRangeUser(0, 0.05);
    hMCmScan->SetXTitle("#Delta#eta");
    hMCmScan->SetYTitle("|v_{1}^{even}|");
    hMCmScan->Draw();
    gN1MCm22SUB2[setcent-1]->Draw("same p");
    gN1MCm18SUB2[setcent-1]->Draw("same p");
    gN1MCm14SUB2[setcent-1]->Draw("same p");
    gN1MCm10SUB2[setcent-1]->Draw("same p");
    gN1MCm06SUB2[setcent-1]->Draw("same p");
    gN1MCm02SUB2[setcent-1]->Draw("same p");
    TLegend * legMCmScan = new TLegend(0.66, 0.67, 0.93, 0.93);
    SetLegend(legMCmScan,18);
    legMCmScan->AddEntry(gN1MCm22SUB2[setcent-1],"N1MCm22SUB2","p");
    legMCmScan->AddEntry(gN1MCm18SUB2[setcent-1],"N1MCm18SUB2","p");
    legMCmScan->AddEntry(gN1MCm14SUB2[setcent-1],"N1MCm14SUB2","p");
    legMCmScan->AddEntry(gN1MCm10SUB2[setcent-1],"N1MCm10SUB2","p");
    legMCmScan->AddEntry(gN1MCm06SUB2[setcent-1],"N1MCm06SUB2","p");
    legMCmScan->AddEntry(gN1MCm02SUB2[setcent-1],"N1MCm02SUB2","p");
    legMCmScan->Draw();
    TPaveText * txMCmScan = new TPaveText(0.2, 0.81, 0.47, 0.92, "NDC");
    SetTPaveTxt(txMCmScan, 20);
    txMCmScan->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txMCmScan->AddText(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    txMCmScan->Draw();
    cMCmScan->Print(Form("figures_MH/EtaGapStudy/MCm_deltaEta_scan_%d_%d.pdf",cminCENT[setcent],cmaxCENT[setcent]),"pdf");


    TCanvas * cMCpScan = new TCanvas("cMCpScan","cMCpScan",650,600);
    TPad * padMCpScan = (TPad *) cMCpScan->cd();
    padMCpScan->SetGrid();
    TH1D * hMCpScan = new TH1D("hMCpScan", "hMCpScan", 100, 0, 4.8);
    hMCpScan->SetTitle("");
    hMCpScan->SetStats(0);
    hMCpScan->GetYaxis()->SetRangeUser(0, 0.05);
    hMCpScan->SetXTitle("#Delta#eta");
    hMCpScan->SetYTitle("|v_{1}^{even}|");
    hMCpScan->Draw();
    gN1MCp22SUB2[setcent-1]->Draw("same p");
    gN1MCp18SUB2[setcent-1]->Draw("same p");
    gN1MCp14SUB2[setcent-1]->Draw("same p");
    gN1MCp10SUB2[setcent-1]->Draw("same p");
    gN1MCp06SUB2[setcent-1]->Draw("same p");
    gN1MCp02SUB2[setcent-1]->Draw("same p");
    TLegend * legMCpScan = new TLegend(0.66, 0.67, 0.93, 0.93);
    SetLegend(legMCpScan,18);
    legMCpScan->AddEntry(gN1MCp22SUB2[setcent-1],"N1MCp22SUB2","p");
    legMCpScan->AddEntry(gN1MCp18SUB2[setcent-1],"N1MCp18SUB2","p");
    legMCpScan->AddEntry(gN1MCp14SUB2[setcent-1],"N1MCp14SUB2","p");
    legMCpScan->AddEntry(gN1MCp10SUB2[setcent-1],"N1MCp10SUB2","p");
    legMCpScan->AddEntry(gN1MCp06SUB2[setcent-1],"N1MCp06SUB2","p");
    legMCpScan->AddEntry(gN1MCp02SUB2[setcent-1],"N1MCp02SUB2","p");
    legMCpScan->Draw();
    TPaveText * txMCpScan = new TPaveText(0.2, 0.81, 0.47, 0.92, "NDC");
    SetTPaveTxt(txMCpScan, 20);
    txMCpScan->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txMCpScan->AddText(Form("%d-%d%%",cminCENT[setcent],cmaxCENT[setcent]));
    txMCpScan->Draw();
    cMCpScan->Print(Form("figures_MH/EtaGapStudy/MCp_deltaEta_scan_%d_%d.pdf",cminCENT[setcent],cmaxCENT[setcent]),"pdf");

}
