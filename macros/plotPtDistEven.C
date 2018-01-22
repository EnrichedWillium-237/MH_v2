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

static const int ncentbins = 13;
static const int cminCENT[] = {0,  5, 10, 15, 20, 25, 30, 35, 40, 50, 60,  0, 20,  60};
static const int cmaxCENT[] = {5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 20, 60, 100};
static const int netabins = 12;
static const double etabins[] = {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0,  2.4};
static const double etaMid[] = {-2.2, -1.8, -1.4, -1.0, -0.6, -0.2,  0.2,  0.6,  1.0,  1.4,  1.8,  2.2};
TString etags[] = {"-2.4", "-2.0", "-1.6", "-1.2", "-0.8", "-0.4", "00.0", "00.4", "00.8", "01.2", "01.6", "02.0", "02.4"};
static const int nptbins = 18;
static const double ptbins[] = {0.30,  0.40,  0.50,  0.60,  0.80,  1.00,  1.25,  1.50,  2.00,  2.50,  3.00,
                     3.50,  4.00,  5.00,  6.00,  7.00,  8.00,  10.00,  12.00};
static int NANALS = 24;
TString ANAL[] = {"N1MCm22SUB2", "N1MCm18SUB2", "N1MCm14SUB2", "N1MCm10SUB2", "N1MCm06SUB2", "N1MCm02SUB2",
                 "N1MCp02SUB2", "N1MCp06SUB2", "N1MCp10SUB2", "N1MCp14SUB2", "N1MCp18SUB2", "N1MCp22SUB2",
                 "N1MCm22SUB3", "N1MCm18SUB3", "N1MCm14SUB3", "N1MCm10SUB3", "N1MCm06SUB3", "N1MCm02SUB3",
                 "N1MCp02SUB3", "N1MCp06SUB3", "N1MCp10SUB3", "N1MCp14SUB3", "N1MCp18SUB3", "N1MCp22SUB3"
};

TH1D * gN1MCm22SUB3_pt[netabins][ncentbins];
TH1D * gN1MCm18SUB3_pt[netabins][ncentbins];
TH1D * gN1MCm14SUB3_pt[netabins][ncentbins];

TH1D * gN1MCp22SUB3_pt[netabins][ncentbins];
TH1D * gN1MCp18SUB3_pt[netabins][ncentbins];
TH1D * gN1MCp14SUB3_pt[netabins][ncentbins];

TH1D * gN1MC22SUB3_pt[netabins][ncentbins];
TH1D * gN1MC18SUB3_pt[netabins][ncentbins];
TH1D * gN1MC14SUB3_pt[netabins][ncentbins];

void plotPtDistEven()
{

    TH1::SetDefaultSumw2();

    for (int ebin = 0; ebin<netabins; ebin++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            gN1MCm22SUB3_pt[ebin][cbin] = new TH1D(Form("N1MCm22SUB3_e%d_c%d_%d",ebin,cminCENT[cbin],cmaxCENT[cbin]), "", nptbins, ptbins);
            gN1MCm18SUB3_pt[ebin][cbin] = new TH1D(Form("N1MCm18SUB3_e%d_c%d_%d",ebin,cminCENT[cbin],cmaxCENT[cbin]), "", nptbins, ptbins);
            gN1MCm14SUB3_pt[ebin][cbin] = new TH1D(Form("N1MCm14SUB3_e%d_c%d_%d",ebin,cminCENT[cbin],cmaxCENT[cbin]), "", nptbins, ptbins);

            gN1MCp22SUB3_pt[ebin][cbin] = new TH1D(Form("N1MCp22SUB3_e%d_c%d_%d",ebin,cminCENT[cbin],cmaxCENT[cbin]), "", nptbins, ptbins);
            gN1MCp18SUB3_pt[ebin][cbin] = new TH1D(Form("N1MCp18SUB3_e%d_c%d_%d",ebin,cminCENT[cbin],cmaxCENT[cbin]), "", nptbins, ptbins);
            gN1MCp14SUB3_pt[ebin][cbin] = new TH1D(Form("N1MCp14SUB3_e%d_c%d_%d",ebin,cminCENT[cbin],cmaxCENT[cbin]), "", nptbins, ptbins);
        }
    }

    for (int i = 0; i<NANALS; i++) {
        for (int ebin = 0; ebin<netabins; ebin++) {
            for (int cbin = 0; cbin<ncentbins; cbin++) {
                TString tag = Form("figures_MH/%s/eta_%s_%s/data/%s_%d_%d_eta_%1.1f_%1.1f_A.dat",ANAL[i].Data(),etags[ebin].Data(),etags[ebin+1].Data(),ANAL[i].Data(),cminCENT[cbin],cmaxCENT[cbin],etabins[ebin],etabins[ebin+1]);
                // cout<<tag<<endl;
                ifstream fin(tag.Data());

                double pt, v1, v1e;
                int npt = 0;
                while (fin >> pt >> v1 >> v1e) {
                    // if (i == 0 && ebin == 0 && cbin == 0) cout<<"pt: "<<pt<<"\tv1: "<<v1<<"\tv1e: "<<v1e<<endl;

                    if (i == 12) {
                        gN1MCm22SUB3_pt[ebin][cbin]->SetBinContent(npt+1, v1);
                        gN1MCm22SUB3_pt[ebin][cbin]->SetBinError(npt+1, v1e);
                    } else if (i == 13) {
                        gN1MCm18SUB3_pt[ebin][cbin]->SetBinContent(npt+1, v1);
                        gN1MCm18SUB3_pt[ebin][cbin]->SetBinError(npt+1, v1e);
                    } else if (i == 14) {
                        gN1MCm14SUB3_pt[ebin][cbin]->SetBinContent(npt+1, v1);
                        gN1MCm14SUB3_pt[ebin][cbin]->SetBinError(npt+1, v1e);
                    } else if (i == 21) {
                        gN1MCp14SUB3_pt[ebin][cbin]->SetBinContent(npt+1, v1);
                        gN1MCp14SUB3_pt[ebin][cbin]->SetBinError(npt+1, v1e);
                    } else if (i == 22) {
                        gN1MCp18SUB3_pt[ebin][cbin]->SetBinContent(npt+1, v1);
                        gN1MCp18SUB3_pt[ebin][cbin]->SetBinError(npt+1, v1e);
                    } else if (i == 23) {
                        gN1MCp22SUB3_pt[ebin][cbin]->SetBinContent(npt+1, v1);
                        gN1MCp22SUB3_pt[ebin][cbin]->SetBinError(npt+1, v1e);
                    } else {continue;}

                    npt++;
                } // end while loop
            } // end cent loop
        } // end eta loop
    } // end analysis loop

    if (!fopen("figures_MH/PtDists","r")) system("mkdir figures_MH/PtDists");

    // find average v1 over +/- eta ranges
    for (int ebin = 0; ebin<netabins; ebin++) {
        if (!fopen(Form("figures_MH/PtDists/eta_%1.1f_%1.1f",etabins[ebin],etabins[ebin+1]),"r")) system(Form("mkdir figures_MH/PtDists/eta_%1.1f_%1.1f",etabins[ebin],etabins[ebin+1]));
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            gN1MC22SUB3_pt[ebin][cbin] = (TH1D *) gN1MCm22SUB3_pt[ebin][cbin]->Clone(Form("N1MC22SUB3_e%d_c%d_%d",ebin,cminCENT[cbin],cmaxCENT[cbin]));
            gN1MC22SUB3_pt[ebin][cbin]->Add(gN1MCp22SUB3_pt[ebin][cbin]);
            gN1MC22SUB3_pt[ebin][cbin]->Scale(0.5);

            gN1MC18SUB3_pt[ebin][cbin] = (TH1D *) gN1MCm18SUB3_pt[ebin][cbin]->Clone(Form("N1MC18SUB3_e%d_c%d_%d",ebin,cminCENT[cbin],cmaxCENT[cbin]));
            gN1MC18SUB3_pt[ebin][cbin]->Add(gN1MCp18SUB3_pt[ebin][cbin]);
            gN1MC18SUB3_pt[ebin][cbin]->Scale(0.5);

            gN1MC14SUB3_pt[ebin][cbin] = (TH1D *) gN1MCm14SUB3_pt[ebin][cbin]->Clone(Form("N1MC14SUB3_e%d_c%d_%d",ebin,cminCENT[cbin],cmaxCENT[cbin]));
            gN1MC14SUB3_pt[ebin][cbin]->Add(gN1MCp14SUB3_pt[ebin][cbin]);
            gN1MC14SUB3_pt[ebin][cbin]->Scale(0.5);
        }
    }

    int setcent(5);

    TH1D * hPtDummy = new TH1D("hPtDummy", "", 100, 0, 12);
    hPtDummy->SetStats(0);
    hPtDummy->SetXTitle("p_{T} (GeV/c)");
    hPtDummy->SetYTitle("v_{1}^{even}");
    hPtDummy->GetYaxis()->SetRangeUser(-0.04, 0.18);


    ///-- minus side comparison

    TCanvas * cPtDist_MCm22_MCm18[netabins];
    for (int ebin = 0; ebin<netabins; ebin++) {
        cPtDist_MCm22_MCm18[ebin] = new TCanvas(Form("cPtDist_MCm22_MCm18_%d",ebin),"",650,600);
        TPad * padPtDist_MCm22_MCm18 = (TPad *) cPtDist_MCm22_MCm18[ebin]->cd();
        padPtDist_MCm22_MCm18->SetGrid();
        TH1D * hPtDist_MCm22_MCm18 = (TH1D *) hPtDummy->Clone(Form("hPtDist_MCm22_MCm18_%d",ebin));
        hPtDist_MCm22_MCm18->Draw();
        gN1MCm22SUB3_pt[ebin][setcent]->SetMarkerColor(kRed);
        gN1MCm22SUB3_pt[ebin][setcent]->SetLineColor(kRed);
        gN1MCm22SUB3_pt[ebin][setcent]->SetMarkerStyle(20);
        gN1MCm22SUB3_pt[ebin][setcent]->SetMarkerSize(1.3);
        gN1MCm22SUB3_pt[ebin][setcent]->Draw("same");
        gN1MCm18SUB3_pt[ebin][setcent]->SetMarkerColor(kBlue);
        gN1MCm18SUB3_pt[ebin][setcent]->SetLineColor(kBlue);
        gN1MCm18SUB3_pt[ebin][setcent]->SetMarkerStyle(21);
        gN1MCm18SUB3_pt[ebin][setcent]->SetMarkerSize(1.2);
        gN1MCm18SUB3_pt[ebin][setcent]->Draw("same");
        TLegend * legPtDist_MCm22_MCm18 = new TLegend(0.59, 0.18, 0.87, 0.28);
        SetLegend(legPtDist_MCm22_MCm18, 18);
        legPtDist_MCm22_MCm18->AddEntry(gN1MCm22SUB3_pt[ebin][setcent],"#Psi_{1}^{trk}: -2.4 < #eta < -2.0","p");
        legPtDist_MCm22_MCm18->AddEntry(gN1MCm18SUB3_pt[ebin][setcent],"#Psi_{1}^{trk}: -2.0 < #eta < -1.6","p");
        legPtDist_MCm22_MCm18->Draw();
        TPaveText * txPtDist_MCm22_MCm18 = new TPaveText(0.59, 0.30, 0.88, 0.40, "NDC");
        SetTPaveTxt(txPtDist_MCm22_MCm18, 18);
        txPtDist_MCm22_MCm18->AddText("PbPb #sqrt{s_{NN}}=5.02 TeV");
        txPtDist_MCm22_MCm18->AddText(Form("%1.1f < #eta < %1.1f (%d - %d%%)",etabins[ebin],etabins[ebin+1],cminCENT[setcent],cmaxCENT[setcent]));
        txPtDist_MCm22_MCm18->Draw();
        cPtDist_MCm22_MCm18[ebin]->Print(Form("figures_MH/PtDists/eta_%1.1f_%1.1f/MCm22_MCm18_%d_%d.png",etabins[ebin],etabins[ebin+1],cminCENT[setcent],cmaxCENT[setcent]),"png");
        cPtDist_MCm22_MCm18[ebin]->Close();
    }


    TCanvas * crat_MCm22_MCm18[netabins];
    TH1D * rat_MCm22_MCm18[netabins];
    for (int ebin = 0; ebin<netabins; ebin++) {
        crat_MCm22_MCm18[ebin] = new TCanvas(Form("crat_MCm22_MCm18_%d",ebin),"",800,300);
        TPad * padrat_MCm22_MCm18 = (TPad *) crat_MCm22_MCm18[ebin]->cd();
        padrat_MCm22_MCm18->SetGrid();
        TH1D * hrat_MCm22_MCm18 = (TH1D *) hPtDummy->Clone(Form("hrat_MCm22_MCm18_%d",ebin));
        hrat_MCm22_MCm18->GetYaxis()->SetRangeUser(-1.0, 3.0);
        hrat_MCm22_MCm18->GetYaxis()->SetNdivisions(507);
        hrat_MCm22_MCm18->SetXTitle("p_{T} (GeV/c)");
        hrat_MCm22_MCm18->SetYTitle("N1MCm18SUB3 / N1MCm22SUB3");
        hrat_MCm22_MCm18->GetXaxis()->SetTitleSize(0.08);
        hrat_MCm22_MCm18->GetXaxis()->SetTitleOffset(0.90);
        hrat_MCm22_MCm18->GetXaxis()->SetLabelSize(0.07);
        hrat_MCm22_MCm18->GetYaxis()->SetTitleSize(0.07);
        hrat_MCm22_MCm18->GetYaxis()->SetTitleOffset(0.80);
        hrat_MCm22_MCm18->GetYaxis()->SetLabelSize(0.08);
        hrat_MCm22_MCm18->Draw();
        rat_MCm22_MCm18[ebin] = (TH1D *) gN1MCm18SUB3_pt[ebin][setcent]->Clone(Form("rat_MCm22_MCm18_%d",ebin));
        rat_MCm22_MCm18[ebin]->Divide(gN1MCm22SUB3_pt[ebin][setcent]);
        rat_MCm22_MCm18[ebin]->SetMarkerColor(kBlack);
        rat_MCm22_MCm18[ebin]->SetLineColor(kBlack);
        rat_MCm22_MCm18[ebin]->SetMarkerStyle(20);
        rat_MCm22_MCm18[ebin]->SetMarkerSize(1.2);
        for (int i = 1; i<=rat_MCm22_MCm18[ebin]->GetNbinsX(); i++) {
            double xin = gN1MCm18SUB3_pt[ebin][setcent]->GetBinContent(i);
            double yin = gN1MCm22SUB3_pt[ebin][setcent]->GetBinContent(i);
            double delxin = gN1MCm18SUB3_pt[ebin][setcent]->GetBinError(i);
            double delyin = gN1MCm22SUB3_pt[ebin][setcent]->GetBinError(i);
            double raterr = ErrorCalc( xin, yin, delxin, delyin );
            rat_MCm22_MCm18[ebin]->SetBinError(i, raterr);
        }
        rat_MCm22_MCm18[ebin]->Draw("same");
        TLine * lnrat_MCm22_MCm18 = new TLine(0.0, 1.0, 12.0, 1.0);
        lnrat_MCm22_MCm18->SetLineWidth(2);
        lnrat_MCm22_MCm18->Draw();
        crat_MCm22_MCm18[ebin]->Print(Form("figures_MH/PtDists/eta_%1.1f_%1.1f/ratio_MCm22_MCm18_%d_%d.png",etabins[ebin],etabins[ebin+1],cminCENT[setcent],cmaxCENT[setcent]),"png");
        crat_MCm22_MCm18[ebin]->Close();
    }


    TCanvas * cPtDist_MCm18_MCm14[netabins];
    for (int ebin = 0; ebin<netabins; ebin++) {
        cPtDist_MCm18_MCm14[ebin] = new TCanvas(Form("cPtDist_MCm18_MCm14_%d",ebin),"",650,600);
        TPad * padPtDist_MCm18_MCm14 = (TPad *) cPtDist_MCm18_MCm14[ebin]->cd();
        padPtDist_MCm18_MCm14->SetGrid();
        TH1D * hPtDist_MCm18_MCm14 = (TH1D *) hPtDummy->Clone(Form("hPtDist_MCm18_MCm14_%d",ebin));
        hPtDist_MCm18_MCm14->Draw();
        gN1MCm18SUB3_pt[ebin][setcent]->SetMarkerColor(kBlue);
        gN1MCm18SUB3_pt[ebin][setcent]->SetLineColor(kBlue);
        gN1MCm18SUB3_pt[ebin][setcent]->SetMarkerStyle(21);
        gN1MCm18SUB3_pt[ebin][setcent]->SetMarkerSize(1.2);
        gN1MCm18SUB3_pt[ebin][setcent]->Draw("same");
        gN1MCm14SUB3_pt[ebin][setcent]->SetMarkerColor(kGreen+2);
        gN1MCm14SUB3_pt[ebin][setcent]->SetLineColor(kGreen+2);
        gN1MCm14SUB3_pt[ebin][setcent]->SetMarkerStyle(20);
        gN1MCm14SUB3_pt[ebin][setcent]->SetMarkerSize(1.3);
        gN1MCm14SUB3_pt[ebin][setcent]->Draw("same");
        TLegend * legPtDist_MCm18_MCm14 = new TLegend(0.59, 0.18, 0.87, 0.28);
        SetLegend(legPtDist_MCm18_MCm14, 18);
        legPtDist_MCm18_MCm14->AddEntry(gN1MCm18SUB3_pt[ebin][setcent],"#Psi_{1}^{trk}: -2.0 < #eta < -1.6","p");
        legPtDist_MCm18_MCm14->AddEntry(gN1MCm14SUB3_pt[ebin][setcent],"#Psi_{1}^{trk}: -1.6 < #eta < -1.2","p");
        legPtDist_MCm18_MCm14->Draw();
        TPaveText * txPtDist_MCm18_MCm14 = new TPaveText(0.59, 0.30, 0.88, 0.40, "NDC");
        SetTPaveTxt(txPtDist_MCm18_MCm14, 18);
        txPtDist_MCm18_MCm14->AddText("PbPb #sqrt{s_{NN}}=5.02 TeV");
        txPtDist_MCm18_MCm14->AddText(Form("%1.1f < #eta < %1.1f (%d - %d%%)",etabins[ebin],etabins[ebin+1],cminCENT[setcent],cmaxCENT[setcent]));
        txPtDist_MCm18_MCm14->Draw();
        cPtDist_MCm18_MCm14[ebin]->Print(Form("figures_MH/PtDists/eta_%1.1f_%1.1f/MCm18_MCm14_%d_%d.png",etabins[ebin],etabins[ebin+1],cminCENT[setcent],cmaxCENT[setcent]),"png");
        cPtDist_MCm18_MCm14[ebin]->Close();
    }


    TCanvas * crat_MCm18_MCm14[netabins];
    TH1D * rat_MCm18_MCm14[netabins];
    for (int ebin = 0; ebin<netabins; ebin++) {
        crat_MCm18_MCm14[ebin] = new TCanvas(Form("crat_MCm18_MCm14_%d",ebin),"",800,300);
        TPad * padrat_MCm18_MCm14 = (TPad *) crat_MCm18_MCm14[ebin]->cd();
        padrat_MCm18_MCm14->SetGrid();
        TH1D * hrat_MCm18_MCm14 = (TH1D *) hPtDummy->Clone(Form("hrat_MCm18_MCm14_%d",ebin));
        hrat_MCm18_MCm14->GetYaxis()->SetRangeUser(0.0, 2.0);
        hrat_MCm18_MCm14->GetYaxis()->SetNdivisions(507);
        hrat_MCm18_MCm14->SetXTitle("p_{T} (GeV/c)");
        hrat_MCm18_MCm14->SetYTitle("N1MCm14SUB3 / N1MCm18SUB3");
        hrat_MCm18_MCm14->GetXaxis()->SetTitleSize(0.08);
        hrat_MCm18_MCm14->GetXaxis()->SetTitleOffset(0.90);
        hrat_MCm18_MCm14->GetXaxis()->SetLabelSize(0.07);
        hrat_MCm18_MCm14->GetYaxis()->SetTitleSize(0.07);
        hrat_MCm18_MCm14->GetYaxis()->SetTitleOffset(0.80);
        hrat_MCm18_MCm14->GetYaxis()->SetLabelSize(0.08);
        hrat_MCm18_MCm14->Draw();
        rat_MCm18_MCm14[ebin] = (TH1D *) gN1MCm14SUB3_pt[ebin][setcent]->Clone(Form("rat_MCm18_MCm14_%d",ebin));
        rat_MCm18_MCm14[ebin]->Divide(gN1MCm18SUB3_pt[ebin][setcent]);
        rat_MCm18_MCm14[ebin]->SetMarkerColor(kBlack);
        rat_MCm18_MCm14[ebin]->SetLineColor(kBlack);
        rat_MCm18_MCm14[ebin]->SetMarkerStyle(20);
        rat_MCm18_MCm14[ebin]->SetMarkerSize(1.2);
        for (int i = 1; i<=rat_MCm18_MCm14[ebin]->GetNbinsX(); i++) {
            double xin = gN1MCm14SUB3_pt[ebin][setcent]->GetBinContent(i);
            double yin = gN1MCm18SUB3_pt[ebin][setcent]->GetBinContent(i);
            double delxin = gN1MCm14SUB3_pt[ebin][setcent]->GetBinError(i);
            double delyin = gN1MCm18SUB3_pt[ebin][setcent]->GetBinError(i);
            double raterr = ErrorCalc( xin, yin, delxin, delyin );
            rat_MCm18_MCm14[ebin]->SetBinError(i, raterr);
        }
        rat_MCm18_MCm14[ebin]->Draw("same");
        TLine * lnrat_MCm18_MCm14 = new TLine(0.0, 1.0, 12.0, 1.0);
        lnrat_MCm18_MCm14->SetLineWidth(2);
        lnrat_MCm18_MCm14->Draw();
        crat_MCm18_MCm14[ebin]->Print(Form("figures_MH/PtDists/eta_%1.1f_%1.1f/ratio_MCm18_MCm14_%d_%d.png",etabins[ebin],etabins[ebin+1],cminCENT[setcent],cmaxCENT[setcent]),"png");
        crat_MCm18_MCm14[ebin]->Close();
    }


    ///-- plus side comparison


    TCanvas * cPtDist_MCp22_MCp18[netabins];
    for (int ebin = 0; ebin<netabins; ebin++) {
        cPtDist_MCp22_MCp18[ebin] = new TCanvas(Form("cPtDist_MCp22_MCp18_%d",ebin),"",650,600);
        TPad * padPtDist_MCp22_MCp18 = (TPad *) cPtDist_MCp22_MCp18[ebin]->cd();
        padPtDist_MCp22_MCp18->SetGrid();
        TH1D * hPtDist_MCp22_MCp18 = (TH1D *) hPtDummy->Clone(Form("hPtDist_MCp22_MCp18_%d",ebin));
        hPtDist_MCp22_MCp18->Draw();
        gN1MCp22SUB3_pt[ebin][setcent]->SetMarkerColor(kRed);
        gN1MCp22SUB3_pt[ebin][setcent]->SetLineColor(kRed);
        gN1MCp22SUB3_pt[ebin][setcent]->SetMarkerStyle(20);
        gN1MCp22SUB3_pt[ebin][setcent]->SetMarkerSize(1.3);
        gN1MCp22SUB3_pt[ebin][setcent]->Draw("same");
        gN1MCp18SUB3_pt[ebin][setcent]->SetMarkerColor(kBlue);
        gN1MCp18SUB3_pt[ebin][setcent]->SetLineColor(kBlue);
        gN1MCp18SUB3_pt[ebin][setcent]->SetMarkerStyle(21);
        gN1MCp18SUB3_pt[ebin][setcent]->SetMarkerSize(1.2);
        gN1MCp18SUB3_pt[ebin][setcent]->Draw("same");
        TLegend * legPtDist_MCp22_MCp18 = new TLegend(0.59, 0.18, 0.87, 0.28);
        SetLegend(legPtDist_MCp22_MCp18, 18);
        legPtDist_MCp22_MCp18->AddEntry(gN1MCp22SUB3_pt[ebin][setcent],"#Psi_{1}^{trk}: 2.0 < #eta < 2.4","p");
        legPtDist_MCp22_MCp18->AddEntry(gN1MCp18SUB3_pt[ebin][setcent],"#Psi_{1}^{trk}: 1.6 < #eta < 2.0","p");
        legPtDist_MCp22_MCp18->Draw();
        TPaveText * txPtDist_MCp22_MCp18 = new TPaveText(0.59, 0.30, 0.88, 0.40, "NDC");
        SetTPaveTxt(txPtDist_MCp22_MCp18, 18);
        txPtDist_MCp22_MCp18->AddText("PbPb #sqrt{s_{NN}}=5.02 TeV");
        txPtDist_MCp22_MCp18->AddText(Form("%1.1f < #eta < %1.1f (%d - %d%%)",etabins[ebin],etabins[ebin+1],cminCENT[setcent],cmaxCENT[setcent]));
        txPtDist_MCp22_MCp18->Draw();
        cPtDist_MCp22_MCp18[ebin]->Print(Form("figures_MH/PtDists/eta_%1.1f_%1.1f/MCp22_MCp18_%d_%d.png",etabins[ebin],etabins[ebin+1],cminCENT[setcent],cmaxCENT[setcent]),"png");
        cPtDist_MCp22_MCp18[ebin]->Close();
    }


    TCanvas * crat_MCp22_MCp18[netabins];
    TH1D * rat_MCp22_MCp18[netabins];
    for (int ebin = 0; ebin<netabins; ebin++) {
        crat_MCp22_MCp18[ebin] = new TCanvas(Form("crat_MCp22_MCp18_%d",ebin),"",800,300);
        TPad * padrat_MCp22_MCp18 = (TPad *) crat_MCp22_MCp18[ebin]->cd();
        padrat_MCp22_MCp18->SetGrid();
        TH1D * hrat_MCp22_MCp18 = (TH1D *) hPtDummy->Clone(Form("hrat_MCp22_MCp18_%d",ebin));
        hrat_MCp22_MCp18->GetYaxis()->SetRangeUser(-1.0, 3.0);
        hrat_MCp22_MCp18->GetYaxis()->SetNdivisions(507);
        hrat_MCp22_MCp18->SetXTitle("p_{T} (GeV/c)");
        hrat_MCp22_MCp18->SetYTitle("N1MCp18SUB3 / N1MCp22SUB3");
        hrat_MCp22_MCp18->GetXaxis()->SetTitleSize(0.08);
        hrat_MCp22_MCp18->GetXaxis()->SetTitleOffset(0.90);
        hrat_MCp22_MCp18->GetXaxis()->SetLabelSize(0.07);
        hrat_MCp22_MCp18->GetYaxis()->SetTitleSize(0.07);
        hrat_MCp22_MCp18->GetYaxis()->SetTitleOffset(0.80);
        hrat_MCp22_MCp18->GetYaxis()->SetLabelSize(0.08);
        hrat_MCp22_MCp18->Draw();
        rat_MCp22_MCp18[ebin] = (TH1D *) gN1MCp18SUB3_pt[ebin][setcent]->Clone(Form("rat_MCp22_MCp18_%d",ebin));
        rat_MCp22_MCp18[ebin]->Divide(gN1MCp22SUB3_pt[ebin][setcent]);
        rat_MCp22_MCp18[ebin]->SetMarkerColor(kBlack);
        rat_MCp22_MCp18[ebin]->SetLineColor(kBlack);
        rat_MCp22_MCp18[ebin]->SetMarkerStyle(20);
        rat_MCp22_MCp18[ebin]->SetMarkerSize(1.2);
        for (int i = 1; i<=rat_MCp22_MCp18[ebin]->GetNbinsX(); i++) {
            double xin = gN1MCp18SUB3_pt[ebin][setcent]->GetBinContent(i);
            double yin = gN1MCp22SUB3_pt[ebin][setcent]->GetBinContent(i);
            double delxin = gN1MCp18SUB3_pt[ebin][setcent]->GetBinError(i);
            double delyin = gN1MCp22SUB3_pt[ebin][setcent]->GetBinError(i);
            double raterr = ErrorCalc( xin, yin, delxin, delyin );
            rat_MCp22_MCp18[ebin]->SetBinError(i, raterr);
        }
        rat_MCp22_MCp18[ebin]->Draw("same");
        TLine * lnrat_MCp22_MCp18 = new TLine(0.0, 1.0, 12.0, 1.0);
        lnrat_MCp22_MCp18->SetLineWidth(2);
        lnrat_MCp22_MCp18->Draw();
        crat_MCp22_MCp18[ebin]->Print(Form("figures_MH/PtDists/eta_%1.1f_%1.1f/ratio_MCp22_MCp18_%d_%d.png",etabins[ebin],etabins[ebin+1],cminCENT[setcent],cmaxCENT[setcent]),"png");
        crat_MCp22_MCp18[ebin]->Close();
    }


    TCanvas * cPtDist_MCp18_MCp14[netabins];
    for (int ebin = 0; ebin<netabins; ebin++) {
        cPtDist_MCp18_MCp14[ebin] = new TCanvas(Form("cPtDist_MCp18_MCp14_%d",ebin),"",650,600);
        TPad * padPtDist_MCp18_MCp14 = (TPad *) cPtDist_MCp18_MCp14[ebin]->cd();
        padPtDist_MCp18_MCp14->SetGrid();
        TH1D * hPtDist_MCp18_MCp14 = (TH1D *) hPtDummy->Clone(Form("hPtDist_MCp18_MCp14_%d",ebin));
        hPtDist_MCp18_MCp14->Draw();
        gN1MCp18SUB3_pt[ebin][setcent]->SetMarkerColor(kBlue);
        gN1MCp18SUB3_pt[ebin][setcent]->SetLineColor(kBlue);
        gN1MCp18SUB3_pt[ebin][setcent]->SetMarkerStyle(21);
        gN1MCp18SUB3_pt[ebin][setcent]->SetMarkerSize(1.2);
        gN1MCp18SUB3_pt[ebin][setcent]->Draw("same");
        gN1MCp14SUB3_pt[ebin][setcent]->SetMarkerColor(kGreen+2);
        gN1MCp14SUB3_pt[ebin][setcent]->SetLineColor(kGreen+2);
        gN1MCp14SUB3_pt[ebin][setcent]->SetMarkerStyle(20);
        gN1MCp14SUB3_pt[ebin][setcent]->SetMarkerSize(1.3);
        gN1MCp14SUB3_pt[ebin][setcent]->Draw("same");
        TLegend * legPtDist_MCp18_MCp14 = new TLegend(0.59, 0.18, 0.87, 0.28);
        SetLegend(legPtDist_MCp18_MCp14, 18);
        legPtDist_MCp18_MCp14->AddEntry(gN1MCp18SUB3_pt[ebin][setcent],"#Psi_{1}^{trk}: 1.6 < #eta < 2.0","p");
        legPtDist_MCp18_MCp14->AddEntry(gN1MCp14SUB3_pt[ebin][setcent],"#Psi_{1}^{trk}: 1.2 < #eta < 1.6","p");
        legPtDist_MCp18_MCp14->Draw();
        TPaveText * txPtDist_MCp18_MCp14 = new TPaveText(0.59, 0.30, 0.88, 0.40, "NDC");
        SetTPaveTxt(txPtDist_MCp18_MCp14, 18);
        txPtDist_MCp18_MCp14->AddText("PbPb #sqrt{s_{NN}}=5.02 TeV");
        txPtDist_MCp18_MCp14->AddText(Form("%1.1f < #eta < %1.1f (%d - %d%%)",etabins[ebin],etabins[ebin+1],cminCENT[setcent],cmaxCENT[setcent]));
        txPtDist_MCp18_MCp14->Draw();
        cPtDist_MCp18_MCp14[ebin]->Print(Form("figures_MH/PtDists/eta_%1.1f_%1.1f/MCp18_MCp14_%d_%d.png",etabins[ebin],etabins[ebin+1],cminCENT[setcent],cmaxCENT[setcent]),"png");
        cPtDist_MCp18_MCp14[ebin]->Close();
    }


    TCanvas * crat_MCp18_MCp14[netabins];
    TH1D * rat_MCp18_MCp14[netabins];
    for (int ebin = 0; ebin<netabins; ebin++) {
        crat_MCp18_MCp14[ebin] = new TCanvas(Form("crat_MCp18_MCp14_%d",ebin),"",800,300);
        TPad * padrat_MCp18_MCp14 = (TPad *) crat_MCp18_MCp14[ebin]->cd();
        padrat_MCp18_MCp14->SetGrid();
        TH1D * hrat_MCp18_MCp14 = (TH1D *) hPtDummy->Clone(Form("hrat_MCp18_MCp14_%d",ebin));
        hrat_MCp18_MCp14->GetYaxis()->SetRangeUser(0.0, 2.0);
        hrat_MCp18_MCp14->GetYaxis()->SetNdivisions(507);
        hrat_MCp18_MCp14->SetXTitle("p_{T} (GeV/c)");
        hrat_MCp18_MCp14->SetYTitle("N1MCp14SUB3 / N1MCp18SUB3");
        hrat_MCp18_MCp14->GetXaxis()->SetTitleSize(0.08);
        hrat_MCp18_MCp14->GetXaxis()->SetTitleOffset(0.90);
        hrat_MCp18_MCp14->GetXaxis()->SetLabelSize(0.07);
        hrat_MCp18_MCp14->GetYaxis()->SetTitleSize(0.07);
        hrat_MCp18_MCp14->GetYaxis()->SetTitleOffset(0.80);
        hrat_MCp18_MCp14->GetYaxis()->SetLabelSize(0.08);
        hrat_MCp18_MCp14->Draw();
        rat_MCp18_MCp14[ebin] = (TH1D *) gN1MCp14SUB3_pt[ebin][setcent]->Clone(Form("rat_MCp18_MCp14_%d",ebin));
        rat_MCp18_MCp14[ebin]->Divide(gN1MCp18SUB3_pt[ebin][setcent]);
        rat_MCp18_MCp14[ebin]->SetMarkerColor(kBlack);
        rat_MCp18_MCp14[ebin]->SetLineColor(kBlack);
        rat_MCp18_MCp14[ebin]->SetMarkerStyle(20);
        rat_MCp18_MCp14[ebin]->SetMarkerSize(1.2);
        for (int i = 1; i<=rat_MCp18_MCp14[ebin]->GetNbinsX(); i++) {
            double xin = gN1MCp14SUB3_pt[ebin][setcent]->GetBinContent(i);
            double yin = gN1MCp18SUB3_pt[ebin][setcent]->GetBinContent(i);
            double delxin = gN1MCp14SUB3_pt[ebin][setcent]->GetBinError(i);
            double delyin = gN1MCp18SUB3_pt[ebin][setcent]->GetBinError(i);
            double raterr = ErrorCalc( xin, yin, delxin, delyin );
            rat_MCp18_MCp14[ebin]->SetBinError(i, raterr);
        }
        rat_MCp18_MCp14[ebin]->Draw("same");
        TLine * lnrat_MCp18_MCp14 = new TLine(0.0, 1.0, 12.0, 1.0);
        lnrat_MCp18_MCp14->SetLineWidth(2);
        lnrat_MCp18_MCp14->Draw();
        crat_MCp18_MCp14[ebin]->Print(Form("figures_MH/PtDists/eta_%1.1f_%1.1f/ratio_MCp18_MCp14_%d_%d.png",etabins[ebin],etabins[ebin+1],cminCENT[setcent],cmaxCENT[setcent]),"png");
        crat_MCp18_MCp14[ebin]->Close();
    }



}
