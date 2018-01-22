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
static int NANALS = 2;
TString ANAL[] = {"N1SUB2", "N1SUB3"};

TH1D * gN1SUB2_tight[ncentbins];
TH1D * gN1SUB3_tight[ncentbins];

TH1D * gN1SUB2_tight2[ncentbins];
TH1D * gN1SUB3_tight2[ncentbins];

TH1D * gN1SUB2_narrow[ncentbins];
TH1D * gN1SUB3_narrow[ncentbins];

TH1D * gN1SUB2_wide[ncentbins];
TH1D * gN1SUB3_wide[ncentbins];

void plotEtaDistOdd()
{

    for (int cbin = 0; cbin<ncentbins; cbin++) {
        gN1SUB2_tight[cbin] = new TH1D(Form("gN1SUB2_tight_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1SUB3_tight[cbin] = new TH1D(Form("gN1SUB3_tight_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);

        gN1SUB2_tight2[cbin] = new TH1D(Form("gN1SUB2_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1SUB3_tight2[cbin] = new TH1D(Form("gN1SUB3_tight2_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);

        gN1SUB2_narrow[cbin] = new TH1D(Form("gN1SUB2_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1SUB3_narrow[cbin] = new TH1D(Form("gN1SUB3_narrow_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);

        gN1SUB2_wide[cbin] = new TH1D(Form("gN1SUB2_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
        gN1SUB3_wide[cbin] = new TH1D(Form("gN1SUB3_wide_%d_%d",cminCENT[cbin],cmaxCENT[cbin]), "", netabins, etabins);
    }

    for (int i = 0; i<NANALS; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            TString tag = Form("figures_MH/%s/EtaDistributions/data/EtaInt_%s_%d_%d.dat",ANAL[i].Data(),ANAL[i].Data(),cminCENT[cbin],cmaxCENT[cbin]);
            // cout<<tag.Data()<<endl;
            ifstream fin(tag.Data());

            double eta, v1, v1e, v1A, v1Ae, v1B, v1Be;
            int neta = 0;
            while (fin >> eta >> v1 >> v1e >> v1A >> v1Ae >> v1B >> v1Be) {
                // if (i == 0 && cbin == 0) cout<<"neta = "<<neta<<"\tdeltaEta = "<<deltaEta[neta]<<"\tv1A = "<<v1A<<"\tv1Ae = "<<v1Ae<<endl;

                if (i == 0) {
                    gN1SUB2_tight[cbin]->SetBinContent(neta+1, v1A);
                    gN1SUB2_tight[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 1) {
                    gN1SUB3_tight[cbin]->SetBinContent(neta+1, v1A);
                    gN1SUB3_tight[cbin]->SetBinError(neta+1, v1Ae);
                }
                neta++;
            }
        }
    }

    for (int i = 0; i<NANALS; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            TString tag = Form("figures_MH_tight2/%s/EtaDistributions/data/EtaInt_%s_%d_%d.dat",ANAL[i].Data(),ANAL[i].Data(),cminCENT[cbin],cmaxCENT[cbin]);
            // cout<<tag.Data()<<endl;
            ifstream fin(tag.Data());

            double eta, v1, v1e, v1A, v1Ae, v1B, v1Be;
            int neta = 0;
            while (fin >> eta >> v1 >> v1e >> v1A >> v1Ae >> v1B >> v1Be) {
                // if (i == 0 && cbin == 0) cout<<"neta = "<<neta<<"\tdeltaEta = "<<deltaEta[neta]<<"\tv1A = "<<v1A<<"\tv1Ae = "<<v1Ae<<endl;

                if (i == 0) {
                    gN1SUB2_tight2[cbin]->SetBinContent(neta+1, v1A);
                    gN1SUB2_tight2[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 1) {
                    gN1SUB3_tight2[cbin]->SetBinContent(neta+1, v1A);
                    gN1SUB3_tight2[cbin]->SetBinError(neta+1, v1Ae);
                }
                neta++;
            }
        }
    }

    for (int i = 0; i<NANALS; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            TString tag = Form("figures_MH_narrow/%s/EtaDistributions/data/EtaInt_%s_%d_%d.dat",ANAL[i].Data(),ANAL[i].Data(),cminCENT[cbin],cmaxCENT[cbin]);
            // cout<<tag.Data()<<endl;
            ifstream fin(tag.Data());

            double eta, v1, v1e, v1A, v1Ae, v1B, v1Be;
            int neta = 0;
            while (fin >> eta >> v1 >> v1e >> v1A >> v1Ae >> v1B >> v1Be) {
                // if (i == 0 && cbin == 0) cout<<"neta = "<<neta<<"\tdeltaEta = "<<deltaEta[neta]<<"\tv1A = "<<v1A<<"\tv1Ae = "<<v1Ae<<endl;

                if (i == 0) {
                    gN1SUB2_narrow[cbin]->SetBinContent(neta+1, v1A);
                    gN1SUB2_narrow[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 1) {
                    gN1SUB3_narrow[cbin]->SetBinContent(neta+1, v1A);
                    gN1SUB3_narrow[cbin]->SetBinError(neta+1, v1Ae);
                }
                neta++;
            }
        }
    }

    for (int i = 0; i<NANALS; i++) {
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            TString tag = Form("figures_MH_wide/%s/EtaDistributions/data/EtaInt_%s_%d_%d.dat",ANAL[i].Data(),ANAL[i].Data(),cminCENT[cbin],cmaxCENT[cbin]);
            // cout<<tag.Data()<<endl;
            ifstream fin(tag.Data());

            double eta, v1, v1e, v1A, v1Ae, v1B, v1Be;
            int neta = 0;
            while (fin >> eta >> v1 >> v1e >> v1A >> v1Ae >> v1B >> v1Be) {
                // if (i == 0 && cbin == 0) cout<<"neta = "<<neta<<"\tdeltaEta = "<<deltaEta[neta]<<"\tv1A = "<<v1A<<"\tv1Ae = "<<v1Ae<<endl;

                if (i == 0) {
                    gN1SUB2_wide[cbin]->SetBinContent(neta+1, v1A);
                    gN1SUB2_wide[cbin]->SetBinError(neta+1, v1Ae);
                } else if (i == 1) {
                    gN1SUB3_wide[cbin]->SetBinContent(neta+1, v1A);
                    gN1SUB3_wide[cbin]->SetBinError(neta+1, v1Ae);
                }
                neta++;
            }
        }
    }

    for (int cbin = 0; cbin<ncentbins; cbin++) {
        gN1SUB3_tight[cbin]->SetMarkerColor(kBlue);
        gN1SUB3_tight[cbin]->SetLineColor(kBlue);
        gN1SUB3_tight[cbin]->SetMarkerStyle(21);
        gN1SUB3_tight[cbin]->SetMarkerSize(1.2);

        gN1SUB3_tight2[cbin]->SetMarkerColor(kRed);
        gN1SUB3_tight2[cbin]->SetLineColor(kRed);
        gN1SUB3_tight2[cbin]->SetMarkerStyle(25);
        gN1SUB3_tight2[cbin]->SetMarkerSize(1.2);

        gN1SUB3_narrow[cbin]->SetMarkerColor(kGreen+2);
        gN1SUB3_narrow[cbin]->SetLineColor(kGreen+2);
        gN1SUB3_narrow[cbin]->SetMarkerStyle(24);
        gN1SUB3_narrow[cbin]->SetMarkerSize(1.3);

        gN1SUB3_wide[cbin]->SetMarkerColor(kMagenta);
        gN1SUB3_wide[cbin]->SetLineColor(kMagenta);
        gN1SUB3_wide[cbin]->SetMarkerStyle(34);
        gN1SUB3_wide[cbin]->SetMarkerSize(1.7);
    }

    if (!fopen("figures_MH/systematics","r")) system("mkdir figures_MH/systematics");

    TCanvas * cTight_Tight2_cent = new TCanvas("cTight_Tight2_cent","",1100,500);
    cTight_Tight2_cent->Divide(4,2,0,0);
    for (int cbin = 0; cbin<8; cbin++) {
        TPad * padTight_Tight2_cent = (TPad *) cTight_Tight2_cent->cd(cbin+1);
        padTight_Tight2_cent->SetGrid();
        TH1D * hTight_Tight2_cent = new TH1D(Form("hTight_Tight2_cent_%d",cbin), "", 100, -2.5, 2.5);
        hTight_Tight2_cent->GetXaxis()->SetNdivisions(509);
        hTight_Tight2_cent->GetYaxis()->SetNdivisions(509);
        hTight_Tight2_cent->SetXTitle("#eta");
        hTight_Tight2_cent->SetYTitle("v_{1}^{odd}");
        hTight_Tight2_cent->GetYaxis()->SetRangeUser(-0.01, 0.01);
        hTight_Tight2_cent->Draw();
        gN1SUB3_tight[cbin]->Draw("same");
        gN1SUB3_tight2[cbin]->Draw("same");
        TPaveText * Tight_Tight2_centscan = new TPaveText(0.73, 0.83, 0.91, 0.95, "NDC");
        SetTPaveTxt(Tight_Tight2_centscan, 18);
        Tight_Tight2_centscan->AddText(Form("%d-%d%%",cminCENT[cbin],cmaxCENT[cbin]));
        Tight_Tight2_centscan->Draw();
    }
    cTight_Tight2_cent->cd(1);
    TPaveText * txTight_Tight2_cent = new TPaveText(0.23, 0.03, 0.77, 0.28, "NDC");
    SetTPaveTxt(txTight_Tight2_cent, 18);
    txTight_Tight2_cent->AddText("PbPb #sqrt{s_{NN}}=5.02 TeV");
    txTight_Tight2_cent->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txTight_Tight2_cent->Draw();
    cTight_Tight2_cent->cd(5);
    TLegend * legTight_Tight2_cent = new TLegend(0.23, 0.21, 0.64, 0.43);
    SetLegend(legTight_Tight2_cent, 18);
    legTight_Tight2_cent->AddEntry(gN1SUB3_tight[0],"v_{1}^{even} (default)","p");
    legTight_Tight2_cent->AddEntry(gN1SUB3_tight2[0],"v_{1}^{even} (tight2)","p");
    legTight_Tight2_cent->Draw();
    cTight_Tight2_cent->Print("figures_MH/systematics/v1odd_Tight_Tight2_centscan.png","png");


    TCanvas * cratTight_Tight2_cent = new TCanvas("cratTight_Tight2_cent","",1100,500);
    cratTight_Tight2_cent->Divide(4,2,0,0);
    for (int cbin = 0; cbin<8; cbin++) {
        TPad * padratTight_Tight2_cent = (TPad *) cratTight_Tight2_cent->cd(cbin+1);
        padratTight_Tight2_cent->SetGrid();
        TH1D * hratTight_Tight2_cent = new TH1D(Form("hratTight_Tight2_cent_%d",cbin), "", 100, -2.5, 2.5);
        hratTight_Tight2_cent->GetXaxis()->SetNdivisions(509);
        hratTight_Tight2_cent->GetYaxis()->SetNdivisions(509);
        hratTight_Tight2_cent->GetYaxis()->SetRangeUser(0.8, 1.2);
        hratTight_Tight2_cent->SetXTitle("#eta");
        hratTight_Tight2_cent->SetYTitle("v_{1}^{odd} {default} / v_{1}^{odd} {tight2}");
        hratTight_Tight2_cent->Draw();
        TLine * lnratTight_Tight2_cent = new TLine(-2.5, 1.0, 2.5, 1.0);
        lnratTight_Tight2_cent->SetLineWidth(2);
        lnratTight_Tight2_cent->Draw();
        TH1D * gratN1SUB3_tight_tight2 = (TH1D *) gN1SUB3_tight[cbin]->Clone(Form("gratN1SUB3_tight_tight2_%d",cbin));
        gratN1SUB3_tight_tight2->Divide(gN1SUB3_tight2[cbin]);
        gratN1SUB3_tight_tight2->SetMarkerColor(kBlack);
        gratN1SUB3_tight_tight2->SetLineColor(kBlack);
        gratN1SUB3_tight_tight2->SetMarkerStyle(20);
        gratN1SUB3_tight_tight2->SetMarkerSize(1.2);
        for (int i = 1; i<=gratN1SUB3_tight_tight2->GetNbinsX(); i++) {
            double xin = gN1SUB3_tight[cbin]->GetBinContent(i);
            double yin = gN1SUB3_tight2[cbin]->GetBinContent(i);
            double delxin = gN1SUB3_tight[cbin]->GetBinError(i);
            double delyin = gN1SUB3_tight2[cbin]->GetBinError(i);
            double raterr = ErrorCalc( xin, yin, delxin, delyin );
            gratN1SUB3_tight_tight2->SetBinError(i, raterr);
        }
        gratN1SUB3_tight_tight2->Draw("same");
        TPaveText * Tight_Tight2_centscan = new TPaveText(0.73, 0.83, 0.91, 0.95, "NDC");
        SetTPaveTxt(Tight_Tight2_centscan, 18);
        Tight_Tight2_centscan->AddText(Form("%d-%d%%",cminCENT[cbin],cmaxCENT[cbin]));
        Tight_Tight2_centscan->Draw();
    }
    cratTight_Tight2_cent->cd(1);
    TPaveText * txratTight_Tight2_cent = new TPaveText(0.23, 0.03, 0.77, 0.28, "NDC");
    SetTPaveTxt(txratTight_Tight2_cent, 18);
    txratTight_Tight2_cent->AddText("PbPb #sqrt{s_{NN}}=5.02 TeV");
    txratTight_Tight2_cent->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txratTight_Tight2_cent->Draw();
    cratTight_Tight2_cent->cd(5);
    cratTight_Tight2_cent->Print("figures_MH/systematics/v1odd_ratio_Tight_Tight2_centscan.png","png");


    TCanvas * cTight_narrow_cent = new TCanvas("cTight_narrow_cent","",1100,500);
    cTight_narrow_cent->Divide(4,2,0,0);
    for (int cbin = 0; cbin<8; cbin++) {
        TPad * padTight_narrow_cent = (TPad *) cTight_narrow_cent->cd(cbin+1);
        padTight_narrow_cent->SetGrid();
        TH1D * hTight_narrow_cent = new TH1D(Form("hTight_narrow_cent_%d",cbin), "", 100, -2.5, 2.5);
        hTight_narrow_cent->GetXaxis()->SetNdivisions(509);
        hTight_narrow_cent->GetYaxis()->SetNdivisions(509);
        hTight_narrow_cent->SetXTitle("#eta");
        hTight_narrow_cent->SetYTitle("v_{1}^{odd}");
        hTight_narrow_cent->GetYaxis()->SetRangeUser(-0.01, 0.01);
        hTight_narrow_cent->Draw();
        gN1SUB3_tight[cbin]->Draw("same");
        gN1SUB3_narrow[cbin]->Draw("same");
        TPaveText * Tight_narrow_centscan = new TPaveText(0.73, 0.83, 0.91, 0.95, "NDC");
        SetTPaveTxt(Tight_narrow_centscan, 18);
        Tight_narrow_centscan->AddText(Form("%d-%d%%",cminCENT[cbin],cmaxCENT[cbin]));
        Tight_narrow_centscan->Draw();
    }
    cTight_narrow_cent->cd(1);
    TPaveText * txTight_narrow_cent = new TPaveText(0.23, 0.03, 0.77, 0.28, "NDC");
    SetTPaveTxt(txTight_narrow_cent, 18);
    txTight_narrow_cent->AddText("PbPb #sqrt{s_{NN}}=5.02 TeV");
    txTight_narrow_cent->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txTight_narrow_cent->Draw();
    cTight_narrow_cent->cd(5);
    TLegend * legTight_narrow_cent = new TLegend(0.23, 0.21, 0.64, 0.43);
    SetLegend(legTight_narrow_cent, 18);
    legTight_narrow_cent->AddEntry(gN1SUB3_tight[0],"v_{1}^{even} (default)","p");
    legTight_narrow_cent->AddEntry(gN1SUB3_narrow[0],"v_{1}^{even} (narrow)","p");
    legTight_narrow_cent->Draw();
    cTight_narrow_cent->Print("figures_MH/systematics/v1odd_Tight_narrow_centscan.png","png");


    TCanvas * cratTight_narrow_cent = new TCanvas("cratTight_narrow_cent","",1100,500);
    cratTight_narrow_cent->Divide(4,2,0,0);
    for (int cbin = 0; cbin<8; cbin++) {
        TPad * padratTight_narrow_cent = (TPad *) cratTight_narrow_cent->cd(cbin+1);
        padratTight_narrow_cent->SetGrid();
        TH1D * hratTight_narrow_cent = new TH1D(Form("hratTight_narrow_cent_%d",cbin), "", 100, -2.5, 2.5);
        hratTight_narrow_cent->GetXaxis()->SetNdivisions(509);
        hratTight_narrow_cent->GetYaxis()->SetNdivisions(509);
        hratTight_narrow_cent->GetYaxis()->SetRangeUser(0.5, 1.5);
        hratTight_narrow_cent->SetXTitle("#eta");
        hratTight_narrow_cent->SetYTitle("v_{1}^{odd} {default} / v_{1}^{odd} {narrow}");
        hratTight_narrow_cent->Draw();
        TLine * lnratTight_narrow_cent = new TLine(-2.5, 1.0, 2.5, 1.0);
        lnratTight_narrow_cent->SetLineWidth(2);
        lnratTight_narrow_cent->Draw();
        TH1D * gratN1SUB3_tight_narrow = (TH1D *) gN1SUB3_tight[cbin]->Clone(Form("gratN1SUB3_tight_narrow_%d",cbin));
        gratN1SUB3_tight_narrow->Divide(gN1SUB3_narrow[cbin]);
        gratN1SUB3_tight_narrow->SetMarkerColor(kBlack);
        gratN1SUB3_tight_narrow->SetLineColor(kBlack);
        gratN1SUB3_tight_narrow->SetMarkerStyle(20);
        gratN1SUB3_tight_narrow->SetMarkerSize(1.2);
        for (int i = 1; i<=gratN1SUB3_tight_narrow->GetNbinsX(); i++) {
            double xin = gN1SUB3_tight[cbin]->GetBinContent(i);
            double yin = gN1SUB3_narrow[cbin]->GetBinContent(i);
            double delxin = gN1SUB3_tight[cbin]->GetBinError(i);
            double delyin = gN1SUB3_narrow[cbin]->GetBinError(i);
            double raterr = ErrorCalc( xin, yin, delxin, delyin );
            gratN1SUB3_tight_narrow->SetBinError(i, raterr);
        }
        gratN1SUB3_tight_narrow->Draw("same");
        TPaveText * Tight_narrow_centscan = new TPaveText(0.73, 0.83, 0.91, 0.95, "NDC");
        SetTPaveTxt(Tight_narrow_centscan, 18);
        Tight_narrow_centscan->AddText(Form("%d-%d%%",cminCENT[cbin],cmaxCENT[cbin]));
        Tight_narrow_centscan->Draw();
    }
    cratTight_narrow_cent->cd(1);
    TPaveText * txratTight_narrow_cent = new TPaveText(0.23, 0.03, 0.77, 0.28, "NDC");
    SetTPaveTxt(txratTight_narrow_cent, 18);
    txratTight_narrow_cent->AddText("PbPb #sqrt{s_{NN}}=5.02 TeV");
    txratTight_narrow_cent->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txratTight_narrow_cent->Draw();
    cratTight_narrow_cent->cd(5);
    cratTight_narrow_cent->Print("figures_MH/systematics/v1odd_ratio_Tight_narrow_centscan.png","png");


    TCanvas * cTight_wide_cent = new TCanvas("cTight_wide_cent","",1100,500);
    cTight_wide_cent->Divide(4,2,0,0);
    for (int cbin = 0; cbin<8; cbin++) {
        TPad * padTight_wide_cent = (TPad *) cTight_wide_cent->cd(cbin+1);
        padTight_wide_cent->SetGrid();
        TH1D * hTight_wide_cent = new TH1D(Form("hTight_wide_cent_%d",cbin), "", 100, -2.5, 2.5);
        hTight_wide_cent->GetXaxis()->SetNdivisions(509);
        hTight_wide_cent->GetYaxis()->SetNdivisions(509);
        hTight_wide_cent->SetXTitle("#eta");
        hTight_wide_cent->SetYTitle("v_{1}^{odd}");
        hTight_wide_cent->GetYaxis()->SetRangeUser(-0.01, 0.01);
        hTight_wide_cent->Draw();
        gN1SUB3_tight[cbin]->Draw("same");
        gN1SUB3_wide[cbin]->Draw("same");
        TPaveText * Tight_wide_centscan = new TPaveText(0.73, 0.83, 0.91, 0.95, "NDC");
        SetTPaveTxt(Tight_wide_centscan, 18);
        Tight_wide_centscan->AddText(Form("%d-%d%%",cminCENT[cbin],cmaxCENT[cbin]));
        Tight_wide_centscan->Draw();
    }
    cTight_wide_cent->cd(1);
    TPaveText * txTight_wide_cent = new TPaveText(0.23, 0.03, 0.77, 0.28, "NDC");
    SetTPaveTxt(txTight_wide_cent, 18);
    txTight_wide_cent->AddText("PbPb #sqrt{s_{NN}}=5.02 TeV");
    txTight_wide_cent->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txTight_wide_cent->Draw();
    cTight_wide_cent->cd(5);
    TLegend * legTight_wide_cent = new TLegend(0.23, 0.21, 0.64, 0.43);
    SetLegend(legTight_wide_cent, 18);
    legTight_wide_cent->AddEntry(gN1SUB3_tight[0],"v_{1}^{even} (default)","p");
    legTight_wide_cent->AddEntry(gN1SUB3_wide[0],"v_{1}^{even} (wide)","p");
    legTight_wide_cent->Draw();
    cTight_wide_cent->Print("figures_MH/systematics/v1odd_Tight_wide_centscan.png","png");


    TCanvas * cratTight_wide_cent = new TCanvas("cratTight_wide_cent","",1100,500);
    cratTight_wide_cent->Divide(4,2,0,0);
    for (int cbin = 0; cbin<8; cbin++) {
        TPad * padratTight_wide_cent = (TPad *) cratTight_wide_cent->cd(cbin+1);
        padratTight_wide_cent->SetGrid();
        TH1D * hratTight_wide_cent = new TH1D(Form("hratTight_wide_cent_%d",cbin), "", 100, -2.5, 2.5);
        hratTight_wide_cent->GetXaxis()->SetNdivisions(509);
        hratTight_wide_cent->GetYaxis()->SetNdivisions(509);
        hratTight_wide_cent->GetYaxis()->SetRangeUser(0.8, 1.2);
        hratTight_wide_cent->SetXTitle("#eta");
        hratTight_wide_cent->SetYTitle("v_{1}^{odd} {default} / v_{1}^{odd} {wide}");
        hratTight_wide_cent->Draw();
        TLine * lnratTight_wide_cent = new TLine(-2.5, 1.0, 2.5, 1.0);
        lnratTight_wide_cent->SetLineWidth(2);
        lnratTight_wide_cent->Draw();
        TH1D * gratN1SUB3_tight_wide = (TH1D *) gN1SUB3_tight[cbin]->Clone(Form("gratN1SUB3_tight_wide_%d",cbin));
        gratN1SUB3_tight_wide->Divide(gN1SUB3_wide[cbin]);
        gratN1SUB3_tight_wide->SetMarkerColor(kBlack);
        gratN1SUB3_tight_wide->SetLineColor(kBlack);
        gratN1SUB3_tight_wide->SetMarkerStyle(20);
        gratN1SUB3_tight_wide->SetMarkerSize(1.2);
        for (int i = 1; i<=gratN1SUB3_tight_wide->GetNbinsX(); i++) {
            double xin = gN1SUB3_tight[cbin]->GetBinContent(i);
            double yin = gN1SUB3_wide[cbin]->GetBinContent(i);
            double delxin = gN1SUB3_tight[cbin]->GetBinError(i);
            double delyin = gN1SUB3_wide[cbin]->GetBinError(i);
            double raterr = ErrorCalc( xin, yin, delxin, delyin );
            gratN1SUB3_tight_wide->SetBinError(i, raterr);
        }
        gratN1SUB3_tight_wide->Draw("same");
        TPaveText * Tight_wide_centscan = new TPaveText(0.73, 0.83, 0.91, 0.95, "NDC");
        SetTPaveTxt(Tight_wide_centscan, 18);
        Tight_wide_centscan->AddText(Form("%d-%d%%",cminCENT[cbin],cmaxCENT[cbin]));
        Tight_wide_centscan->Draw();
    }
    cratTight_wide_cent->cd(1);
    TPaveText * txratTight_wide_cent = new TPaveText(0.23, 0.03, 0.77, 0.28, "NDC");
    SetTPaveTxt(txratTight_wide_cent, 18);
    txratTight_wide_cent->AddText("PbPb #sqrt{s_{NN}}=5.02 TeV");
    txratTight_wide_cent->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txratTight_wide_cent->Draw();
    cratTight_wide_cent->cd(5);
    cratTight_wide_cent->Print("figures_MH/systematics/v1odd_ratio_Tight_wide_centscan.png","png");


    TCanvas * cv1odd_syst = new TCanvas("cv1odd_syst","",1100,500);
    cv1odd_syst->Divide(4,2,0,0);
    for (int cbin = 0; cbin<8; cbin++) {
        TPad * padv1odd_syst = (TPad *) cv1odd_syst->cd(cbin+1);
        padv1odd_syst->SetGrid();
        TH1D * hv1odd_syst = new TH1D(Form("hv1odd_syst_%d",cbin), "", 100, -2.5, 2.5);
        hv1odd_syst->GetXaxis()->SetNdivisions(509);
        hv1odd_syst->GetYaxis()->SetNdivisions(509);
        hv1odd_syst->SetXTitle("#eta");
        hv1odd_syst->SetYTitle("v_{1}^{odd}");
        hv1odd_syst->GetYaxis()->SetRangeUser(-0.01, 0.01);
        hv1odd_syst->Draw();
        gN1SUB3_tight[cbin]->Draw("same");
        gN1SUB3_wide[cbin]->Draw("same");
        TPaveText * v1odd_systscan = new TPaveText(0.73, 0.83, 0.91, 0.95, "NDC");
        SetTPaveTxt(v1odd_systscan, 18);
        v1odd_systscan->AddText(Form("%d-%d%%",cminCENT[cbin],cmaxCENT[cbin]));
        v1odd_systscan->Draw();
    }
    cv1odd_syst->cd(1);
    TPaveText * txv1odd_syst = new TPaveText(0.23, 0.03, 0.77, 0.28, "NDC");
    SetTPaveTxt(txv1odd_syst, 18);
    txv1odd_syst->AddText("PbPb #sqrt{s_{NN}}=5.02 TeV");
    txv1odd_syst->AddText("0.3 < p_{T} < 3.0 GeV/c");
    txv1odd_syst->Draw();
    cv1odd_syst->cd(5);
    TLegend * legv1odd_syst = new TLegend(0.23, 0.21, 0.64, 0.43);
    SetLegend(legv1odd_syst, 18);
    legv1odd_syst->AddEntry(gN1SUB3_tight[0],"v_{1}^{even} (default)","p");
    legv1odd_syst->AddEntry(gN1SUB3_wide[0],"v_{1}^{even} (wide)","p");
    legv1odd_syst->Draw();
    cv1odd_syst->Print("figures_MH/systematics/v1odd_v1odd_systscan.png","png");

}
