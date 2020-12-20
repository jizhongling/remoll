#ifndef __BEAMLINEDETHISTO_H
#define __BEAMLINEDETHISTO_H

/* NOT MEANT TO BE USED INDEPENDENTLY

   INPUT: output file and  detector ID
*/
#include "anaConst.h"

class beamLineDetHistos {
  private:
    std::vector<TH2F*> xy[nSpecies][nDmg][nFB];
    std::vector<TH1F*> r[nSpecies][nDmg][nFB];
    std::vector<TH1F*> energy[nSpecies][nFB];
    std::map<int, int> ID2entry;

  public:
    beamLineDetHistos() {}
    ~beamLineDetHistos() {}
    void initHisto(TFile*,const int,const string,const string,const double,const int);
    void fillHisto(int detID, int sp, double rdDmg[3], double pz,
		   double xx, double yy, double kinE, const int subDet);
    void writeOutput(TFile *fout, int detID, double scaleFactor,const int subDet);
};

void beamLineDetHistos::initHisto(TFile *fout, const int detID, const string detNm, const string postfix = "", const double range=2000, const int subDet=0){

  fout->cd();
  if(!fout->GetDirectory(Form("det%d",detID)))
    fout->mkdir(Form("det%d",detID), detNm.c_str());
  fout->cd(Form("det%d",detID));

  ID2entry.insert(std::pair<int, int>(detID+10000*subDet,ID2entry.size()));
  for(int k=0;k<nFB;k++)
    for(int i=0;i<nSpecies;i++){
      energy[i][k].push_back(new TH1F(Form("d%d_%s_energy_%s_%s", detID, postfix.c_str(), fbH[k].c_str(), spH[i].c_str()),
				      Form("energy distribution %s %s", fbH[k].c_str(), spH[i].c_str()),
				      121, -8, 4.1));
      niceLogXBins(energy[i][k][ID2entry[detID+10000*subDet]]);
      
      for(int j=0;j<nDmg;j++){
	xy[i][j][k].push_back(new TH2F(Form("d%d_%s_xy_%s_%s_Dmg%d",detID, postfix.c_str(), spH[i].c_str(),fbH[k].c_str(),j),
				       Form("%s for %s %s;x[mm];y[mm]",dmgTit[j].c_str(),fbH[k].c_str(),spTit[i].c_str()),
				       1000, -range, range,
				       1000, -range, range));
	r[i][j][k].push_back(new TH1F(Form("d%d_%s_r_%s_%s_Dmg%d",detID, postfix.c_str(), spH[i].c_str(),fbH[k].c_str(),j),
				      Form("%s for %s %s;r[mm];",dmgTit[j].c_str(),fbH[k].c_str(),spTit[i].c_str()),
				      1000, 0, range));
      }
    }
}

void beamLineDetHistos::fillHisto(int detID, int sp, double rdDmg[3], double pz,
				  double xx, double yy, double kinE, const int subDet=0) {
  if (ID2entry.find(detID+10000*subDet) == ID2entry.end()) 
    return;

  int det = ID2entry[detID+10000*subDet];
  double rr = sqrt(xx*xx+yy*yy);

  energy[sp][0][det]->Fill(kinE);
  if(pz<0)
    energy[sp][2][det]->Fill(kinE);
  else
    energy[sp][1][det]->Fill(kinE);

  for(int kk=0;kk<nDmg;kk++){
    xy[sp][kk][0][det]->Fill(xx,yy,rdDmg[kk]);
    if(pz<0)
      xy[sp][kk][2][det]->Fill(xx,yy,rdDmg[kk]);
    else
      xy[sp][kk][1][det]->Fill(xx,yy,rdDmg[kk]);

    r[sp][kk][0][det]->Fill(rr,rdDmg[kk]);
    if(pz<0)
      r[sp][kk][2][det]->Fill(rr,rdDmg[kk]);
    else
      r[sp][kk][1][det]->Fill(rr,rdDmg[kk]);
  }
}

void beamLineDetHistos::writeOutput(TFile *fout, int detID, double scaleFactor,const int subDet=0) {
  if (ID2entry.find(detID+10000*subDet) == ID2entry.end())
    return;

  int det = ID2entry[detID+10000*subDet];
  fout->cd();
  fout->cd(Form("det%d",detID));
  for(int k=0;k<nFB;k++)
    for(int i=0;i<nSpecies;i++){
      energy[i][k][det]->Scale(scaleFactor);
      energy[i][k][det]->Write();
      
      for(int j=0;j<nDmg;j++){
        xy[i][j][k][det]->Scale(scaleFactor);
        xy[i][j][k][det]->Write();
        r[i][j][k][det]->Scale(scaleFactor);
        r[i][j][k][det]->Write();
      }
    }
}
#endif //__BEAMLINEDETHISTO_H
