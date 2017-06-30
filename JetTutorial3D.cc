/////////////////////////////////////////////////////////////////////
///                  3D JetTutorial Using Vector                  ///
///                  Zachary Montague 27.06.17                    ///
///                                                               ///
///  Models a three dimensional parton shower and displays plots  ///
///  of the cumulative y-momentum vs. branch level of the shower. ///
///  Overcomes limitations imposed by map by using a vector.      ///
/////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <iterator>
#include <tuple>
#include <ctime>

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TMultiGraph.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <TString.h>
#include <TView.h>
#include <TLorentzVector.h>

//  Each parton has some energy, angle at which
//  it is deflected, and index of its parent parton.
struct parton {
  TLorentzVector* lv = new TLorentzVector();
  double theta;
  double eta;
  int parentId;
};

const double PI  = 3.141592653589793238463;

//  Creates an emitted parton from the parent parton.
//  The emitted parton's properties are extracted
//  from a pseudorandom number generator.
parton emittedParton(parton parent, int parentId, TRandom3* r, TH1F* h_e, TH1F* h_t, TH1F* h_et) {
  parton emit;
  //  Generate random numbers for the fraction of energy from
  //  the parent parton that the emitted parton has and the angle
  //  at which it shoots off. Because the radtion coming off
  //  the splittings is soft, the energy fraction follows a form
  //  of 1 / (1 + z), where there is an additional 1 in the denominator
  //  to prevent an undefined fraction of energy. Similarly,
  //  the splittings are collinear and thus also follow a 1 / (1 + z) form.
  double energy_fraction = 1 / (1 + r -> Uniform(100));
  //  Polar angle.
  double theta = PI  / 2 / (1 + r -> Uniform(100));
  //  Azimuthal angle.
  double eta = PI * r -> Uniform(1);

  //double energy_fraction = r -> Rndm();
  //double angle = r -> Rndm();
  h_e -> Fill(energy_fraction);
  h_t -> Fill(theta);
  h_et -> Fill(eta);

  //  Emitted parton properties.
  emit.lv.SetE(parent.energy * energy_fraction);
  //  Note that the angle is completely independent.
  emit.theta = theta;
  emit.eta = eta;
  emit.lv.SetPx(emit.energy * sin(theta) * cos(eta));
  emit.lv.SetPy(emit.energy * sin(theta) * sin(eta));
  emit.lv.SetPz(emit.energy * cos(theta));
  emit.parentId = parentId;

  return emit;
}

//  Creates a recoiled parton from the parent-emitted parton pair.
//  The recoiled parton's properties are constrained by energy-
//  momentum conservation.
parton recoiledParton(parton parent, parton emit, int parentId) {
  parton recoil;
  
  //  Set recoiled parton properties.
  recoil.lv.SetE(parent.energy - emit.energy);
  recoil.lv.SetPx(parent.px - emit.px);
  recoil.lv.SetPy(parent.py - emit.py);
  recoil.lv.SetPz(parent.pz - emit.pz);
  //  Note the angle is dependent on the child's angle through
  //  momentum conservation.
  std::cout << recoil.lv -> Eta() << std::endl;
  recoil.eta = atan(recoil.lv.Py() / recoil.lv.Px());
  recoil.theta = acos(recoil.lv.Pz() / sqrt(recoil.lv.Px() * recoil.lv.Px() + recoil.lv.Py() * recoil.lv.Py() + recoil.lv.Pz() * recoil.lv.Pz()));
  recoil.parentId = parentId;
  return recoil;
}

//  Runs the 2D parton shower simulation and stores all the
//  parton information into a map.
std::vector<parton> create3DShower(double energy_thres, TRandom3* r, TH1F* h_e, TH1F* h_t, TH1F* h_et) {
  //  Vector of particles.
  std::vector<parton> shower;
  //  Initialize initial parent (root) parton.
  parton parton0 = {{0.0, 0.0, 1.0, 1.0}, 0.0, 0.0, -1};
  shower.push_back(parton0);

  //  Parton simulation loop.
  for (unsigned int i = 0; i < shower.size(); ++i) {
    if (shower[i].energy > energy_thres) {
      // Create emitted and recoiled parton from parent parton.
      shower.push_back(emittedParton(shower[i], i, r, h_e, h_t, h_et));
      shower.push_back(recoiledParton(shower[i], shower.back(), i));
    }
  }

  return shower;
}

//  Prints out information in the dictionary to ensure energy-momentum conservation.
void printShower(std::vector<parton> shower, double energy_thres) {
  std::cout << "3-D Toy Jet Monte Carlo" << std::endl;
  std::cout << "Minimum energy to decay: " << energy_thres << std::endl;
  printf("%8s %13s %11s %12s %11s %12s \n", "index", "parent index", "energy", "theta", "px", "py");
  printf("%8s %13s %12s %13s %13s %12s \n", "-----", "--------------", "-----------", "-----------", "------------", "------------");
  
  double px, py;

  for (unsigned int i = 0; i < shower.size(); ++i) {
    printf("%6d %10d %18.4e %13.4e %12.4e %12.4e \n", i, shower[i].parentId, shower[i].lv(3), shower[i].theta, shower[i].lv(0), shower[i].lv(1));
  }
}

//  Main function.
void JetTutorial3D() {
  //  Set energy threshold inidicating when particles are incapable
  //  of decaying.
  double energy_thres = 0.01;

  //  Initialize random number generator
  TRandom3* r = new TRandom3();

  //  Histograms to check inverse distributions for energy
  //  and angle values.
  TH1F* h_e = new TH1F("rndm_energy_frac", "Distribution of Emitted Parton Fractional Energy; energy; frequency", 50, 0, 1);
  TH1F* h_t = new TH1F("rndm_theta", "Distribution of Emitted Parton Theta; angle; frequency", 50, 0, 2);
  TH1F* h_et = new TH1F("rndm_eta", "Distribution of Emitted Parton Eta; eta; frequency", 50, 0, 3.15);
  //  Get vector of partons from 2D simulation.
  std::vector<parton> shower = create3DShower(energy_thres, r, h_e, h_t, h_et);
  std::cout << "------------------ 2D Simulation Successful ---------------------" << std::endl;
  
  //  Check parton simulation.
  printShower(shower, energy_thres);
  
  //  Plot random variables
  TCanvas* c1 = new TCanvas();
  h_e -> Draw();
  TCanvas* c2 = new TCanvas();
  h_t -> Draw();
  TCanvas* c3 = new TCanvas();
  h_et -> Draw();

  std::vector<double> b;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  
  TCanvas* c = new TCanvas("c", "c");
  TMultiGraph *mg = new TMultiGraph();
  //  TView* v = gPad -> GetView();
  //v->SetRange(-800.0, -800.0, -520.0, 800.0, 800.0, 520.0); 
  //  Loop through shower starting at the last key
  //  to trace back from leaf to root.
  for (unsigned int i = (shower.size() - 1); i != 0; --i) {
    //  Since the shower is being traced backwards, all
    //  unstable particles will already be accounted for.
    if (shower[i].energy > energy_thres) {
      continue;
    }
    
    b.push_back(0);
    x.push_back(shower[i].lv(0));
    y.push_back(shower[i].lv(1));
    z.push_back(shower[i].lv(2));
    
    int parentId = shower[i].parentId;
    //  Go backwards, tracing the child to parent to child, etc.
    while (parentId != -1) {
      // If there is a parent, add a branch level to previous x entries.
      // If there is a parent, add its py to previous y entries.
      for (unsigned int i = 0; i < x.size(); i++) {
	b[i] = b[i] + 1;
	x[i] = x[i] + shower[parentId].lv(0);
	y[i] = y[i] + shower[parentId].lv(1);
	z[i] = z[i] + shower[parentId].lv(2);
      }
      b.push_back(0);
      x.push_back(shower[parentId].lv(0));
      y.push_back(shower[parentId].lv(1));
      z.push_back(shower[parentId].lv(2));
      
      parentId = shower[parentId].parentId;
    }
    
    //  Initialize arrays for coordinates for plotting.
    double b_arr[b.size()];
    double x_arr[x.size()];
    double y_arr[y.size()];
    double z_arr[z.size()];
    
    //  Convert vector to array.
    std::copy(b.begin(), b.end(), b_arr);
    std::copy(x.begin(), x.end(), x_arr);
    std::copy(y.begin(), y.end(), y_arr);
    std::copy(z.begin(), z.end(), z_arr);
    
    TGraph* g = new TGraph(b.size(), b_arr, y_arr);
    //    g -> Draw();
    mg -> Add(g);
    
    //  Clear vector of coordinates to save tracing back
    //  of the next stable parton.
    b.clear();
    x.clear();
    y.clear();
    z.clear();
    gPad -> Modified();
    gPad -> Update();
  }
  mg -> Draw("APL*");
  mg -> SetTitle("2D Parton Shower Model Y-Momentum Paths with Threshold " + TString(std::to_string(energy_thres)));
  mg -> GetXaxis() -> SetTitle("branch level");
  mg -> GetYaxis() -> SetTitle("y");
  
}
