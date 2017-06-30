// ToyJet Tutorial
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
#include <TMultiGraph.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <TString.h>

//  Each parton has some energy, angle at which
//  it is deflected, and index of its parent parton.
struct parton {
  double energy;
  double angle;
  double px;
  double py;
  int parentId;
};

const double PI  = 3.141592653589793238463;

//  Creates an emitted parton from the parent parton.
//  The emitted parton's properties are extracted
//  from a pseudorandom number generator.
parton emittedParton(parton parent, int parentId, TRandom3* r, TH1F* h_e, TH1F* h_a) {
  parton emit;
  //  Generate random numbers for the fraction of energy from
  //  the parent parton that the emitted parton has and the angle
  //  at which it shoots off. Because the radtion coming off
  //  the splittings is soft, the energy fraction follows a form
  //  of 1 / (1 + z), where there is an additional 1 in the denominator
  //  to prevent an undefined fraction of energy. Similarly,
  //  the splittings are collinear and thus also follow a 1 / (1 + z) form.
  double energy_fraction = 1 / (1 + r -> Uniform(100));
  double angle = PI / 2 / (1 + r -> Uniform(100));
  //double energy_fraction = r -> Rndm();
  //double angle = r -> Rndm();
  h_e -> Fill(energy_fraction);
  h_a -> Fill(angle);

  //  Emitted parton properties.
  emit.energy = parent.energy * energy_fraction;
  //  Note that the angle is completely independent.
  emit.angle = angle;
  emit.px = emit.energy * cos(angle);
  emit.py = emit.energy * sin(angle);
  emit.parentId = parentId;

  return emit;
}

//  Creates a recoiled parton from the parent-emitted parton pair.
//  The recoiled parton's properties are constrained by energy-
//  momentum conservation.
parton recoiledParton(parton parent, parton emit, int parentId) {
  parton recoil;
  
  //  Set recoiled parton properties.
  recoil.energy = parent.energy - emit.energy;
  recoil.px = parent.px - emit.px;
  recoil.py = parent.py - emit.py;
  //  Note the angle is dependent on the child's angle through
  //  momentum conservation.
  recoil.angle = atan(recoil.py / recoil.px);
  recoil.parentId = parentId;
  return recoil;
}

//  Runs the 2D parton shower simulation and stores all the
//  parton information into a map.
std::vector<parton> create2DShower(double energy_thres, TRandom3* r, TH1F* h_e, TH1F* h_a) {
  //  Vector of particles.
  std::vector<parton> shower;
  //  Initialize initial parent (root) parton.
  parton parton0 = {1.0, 0.0, 1.0, 0.0, -1};
  shower.push_back(parton0);

  //  Parton simulation loop.
  for (unsigned int i = 0; i < shower.size(); ++i) {
    if (shower[i].energy > energy_thres) {
      // Create emitted and recoiled parton from parent parton.
      shower.push_back(emittedParton(shower[i], i, r, h_e, h_a));
      shower.push_back(recoiledParton(shower[i], shower.back(), i));
    }
  }

  return shower;
}

//  Prints out information in the dictionary to ensure energy-momentum conservation.
void printShower(std::vector<parton> shower, double energy_thres) {
  std::cout << "2-D Toy Jet Monte Carlo" << std::endl;
  std::cout << "Minimum energy to decay: " << energy_thres << std::endl;
  printf("%8s %13s %11s %12s %11s %12s \n", "index", "parent index", "energy", "angle", "px", "py");
  printf("%8s %13s %12s %13s %13s %12s \n", "-----", "--------------", "-----------", "-----------", "------------", "------------");
  
  double px, py;

  for (unsigned int i = 0; i < shower.size(); ++i) {
    printf("%6d %10d %18.4e %13.4e %12.4e %12.4e \n", i, shower[i].parentId, shower[i].energy, shower[i].angle, shower[i].px, shower[i].py);
  }
}

//  Main function.
void JetTutorial2DVec() {
  //  Set energy threshold inidicating when particles are incapable
  //  of decaying.
  double energy_thres = 0.001;

  //  Initialize random number generator
  TRandom3* r = new TRandom3();

  //  Histograms to check inverse distributions for energy
  //  and angle values.
  TH1F* h_e = new TH1F("rndm_energy_frac", "Distribution of Emitted Parton Fractional Energy; energy; frequency", 50, 0, 1);
  TH1F* h_a = new TH1F("rndm_angle", "Distribution of Emitted Parton Angle; angle; frequency", 50, 0, 2);

  //  Get vector of partons from 2D simulation.
  std::vector<parton> shower = create2DShower(energy_thres, r, h_e, h_a);
  std::cout << "------------------ 2D Simulation Successful ---------------------" << std::endl;
  //  Check parton simulation.
  printShower(shower, energy_thres);
  TCanvas* c1 = new TCanvas();
  h_e -> Draw();
  TCanvas* c2 = new TCanvas();
  h_a -> Draw();

  std::vector<double> x;
  std::vector<double> y;
  
  TCanvas* c = new TCanvas("c", "c");
  TMultiGraph *mg = new TMultiGraph();

  //  Loop through shower starting at the last key
  //  to trace back from leaf to root.
  for (unsigned int i = (shower.size() - 1); i != 0; --i) {
    //  Since the shower is being traced backwards, all
    //  unstable particles will already be accounted for.
    if (shower[i].energy > energy_thres) {
      continue;
    }
    
    x.push_back(0);
    y.push_back(shower[i].py);
    
    int parentId = shower[i].parentId;
    //  Go backwards, tracing the child to parent to child, etc.
    while (parentId != -1) {
      // If there is a parent, add a branch level to previous x entries.
      // If there is a parent, add its py to previous y entries.
      for (unsigned int i = 0; i < x.size(); i++) {
	x[i] = x[i] + 1;
	y[i] = y[i] + shower[parentId].py;
      }
      x.push_back(0);
      y.push_back(shower[parentId].py);
      
      parentId = shower[parentId].parentId;
    }
    
    //  Initialize arrays for coordinates for plotting.
    double x_arr[x.size()];
    double y_arr[y.size()];
    
    //  Convert vector to array.
    std::copy(x.begin(), x.end(), x_arr);
    std::copy(y.begin(), y.end(), y_arr);
    
    TGraph* g = new TGraph(x.size(), x_arr, y_arr);
    mg -> Add(g);
    mg -> Draw("APL*");
    
    //  Clear vector of coordinates to save tracing back
    //  of the next stable parton.
    x.clear();
    y.clear();
  }
  mg -> Draw("APL*");
  mg -> SetTitle("2D Parton Shower Model Y-Momentum Paths with Threshold " + TString(std::to_string(energy_thres)));
  mg -> GetXaxis() -> SetTitle("branch level");
  mg -> GetYaxis() -> SetTitle("y");
  gPad -> Modified();
}

