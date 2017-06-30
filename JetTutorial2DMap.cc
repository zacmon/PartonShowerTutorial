/////////////////////////////////////////////////////////////////////
///                  2D JetTutorial Using Map                     ///
///                  Zachary Montague 27.06.17                    ///
///                                                               ///
///  Models a two dimensional parton shower and displays a graph  ///
///  of the cumulative y-momentum vs. branch level of the shower. /// 
///  Limited use since map key function exceeds in boundary.      ///
/////////////////////////////////////////////////////////////////////


#include <iostream>
#include <map>
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
//  it is deflected, x-momentum, and y-momentum.
struct parton {
  double energy;
  double angle;
  double px;
  double py;
};

typedef std::map<long long int, parton> Dict;
typedef Dict::const_iterator It;
typedef Dict::const_reverse_iterator Rit;

const double PI  = 3.141592653589793238463;

//  Creates a child parton from the parent parton.
//  The child parton's properties are extracted
//  from a random number generator.
parton childParton(parton parent, TRandom3* r, TH1F* h_e, TH1F* h_a) {
  parton child;

  //  Generate random numbers for the fraction of energy from
  //  the parent parton that the child parton has and the angle
  //  at which it shoots off. Because the radtion coming off from
  //  the splittings is soft, the energy fraction follows a form
  //  of 1 / (1 + z), where there is an additional 1 in the denominator
  //  to prevent an undefined fraction of energy. Similarly,
  //  the splittings are collinear and thus also follow a 1 / (1 + z) form.
  
  double energy_fraction = 1 / (1 + r -> Uniform(10));
  double angle = r -> Rndm() * 0.5 + 0.25; //1 / (1 + r -> Uniform(30));

  //double energy_fraction = r -> Rndm();
  //double angle = r -> Rndm();
  h_e -> Fill(energy_fraction);
  h_a -> Fill(angle);

  //  Child parton properties.
  child.energy = parent.energy * energy_fraction;
  //  Note that the angle is completely independent.
  child.angle = angle;
  child.px = child.energy * cos(angle);
  child.py = child.energy * sin(angle);
  return child;
}

//  Creates a sibling parton from the parent-child parton pair.
//  The sibling parton's properties are constrained by energy-
//  momentum conservation.
parton siblingParton(parton parent, parton child) {
  parton sibling;
  
  //  Set sibling parton properties.
  sibling.energy = parent.energy - child.energy;
  sibling.px = parent.px - child.px;
  sibling.py = parent.py - child.py;
  //  Note the angle is dependent on the child's angle through
  //  momentum conservation.
  sibling.angle = atan(sibling.py / sibling.px);
  
  return sibling;
}

//  Runs the 2D parton shower simulation and stores all the
//  parton information into a map.
Dict create2DShower(double energy_thres, TRandom3* r, TH1F* h_e, TH1F* h_a) {
  //  Dictionary of particles.
  Dict shower;
  
  //  Initialize originating parton.
  shower[0] = {1.0, 0.0, 1.0, 0.0};
  
  //  Parton simulation loop.
  for (It it = shower.begin(); it!= shower.end(); ++it) {
    if (it -> second.energy > energy_thres) {
      //  Create a one-to-one correspondence between
      //  parent and children.
      long long int child_num = (it -> first) * 2 + 1;
      long long int sibling_num = (it -> first) * 2 + 2;
      shower[child_num] = childParton(shower[it -> first], r, h_e, h_a);
      shower[sibling_num] = siblingParton(shower[it -> first], shower[child_num]);
      std::cout << "Energy conserved: " << (it -> second.energy - (shower[child_num].energy + shower[sibling_num].energy) < 0.001); 
      std::cout << "\t px conserved: " << (it -> second.px - (shower[child_num].px + shower[sibling_num].px) < 0.001);
      std::cout << "\t py conserved: " << (it -> second.py - (shower[child_num].py + shower[sibling_num].py) < 0.001) << std::endl;
    }
  }

  return shower;
}

//  Prints out information in the dictionary to ensure energy-momentum conservation.
void printDictionary(Dict dict, double energy_thres) {
  std::cout << "2-D Toy Jet Monte Carlo" << std::endl;
  std::cout << "Minimum energy to decay: " << energy_thres << std::endl;
  printf("%8s %9s %13s %11s %13s \n", "index", "energy", "angle", "px", "py");
  printf("%8s %12s %13s %13s %12s \n", "-----", "-----------", "-----------", "------------", "------------");
  
  for (It it(dict.begin()); it != dict.end(); ++it) {
    printf("%6lld %13.4e %13.4e %13.4e %13.4e \n", it -> first, it -> second.energy, it -> second.angle, it -> second.px, it -> second.py);
  }
}

//  Main function.
void JetTutorial2DMap() {
  //  Set energy threshold inidicating when particles are incapable
  //  of decaying.
  double energy_thres = 0.001;

  //  Initialize random number generator
  TRandom3* r = new TRandom3();

  //  Histograms to check inverse distributions for energy
  //  and angle values.
  TH1F* h_e = new TH1F("rndm_energy_frac", "Distribution of Emitted Parton Fractional Energy; energy; frequency", 50, 0, 1);
  TH1F* h_a = new TH1F("rndm_angle", "Distribution of Emitted Parton Angle; angle; frequency", 50, 0, 2);

  //  Get dictionary of partons from 2D simulation.
  Dict shower = create2DShower(energy_thres, r, h_e, h_a);

  //  Check parton simulation.
  printDictionary(shower, energy_thres);
  TCanvas* c1 = new TCanvas();
  h_e -> Draw();
  TCanvas* c2 = new TCanvas();
  h_a -> Draw();

  std::vector<double> x;
  std::vector<double> y;
  
  TCanvas* c = new TCanvas("c", "c");
  TMultiGraph *mg = new TMultiGraph();

  //  Loop through shower starting at the last key.
  for (Rit rit(shower.rbegin()); rit != shower.rend(); ++rit) {
    long long int key = rit -> first;

    // Plotting variables;
    double var_x = 0;
    double var_y = 0;

    //  Plotting variable.
    var_x = ceil(log2(key));
    if (key == 1 || key == 0) {
      var_x = 1;
    }
    var_y = shower[key].py;
    //  Since the shower is being traced backwards, all 
    //  unstable particles will already be accounted for.
    if (shower[key].energy > energy_thres) {
      continue;
    }
    
    //  Save external parton's coordinates.
    x.push_back(var_x);
    y.push_back(var_y);
    
    //  Go backwards, tracing the child to parent to child, etc.
    while (key != 0) {
      if (key % 2 == 0) {
	key = (key - 2) / 2;
      }
      else { 
	key = (key - 1) / 2;
      }
      
      if (key > 1) {
	var_x = (int) log2(key);
      }
      else if (key == 1) {
	var_x = 1;
      }
      else if (key == 0) {
	var_x = 0;
      }
      var_y = shower[key].py;
      
    
      //  Save parents at the beginning for easier processing.
      x.insert(x.begin(), var_x);
      y.insert(y.begin(), var_y);

      //  Add the parent's variable to the child's so it is plotting
      //  going forward.
      for (unsigned int i = 1; i < x.size(); i++) {
	//	x[i] = x[i] + x[0];
	y[i] = y[i] + y[0];
      }
    }
  
    /* for (unsigned int i = 0; i < x.size(); i++) {
      std::cout << "x: " << x[i] << "    y: " << y[i] << std::endl;
      }*/
    
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
