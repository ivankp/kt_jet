#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include "fjcore.hh"
#include "Timer.h"

// #define __cluster_debug
#include "cluster.hh"

using namespace std;

istream& operator>>(istream &in, vector<fjcore::PseudoJet>& p) {
  static double px, py, pz, E;
  in >> px >> py >> pz >> E;
  if (in) p.push_back(fjcore::PseudoJet(px,py,pz,E));
  return in;
}

typedef vector<fjcore::PseudoJet>::const_iterator iter_t;

int main()
{
  const double R = 0.6;
  
  // Read input
  // ****************************************************************

  vector<fjcore::PseudoJet> particles;
  while ( cin >> particles );
  
  // pgJet
  // ****************************************************************

  Timer pg_tm;
  pg_tm.start();
  
  const vector<fjcore::PseudoJet> pg_jets = fjcore::sorted_by_pt(
    clustering::cluster<clustering::antikt_alg>(particles, R)
  );
    
  pg_tm.stop();

  // FastJet
  // ****************************************************************
  
  Timer fj_tm;
  fj_tm.start();
  
  fjcore::JetDefinition jet_def(fjcore::antikt_algorithm, R);
  fjcore::ClusterSequence seq(particles, jet_def);
  const vector<fjcore::PseudoJet> fj_jets = fjcore::sorted_by_pt( seq.inclusive_jets() );
  
  fj_tm.stop();
  
  // Compare
  // ****************************************************************

  const size_t nfj = fj_jets.size(), npg = pg_jets.size();
  bool ok = (nfj==npg);
  
  if (ok) {
    for (size_t i=0;i<nfj;++i)
      if ( fj_jets[i].pt() != pg_jets[i].pt() ) {
        ok = false;
        break;
      }
  }
  
  if (!ok) {
    cout << fixed << setprecision(8);
    for (size_t i=0,n=max(nfj,npg);i<n;++i) {
      cout << setw(12);
      if (i<nfj) cout << fj_jets[i].pt();
      else cout << ' ';
      cout << setw(12);
      if (i<npg) cout << pg_jets[i].pt();
      else cout << ' ';
      cout << endl;
    }
    cout << endl;
  } else {
    cout << "Complete agreement" << endl;
  }

  cout << "pg: " << pg_tm.duration() << " ms" << endl;
  cout << "fj: " << fj_tm.duration() << " ms" << endl;

  return 0;
}
