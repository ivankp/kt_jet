#include <iostream>
#include <list>

#include <TLorentzVector.h>

// #define __cluster_debug
#define __cluster_time
#include "cluster.hh"

using namespace std;

istream& operator>>(istream &in, list<TLorentzVector>& p) {
  static double px, py, pz, E;
  in >> px >> py >> pz >> E;
  if (in) p.push_back(TLorentzVector(px,py,pz,E));
  return in;
}

namespace clustering {
  template<> inline double __px (const TLorentzVector& p) { return p.Px(); }
  template<> inline double __py (const TLorentzVector& p) { return p.Py(); }
  template<> inline double __pz (const TLorentzVector& p) { return p.Pz(); }
  template<> inline double __E  (const TLorentzVector& p) { return p.E (); }
}

int main()
{
  list<TLorentzVector> particles;
  while ( cin >> particles );

  list<TLorentzVector> jets =
    clustering::cluster<clustering::antikt_alg>(particles,0.6);

  for (list<TLorentzVector>::iterator it=jets.begin(),
       end=jets.end(); it!=end; ++it) cout << it->Pt() << endl;

  return 0;
}
