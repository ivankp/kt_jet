#include <iostream>
#include <list>

#include <TLorentzVector.h>

#include "cluster.hh"

using namespace std;

istream& operator>>(istream &in, list<TLorentzVector>& p) {
  static double px, py, pz, E;
  in >> px >> py >> pz >> E;
  if (in) p.push_back(TLorentzVector(px,py,pz,E));
  return in;
}

namespace clustering {
  template<> inline double __pt (const TLorentzVector& p) { return p.Pt(); }
  template<> inline double __rap(const TLorentzVector& p) { return p.Rapidity(); }
  template<> inline double __phi(const TLorentzVector& p) { return p.Phi(); }
}

int main()
{
  list<TLorentzVector> particles;
  while ( cin >> particles );
  
  list<TLorentzVector> jets = clustering::cluster<clustering::antikt_alg>(particles,0.6);
  
  for (list<TLorentzVector>::iterator it=jets.begin(), end=jets.end(); it!=end; ++it)
    cout << it->Pt();

  return 0;
}
