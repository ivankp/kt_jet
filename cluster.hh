#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

#ifdef __cluster_debug
#include <iomanip>
#endif

namespace clustering {

// helpful functions
// ******************************************************************

template<typename T> inline T sq(T x) { return x*x; }

// functions to get variables from particle
// ******************************************************************

template<class P> inline double __px (const P& p) { return p.px(); }
template<class P> inline double __py (const P& p) { return p.py(); }
template<class P> inline double __pz (const P& p) { return p.pz(); }
template<class P> inline double __E  (const P& p) { return p.E (); }

// clustering type enum
// ******************************************************************

enum algorithm { kt_alg, antikt_alg };

// 4-momentum structure
// ******************************************************************

template<algorithm alg>
struct p4 {
  double px, py, pz, E, d, rap, phi;
  int id;

  p4() { }
  p4(double px, double py, double pz, double E)
  : px(px), py(py), pz(pz), E(E),
    d( alg==kt_alg ? pt2() : 1./pt2() ),
    rap(0.5*log((E+pz)/(E-pz))),
    phi(px == 0. && py == 0. ? 0. : std::atan2(py,px)),
    id(num++) { }

  p4<alg> operator+(const p4<alg>& p) const {
    return p4<alg>(px+p.px,py+p.py,pz+p.pz,E+p.E);
  }

  // distance between two particles
  inline double dij(const p4<alg>& p) const {
    double deltaPhi = std::fabs(phi-p.phi);
    if (deltaPhi > M_PI) deltaPhi = 2*M_PI - deltaPhi;
    return std::min(d,p.d) * ( sq(rap-p.rap) + sq(deltaPhi) );
  }

  // comparison operator for sorting by distance in ascending order
  bool operator<(const p4<alg>& p) const { return ( d < p.d ); }

private:
  inline double pt2() const { return px*px + py*py; }
  static int num;
};
template<algorithm alg> int p4<alg>::num = 0;

// clustering function
// ******************************************************************
template<algorithm alg, class InputContainer>
std::vector<typename InputContainer::value_type>
cluster(const InputContainer& pp, double R)
{
  // input particle type
  typedef typename InputContainer::value_type pp_type;
  typedef typename InputContainer::const_iterator pp_iter;

  const size_t n = pp.size();
  size_t n_ok = n;
  p4<alg> particles[n]; // particles
  bool ok[n];           // particle exists
  double dij[n][n];     // cache of pairwise distances

  // read input particles
  for (pp_iter it=pp.begin(), end=pp.end(); it!=end; ++it) {
    static size_t i=0;
    particles[i] = p4<alg>( __px(*it),__py(*it),__pz(*it),__E(*it) );
    ok[i] = true;
    ++i;
  }

  // cache pairwise distances
  for (size_t i=0, ni=n-1; i<ni; ++i)
    for (size_t j=i+1; j<n; ++j)
      dij[i][j] = particles[i].dij( particles[j] );

  #ifdef __cluster_debug
  for (size_t i=0; i<n; ++i)
    std::cout << particles[i].id << ": " << particles[i].d << std::endl;
  std::cout << std::endl;
  #endif

  // output list of jets with constituents
  std::vector<pp_type> jets;
  jets.reserve(n/10+1);

  // perform clustering iterations until no more particles left
  while (n_ok) {
    static size_t i1, i2;

    double dist = std::numeric_limits<double>::max();
    bool merged = false;

    // find smallest single distance
    for (size_t i=0; i<n; ++i) {
      if (!ok[i]) continue;
      double d = particles[i].d*R*R;
      if (d < dist) {
        dist = d;
        i1 = i;
      }
    }

    // find closest pair
    for (size_t i=0, ni=n-1; i<ni; ++i) {
      if (!ok[i]) continue;
      for (size_t j=i+1; j<n; ++j) {
        if (!ok[j]) continue;

        double d = dij[i][j];
        if (d < dist) {
          dist = d;
          i1 = i;
          i2 = j;
          if (!merged) merged = true;
        }
      }
    }

    if (merged) {
      // merge particles
      p4<alg> p( particles[i1] + particles[i2] );

      // print clustering step
      #ifdef __cluster_debug
      std::cout << std::setw(3) << p.id << ": merged "
                << std::setw(3) << particles[i1].id << " & "
                << std::setw(3) << particles[i2].id
                << " | d = " << dist << std::endl;
      #endif

      // "remove" merged particles
      particles[i1] = p;
      ok[i2] = false;

      // cache new pairwise distances
      for (size_t i=0; i<i1; ++i) {
        if (!ok[i]) continue;
        dij[i][i1] = particles[i].dij( particles[i1] );
      }
      for (size_t i=i1+1; i<n; ++i) {
        if (!ok[i]) continue;
        dij[i1][i] = particles[i].dij( particles[i1] );
      }

    } else {
      const p4<alg>& p = particles[i1];

      // identify as jet
      jets.push_back( pp_type( p.px, p.py, p.pz, p.E ) );

      // print clustering step
      #ifdef __cluster_debug
      std::cout << std::setw(3) << p.id
                << " is a Jet | d = " << dist << std::endl;
      #endif

      // "remove"
      ok[i1] = false;
    }

    --n_ok;
  } // end while

  #ifdef __cluster_debug
  std::cout << std::endl;
  #endif

  return jets;
}

} // end namespace
