#include <iostream>
#include <iomanip>
#include <list>
#include <cmath>

#ifdef __cluster_time
#include "Timer.h"
#endif

namespace clustering {

// helpful functions
// ******************************************************************

template<typename T> inline T sq(T x) { return x*x; }

#if __cplusplus == 199711L
template<class ForwardIt> inline ForwardIt next(ForwardIt it) { return ++it; }
#endif

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

template<class Container>
class sorted: public Container {
public:
  typedef typename Container::iterator   iterator;
  typedef typename Container::value_type value_type;
  sorted(): Container() { }
  iterator insert(const value_type& val) {
    return Container::insert(
      std::lower_bound(Container::begin(), Container::end(), val),
      val);
  }
};

// clustering function
// ******************************************************************
template<algorithm alg, class InputIterator>
std::list<typename InputIterator::value_type>
cluster(InputIterator first, InputIterator last, double R)
{
  // input particle type
  typedef typename InputIterator::value_type in_type;
  // internal list iterator type
  typedef typename std::list< p4<alg> >::iterator iter_t;
  
  #ifdef __cluster_time
  Timer tm;
  tm.start();
  #endif

  // collect initial particles into a list
  // sorted by distance to the beam
  sorted< std::list< p4<alg> > > particles;
  for (InputIterator it=first; it!=last; ++it)
    particles.insert(p4<alg>(__px(*it),__py(*it),__pz(*it),__E(*it)));

  #ifdef __cluster_debug
  for (iter_t it=particles.begin(), end=particles.end(); it!=end; ++it)
    std::cout << it->id << ": " << it->d << std::endl;
  std::cout << std::endl;
  #endif

  // output list of jets with constituents
  std::list<in_type> jets;

  // perform clustering iterations until no more particles left
  while (particles.size()) {
    iter_t begin = particles.begin();
    iter_t end   = particles.end();
    iter_t it = begin, it1 = end, it2 = end;

    double _dist = it->d*R*R;
    bool merged = false;

    // find closest pair
    for (; it!=end; ++it) {
      for (iter_t jt=next(it); jt!=end; ++jt) {
        double d = it->dij(*jt);
        if (d < _dist) {
          _dist = d;
          it1 = it;
          it2 = jt;
          if (!merged) merged = true;
        }
      }
    }

    if (merged) {
      // merge particles
      #ifdef __cluster_debug
      iter_t p =
      #endif
      particles.insert( *it1 + *it2 );

      // print clustering step
      #ifdef __cluster_debug
      std::cout << std::setw(3) << p->id << ": merged "
                << std::setw(3) << it1->id << " & "
                << std::setw(3) << it2->id
                << " | d = " << _dist << std::endl;
      #endif

      // remove merged particles from set
      particles.erase(it1);
      particles.erase(it2);

    } else {
      // identify as jet
      jets.push_back( in_type(begin->px,begin->py,begin->pz,begin->E) );

      // print clustering step
      #ifdef __cluster_debug
      std::cout << std::setw(3) << begin->id
                << " is a Jet | d = " << _dist << std::endl;
      #endif

      // remove from particles set
      particles.erase(begin);
    }
  } // end while

  #ifdef __cluster_debug
  std::cout << std::endl;
  #endif
  
  #ifdef __cluster_time
  tm.stop();
  std::cout << "Clustering time: " << tm.duration() << " ms\n" << std::endl;
  #endif

  return jets;
}

} // end namespace
