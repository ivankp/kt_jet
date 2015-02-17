#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <cmath>

namespace clustering {

// helpful functions
// ******************************************************************

template<typename T> inline T sq(T x) { return x*x; }

#if __cplusplus == 199711L
template<class ForwardIt> inline ForwardIt next(ForwardIt it) { return ++it; }
#endif

// functions to get variables from particle
// ******************************************************************

template<class P> inline double __pt (const P& p) { return p.pt(); }
template<class P> inline double __rap(const P& p) { return p.rap(); }
template<class P> inline double __phi(const P& p) { return p.phi(); }

// clustering type enum
// ******************************************************************

enum algorithm { kt_alg, antikt_alg };

// distance calculations
// ******************************************************************

// distance to the beam
template<algorithm alg, class P>
inline double dist(const P &p) {
  return ( alg==kt_alg ? sq(__pt(p)) : sq(1.0/__pt(p)) );
}

// distance between particles
template<algorithm alg, class P>
inline double dist(const P &a, const P &b) {
  double deltaPhi = __phi(a) - __phi(b);
  if (deltaPhi >  M_PI) deltaPhi = 2*M_PI - deltaPhi;
  if (deltaPhi < -M_PI) deltaPhi = 2*M_PI + deltaPhi;

  return ( alg==kt_alg ? sq( std::min(   __pt(a),   __pt(b)) )
                       : sq( std::min(1./__pt(a),1./__pt(b)) )
         ) * ( sq(__rap(a)-__rap(b)) + sq(deltaPhi) );
}

// structure that ties distance to a particle
template<algorithm alg, class P>
struct pwrap {
  P *p;     // particle
  double d; // distance
  int id;

  pwrap(P *p): p(p), d( dist<alg>(*p) ), id(num++), destroy(true) { }
  pwrap(const pwrap<alg,P> &other)
  : p(other.p), d(other.d), id(other.id), destroy(other.destroy)
  {
    other.destroy = false;
  }
  ~pwrap() {
    //std::cout << "Deleting " << id << std::endl;
    if (destroy) delete p;
  }

  // if equal distances, sort by pt
  bool operator<(const pwrap<alg,P>& rhs) const {
    if (d == rhs.d) return ( __pt(*p) > __pt(*rhs.p) );
    else return ( d < rhs.d );
  }

  void keep() const { destroy = false; }

private:
  mutable bool destroy;
  static int num;
};
template<algorithm alg, class P> int pwrap<alg,P>::num = 0;

template<class T>
struct ordered_pair: public std::pair<T,T> {
  typedef std::pair<T,T> base;
  ordered_pair(const T &a, const T &b): base(a,b) {
    if (base::first > base::second) std::swap(base::first,base::second);
  }
  ordered_pair(const std::pair<T,T> &other): std::pair<T,T>(other) {
    if (base::first > base::second) std::swap(base::first,base::second);
  }
};

// clustering function
// ******************************************************************
template<algorithm alg, class Container>
std::list<typename Container::value_type>
cluster(const Container& pp, double R, bool print_steps=false) {

  typedef typename Container::value_type particle_t;
  typedef typename Container::const_iterator c_iter_t;
  typedef typename std::set< pwrap<alg,const particle_t> >::iterator s_iter_t;

  // collect initial particles into a set
  // sorted by distance to the beam
  std::set< pwrap<alg,const particle_t> > particles;
  for (c_iter_t it=pp.begin(), end=pp.end(); it!=end; ++it)
    particles.insert( &(*it) ).first->keep();

  // map for cashing pairwise distances
  std::map< ordered_pair<const particle_t*>, double > dij;
  for (c_iter_t it=pp.begin(), endi=(--pp.end()); it!=endi; ++it)
    for (c_iter_t jt=next(it), endj=pp.end(); jt!=endj; ++jt)
      dij[std::make_pair( &(*it), &(*jt) )] = dist<alg>(*it,*jt);

  // output list of jets with constituents
  std::list<particle_t> jets;

  // perform clustering iterations until no more particles left
  while (particles.size()) {
    s_iter_t begin = particles.begin();
    s_iter_t end   = particles.end();
    s_iter_t it = begin, it1 = end, it2 = end;

    double _dist = it->d*R*R;
    bool merged = false;

    // find closest pair
    for (; it!=end; ++it) {
      for (s_iter_t jt=next(it); jt!=end; ++jt) {
        double d = dij[std::make_pair( it->p, jt->p )];
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
      const particle_t *p = new particle_t( (*it1->p) + (*it2->p) );

      // print clustering step
      if (print_steps)
        std::cout << "Merged "
                  << std::setw(3) << it1->id << " & "
                  << std::setw(3) << it2->id
                  << ".  d = " << _dist << std::endl;

      // remove merged particles from set
      particles.erase(it1);
      particles.erase(it2);

      // precompute distances
      for (s_iter_t it=particles.begin(), end=particles.end(); it!=end; ++it)
        dij[std::make_pair( it->p, p )] = dist<alg>( *it->p, *p );

      // insert combined pseudo-particle
      particles.insert(p);

    } else {
      // identify as jet
      jets.push_back(*begin->p);

      // print clustering step
      if (print_steps)
        std::cout << begin->id << " is a Jet.  d = " << _dist << std::endl;

      // remove from particles set
      particles.erase(begin);
    }
  } // end while

  return jets;
}

} // end namespace
