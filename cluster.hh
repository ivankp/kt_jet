#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <cmath>

namespace clustering {

// helpful functions
// ******************************************************************

template<typename T> inline T sq(T x) { return x*x; }

template<class ForwardIt> ForwardIt increment(ForwardIt it) { return ++it; }

template<class S> inline void tmp_erase(S& s, typename S::value_type x) {
  typename S::iterator it = s.find(x);
  if (it!=s.end()) {
    delete x;
    s.erase(it);
  }
}

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
struct add_dist {
  P *p;     // particle
  double d; // distance
  int id;
  
  add_dist(P *p): p(p), d( dist<alg>(*p) ), id(num++) { }
  
  // if equal distances, sort by pt
  bool operator<(const add_dist<alg,P>& rhs) const {
    if (d == rhs.d) return ( __pt(*p) > __pt(*rhs.p) );
    else return ( d < rhs.d );
  }
  
private:
  static int num;
};
template<algorithm alg, class P> int add_dist<alg,P>::num = 0;

// clustering function
// ******************************************************************
template<algorithm alg, class Container>
std::list<typename Container::value_type>
cluster(const Container& pp, double R) {

  typedef typename Container::value_type particle_t;
  typedef typename Container::const_iterator c_iter_t;
  typedef typename std::set< add_dist<alg,const particle_t> >::iterator s_iter_t;
  
  // collect initial particles into a set
  // sorted by distance to the beam
  std::set< add_dist<alg,const particle_t> > particles;
  for (c_iter_t it=pp.begin(), end=pp.end(); it!=end; ++it)
    particles.insert( &(*it) );

  // map for cashing pairwise distances
  std::map< std::pair<const particle_t*,const particle_t*>, double > dij;
  for (c_iter_t it=pp.begin(), end=pp.end(); it!=end; ++it)
    for (c_iter_t jt=increment(it); jt!=end; ++jt)
      dij[std::make_pair( &(*it), &(*jt) )] = dist<alg>(*it,*jt);
      
  // output list of jets with constituents
  std::list<particle_t> jets;
      
  // collect intermediate particles to prevent memory leaks
  std::set<const particle_t*> tmp;
  
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
      // insert combined pseudo-particle
      std::cout << particles.insert( 
        *tmp.insert( new particle_t( (*it1->p) + (*it2->p) ) ).first
      ).first->d << std::endl;
      
      // test
      std::cout << "Merged " << it1->id << " & " << it2->id << std::endl;
      
      // prevent memory leaks
      tmp_erase(tmp,it1->p);
      tmp_erase(tmp,it2->p);
      
      // remove from particles set
      particles.erase(it1);
      particles.erase(it2);
      
    } else {
      // identify as jet
      jets.push_back(*begin->p);

      // prevent memory leaks
      tmp_erase(tmp,begin->p);

      // remove from particles set
      particles.erase(begin);
    }
  } // end while
  
  return jets;
}
  
} // end namespace
