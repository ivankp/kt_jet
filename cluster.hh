#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <cmath>

namespace clustering {

template<typename T> inline T sq(T x) { return x*x; }

template<class ForwardIt> ForwardIt next(ForwardIt it) { return ++it; }

enum algorithm { kt_alg, antikt_alg };

// distance to the beam
template<algorithm alg, class P>
inline double dist(const P &p) {
  return ( alg==kt_alg ? sq(p->pt()) : sq(1.0/p->pt()) );
}

// distance between particles
template<algorithm alg, class P>
inline double dist(const P &a, const P &b) {
  double deltaPhi = a.phi() - b.phi();
  if (deltaPhi >  M_PI) deltaPhi = 2*M_PI - deltaPhi;
  if (deltaPhi < -M_PI) deltaPhi = 2*M_PI + deltaPhi;
  
  return ( alg==kt_alg ? sq( min(   a.pt(),   b.pt()) )
                       : sq( min(1./a.pt(),1./b.pt()) )
         ) * ( sq(a.rap()-b.rap()) + sq(deltaPhi) );
}

// structure that ties distance to a particle
template<algorithm alg, class P>
struct add_dist {
  P *p;     // particle
  double d; // distance
  
  dist(P *p): p(p), d( dist<alg>(p) ) { }
  
  // if equal distances, sort by pt
  bool operator(const add_dist<P,kt>& rhs) const {
    if (d == rhs.d) return ( p->pt() > rhs.p->pt() );
    else return ( d < rhs.d );
  }
};

template<class S> inline void tmp_erase(S& s, S::value_type x) {
  S::iterator it = s.find(x);
  if (it!=s.end()) {
    delete *it;
    s.erase(it);
  }
}

// clustering function
// ******************************************************************
template<algorithm alg, class Container>
std::list<Container::value_type> cluster(const Container& pp, double R) {

  typedef Container::value_type particle_t;
  typedef Container::iterator c_iter_t;
  typedef std::set< add_dist<particle_t,alg> >::iterator s_iter_t;
  
  // collect initial particles into a set
  // sorted by distance to the beam
  std::set< add_dist<particle_t,alg> > particles;
  for (c_iter_t it=pp.begin(), end=pp.end(); it!=end; ++it)
    particles.insert( &(*it) );

  // map for cashing pairwise distances
  std::map< std::pair<particle_t*,particle_t*>, double > dij;
  for (c_iter_t it=pp.begin(), end=pp.end(); it!=end; ++it)
    for (c_iter_t jt=next(it); jt!=end; ++jt)
      dij[std::make_pair( &(*it), &(*jt) )] = dist<alg>( &(*it), &(*jt) );
      
  // output list of jets with constituents
  std::list<particle_t> jets;
      
  // collect intermediate particles to prevent memory leaks
  std::set<particle_t*> tmp;
  
  // perform clustering iterations until no more particles left
  while (particles.size()) {
    s_iter_t begin = particles.begin();
    s_iter_t end   = particles.end();
    s_iter_t it = begin, it1 = end, it2 = end;
    
    double _dist = it->d*R*R;
    bool merged = false;
    
    for (; it!=end; ++it) {
      for (s_iter_t jt=next(it); jt!=end; ++jt) {
        double d = dij[std::make_pair( &(*it), &(*jt) )];
        if (d < _dist) {
          dist = d;
          it1 = it;
          it2 = jt;
          if (!merged) merged = true;
        }
      }
    }
    
    if (merged) {
      // insert combined pseudo-particle
      particles.insert( 
        *tmp.insert( new particle_t( (*it1->p) + (*it2->p) ) ).first
      );
      
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
  
}
