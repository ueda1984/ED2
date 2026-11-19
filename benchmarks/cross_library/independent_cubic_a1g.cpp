// Independent C++17 physical oracle, NOT ED2, QuSpin or XDiag.
// Reconstructs the Gamma/O_h A1g (trivial representation) 3-magnon block
// for the exact 9N-bond spin-half J1-J2 XXZ cubic model in ED2 Draft PR16.
// Output is Matrix Market coordinate real GENERAL (indices 1-based).
//
// Two independent checks are performed before writing the matrix:
//   sum_{orbits} |orbit| == C(N,3), and max|H_ij-H_ji| < 1e-9.
// No eigenvalues are claimed here; use independent_cubic_solver.py.
// Source-of-physics: ED2 PR16 pinned canonical_bonds, 2026-11-18.
// Do not treat source creation time as a cross-library speed benchmark.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

struct Pos { int x,y,z; };
struct Bond { int i,j; double xy,zz; };
struct Point { std::array<int,3> p,s; };
struct Edge { int j; double xy,zz; };
struct Orb { int a,b; long long size; int anchored; };

struct Builder {
  int L,N;
  std::vector<Pos> coords;
  std::vector<Point> points;
  std::vector<Bond> bonds;
  std::vector<std::vector<Edge>> neighbors;
  std::vector<std::unordered_map<int,double>> diag_map;
  std::vector<int> owner;
  std::vector<Orb> orbs;
  double allup=0;

  explicit Builder(int side): L(side), N(L*L*L), coords(N),
      neighbors(N), diag_map(N), owner((N-1)*(N-2)/2,-1) {
    if(L<4 || L>10 || L%2)
      throw std::runtime_error("side must be 4,6,8,10");
    for(int i=0;i<N;++i) coords[i] = {i%L,(i/L)%L,i/(L*L)};
    std::array<int,3> p={0,1,2};
    do {
      for(int s0:{-1,1}) for(int s1:{-1,1}) for(int s2:{-1,1})
        points.push_back({p,{s0,s1,s2}});
    } while(std::next_permutation(p.begin(),p.end()));
    if(points.size()!=48) throw std::runtime_error("point group size");
    auto add=[&](int i,int j,double xy,double zz) {
      bonds.push_back({i,j,xy,zz});
      neighbors[i].push_back({j,xy,zz});
      neighbors[j].push_back({i,xy,zz});
      diag_map[i][j]=zz;
      diag_map[j][i]=zz;
      allup+=zz*.25;
    };
    const std::array<Pos,3> ax={{{1,0,0},{0,1,0},{0,0,1}}};
    const std::array<Pos,6> f={{{1,1,0},{1,-1,0},{1,0,1},
                                {1,0,-1},{0,1,1},{0,1,-1}}};
    for(int z=0;z<L;++z) for(int y=0;y<L;++y) for(int x=0;x<L;++x) {
      int i=site(x,y,z);
      for(auto d:ax) add(i,site(x+d.x,y+d.y,z+d.z),1.,.77);
      for(auto d:f) add(i,site(x+d.x,y+d.y,z+d.z),.31,.19);
    }
    if((int)bonds.size()!=9*N) throw std::runtime_error("bond count");
    for(int i=0;i<N;++i)
      if(neighbors[i].size()!=18) throw std::runtime_error("site degree");
  }
  int wrap(int x) const { x%=L; return x<0?x+L:x; }
  int site(int x,int y,int z)const {
    return wrap(x)+L*(wrap(y)+L*wrap(z));
  }
  int point_site(int i,const Point &p)const {
    auto c=coords[i];int a[3]={c.x,c.y,c.z};
    return site(p.s[0]*a[p.p[0]],p.s[1]*a[p.p[1]],p.s[2]*a[p.p[2]]);
  }
  int diff_site(int i,int j)const {
    auto a=coords[i],b=coords[j];
    return site(b.x-a.x,b.y-a.y,b.z-a.z);
  }
  // 0-based triangular indexing for 1 <= a < b <= N-1.
  int pair_slot(int a,int b)const {
    if(a>b)std::swap(a,b);
    if(a<1 || b>=N || a==b)throw std::runtime_error("invalid anchor pair");
    return ((b-1)*(b-2))/2+(a-1);
  }
  int orbit_of(std::array<int,3> st) const {
    int a=diff_site(st[0],st[1]), b=diff_site(st[0],st[2]);
    int id=owner[pair_slot(a,b)];
    if(id<0)throw std::runtime_error("missing orbit lookup");
    return id;
  }
  void build_orbits() {
    // Quotient translations by anchoring any of 3 occupied sites at origin;
    // there are C(N-1,2) anchored pairs, not C(N,3) raw configurations.
    for(int a=1;a<N-1;++a) for(int b=a+1;b<N;++b) {
      if(owner[pair_slot(a,b)]>=0)continue;
      int id=(int)orbs.size();
      std::vector<int> slots;slots.reserve(144);
      std::array<int,3> st={0,a,b};
      for(auto const &p:points) {
        int t[3]={point_site(st[0],p),point_site(st[1],p),
                  point_site(st[2],p)};
        for(int anchor=0;anchor<3;++anchor) {
          int x=diff_site(t[anchor],t[(anchor+1)%3]);
          int y=diff_site(t[anchor],t[(anchor+2)%3]);
          slots.push_back(pair_slot(x,y));
        }
      }
      std::sort(slots.begin(),slots.end());
      slots.erase(std::unique(slots.begin(),slots.end()),slots.end());
      for(int slot:slots) {
        if(owner[slot]>=0 && owner[slot]!=id)
          throw std::runtime_error("inconsistent orbit partitions");
        owner[slot]=id;
      }
      // Counting (state, chosen anchor) pairs gives N * anchored_size =
      // 3 * full_orbit_size, even for translation-stabilized triples.
      long long raw_count=(long long)N*(long long)slots.size();
      if(raw_count%3)throw std::runtime_error("orbit normalization noninteger");
      orbs.push_back({a,b,raw_count/3,(int)slots.size()});
    }
    for(int id:owner)if(id<0)throw std::runtime_error("unmapped pair");
    long long raw=0;for(auto o:orbs)raw+=o.size;
    long long target=(long long)N*(N-1)*(N-2)/6;
    std::cerr<<"ORBIT_CHECK side="<<L<<" D="<<orbs.size()
             <<" raw="<<raw<<" target="<<target
             <<" Eall="<<std::setprecision(15)<<allup<<"\n";
    if(raw!=target)throw std::runtime_error("orbit size sum mismatch");
  }
  double diagonal(const std::array<int,3>&s)const {
    double e=allup;
    for(int i:s) for(auto edge:neighbors[i]) e-=.5*edge.zz;
    for(int k=0;k<3;++k)for(int l=k+1;l<3;++l) {
      auto it=diag_map[s[k]].find(s[l]);
      if(it!=diag_map[s[k]].end())e+=it->second;
    }
    return e;
  }
  std::vector<std::unordered_map<int,double>> matrix()const {
    int d=(int)orbs.size();
    std::vector<std::unordered_map<int,double>> columns(d);
    for(int col=0;col<d;++col) {
      auto o=orbs[col];
      std::array<int,3> s={0,o.a,o.b};
      columns[col][col]+=diagonal(s);
      // Applying H to one configuration per orbit is sufficient by
      // symmetry; the normalization coefficient is sqrt(|O_col|/|O_row|).
      for(int pos=0;pos<3;++pos) {
        int i=s[pos];
        for(auto e:neighbors[i]) {
          if(e.j==s[(pos+1)%3] || e.j==s[(pos+2)%3])continue;
          auto t=s;t[pos]=e.j;
          int row=orbit_of(t);
          columns[col][row]+=.5*e.xy*
              std::sqrt(double(o.size)/double(orbs[row].size));
        }
      }
    }
    double err=0;size_t nnz=0;
    for(int col=0;col<d;++col) {
      nnz+=columns[col].size();
      for(auto [row,v]:columns[col]) {
        auto it=columns[row].find(col);
        double other=(it==columns[row].end()?0.:it->second);
        err=std::max(err,std::abs(v-other));
      }
    }
    std::cerr<<"MATRIX_CHECK D="<<d<<" entries="<<nnz
             <<" max_asym="<<std::setprecision(16)<<err<<"\n";
    if(err>1e-9)throw std::runtime_error("non-Hermitian projection");
    return columns;
  }
};
int main(int argc,char**argv)try {
  if(argc!=3) {
    std::cerr<<"Usage: independent_cubic_a1g SIDE output.mtx\n";
    return 2;
  }
  int L=std::stoi(argv[1]);
  Builder b(L);b.build_orbits();auto cols=b.matrix();
  std::ofstream o(argv[2]);
  if(!o)throw std::runtime_error("cannot write matrix");
  size_t nnz=0;for(auto const &col:cols)nnz+=col.size();
  o<<"%%MatrixMarket matrix coordinate real general\n"
    <<"% independent cubic J1-J2 XXZ Gamma Oh A1g\n";
  o<<cols.size()<<" "<<cols.size()<<" "<<nnz<<"\n";
  o<<std::setprecision(17);
  for(size_t j=0;j<cols.size();++j)
    for(auto [i,v]:cols[j]) o<<i+1<<" "<<j+1<<" "<<v<<"\n";
  std::cout<<"D="<<cols.size()<<" nnz="<<nnz
           <<" allup="<<std::setprecision(16)<<b.allup<<"\n";
} catch(std::exception const &e) {
  std::cerr<<"ERROR: "<<e.what()<<"\n";
  return 3;
}
