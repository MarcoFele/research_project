#ifndef TORUS_HRTREE_HPP_INCLUDED
#define TORUS_HRTREE_HPP_INCLUDED

// basic Hilbert Rtree for torus::aabb_t
//
// all bugs are mine: Hanno 2021


#include <vector>
#include <hrtree/isfc/hilbert.hpp>
#include <hrtree/isfc/key_gen.hpp>
#include <hrtree/sorting/radix_sort.hpp>
#include <hrtree/rtree.hpp>
#include "torus.hpp"


// this macro specialize a traits class required by the hrtree library
// defined in hrtree/adapt_point.hpp.
HRTREE_ADAPT_POINT_FUNCTION(torus::vec_t, float, 2, (float*)std::addressof)


namespace torus {


  namespace detail {
  
    // used internally by hrtree_t.
    // the rtree needs to know how to generate bounding volumes over
    // bounding volumes. here, we want to use torus::aabb_t as bonding 
    // volumes. Thus, we use torus::include as the basic primitive. 
    struct aabb_build_policy
    {
      template <typename IT, typename OIT>
      void operator()(IT first, IT last, OIT out) const
      {
        torus::aabb_t bbox = *out;
        for (; first != last; ++first) {
          bbox = torus::include(bbox, *first);
        }
        *out = bbox;
      }
    };
  
    // Hilbert values and radix-sort stuff
    using key_t = hrtree::hilbert<2, 15>::type;          // 2D 'Hilbert value' of order 15
    using keygen_t = hrtree::key_gen_01<key_t, vec_t>;   // pick a generator
    using keyidx_t = std::pair<key_t::word_type, uint32_t>;

    // radix sort is a bit of a pain (but fast)
    // works with the bytes of the memory representation and
    // needs to be told where to find them 
    struct keyidx_conv_t
    {
      static constexpr int key_bytes = key_t::key_bytes;
      const std::uint8_t* operator()(const keyidx_t& x) const 
      {
        return reinterpret_cast<const std::uint8_t*>(std::addressof(x.first));
      }
    };

  }


  class hrtree_t
  {
  public:
    using index_t = int32_t;

    hrtree_t() {}

    template <typename RaIt, typename Conv>
    void build(RaIt first, RaIt last, Conv conv);

    template <typename Fun>
    void query(const aabb_t& bbox, Fun fun) const; // how can it be constant if fun changes the resources inside the prey??

  private:
    hrtree::rtree<aabb_t, detail::aabb_build_policy, 8> hrtree_;
    std::vector<detail::keyidx_t> ki_;      
    std::vector<detail::keyidx_t> ki_buf_;  // some more that is needed by radix-sort
  };


  template <typename RaIt, typename Conv>
  void hrtree_t::build(RaIt first, RaIt last, Conv conv)
  {
    const auto N = static_cast<index_t>(std::distance(first, last));
    ki_.resize(N);
    ki_buf_.resize(N);
    hrtree_.build_index(N);   // i.e. allocate memory for our leaves
    if (N) {
      // generate <Hilbert value, index> pairs
      detail::keygen_t keygen{};
      for (index_t i = 0; i < N; ++i) {
        ki_[i] = { keygen(conv(first[i]).center).asWord(), i };
      }
      // sort by Hilbert values
      if (hrtree::radix_sort(ki_.begin(), ki_.end(), ki_buf_.begin(), detail::keyidx_conv_t{})) {
        ki_.swap(ki_buf_);
      }
      // store leaves in Hilbert value order
      auto dummy = hrtree_.level_begin(0);
      for (index_t i = 0; i < N; ++i) {
        hrtree_.leaf_bv(i) = conv(first[ki_[i].second]);
      }
      hrtree_.build_hierarchy();
    }
  }


  template <typename Fun>
  void hrtree_t::query(const aabb_t& bbox, Fun fun) const
  {
    auto wfun = [fun=fun, it=ki_.cbegin()](size_t i) { fun((it + i)->second); };
    hrtree_.query( 
      [bbox=bbox](const aabb_t& rhs) { 
        return intersects(bbox, rhs); 
      },
      wfun
    );            
  }


  // for comparison ;)
  class brute_force_t
  {
  public:
    brute_force_t() {}

    template <typename RaIt, typename Conv>
    void build(RaIt first, RaIt last, Conv conv)
    {
      leaves_.clear();
      for (; first != last; ++first) {
        leaves_.push_back(conv(*first));
      }
    }

    template <typename Fun>
    void query(const aabb_t& bbox, Fun fun) const
    {
      for (size_t i = 0; i < leaves_.size(); ++i) {
        if (intersects(leaves_[i], bbox)) {
          fun(i);
        }
      }
    }

  private:
    std::vector<aabb_t> leaves_;
  };
}

#endif
