#ifndef TORUS_GRID_HPP_INCLUDED
#define TORUS_GRID_HPP_INCLUDED

// simple grid class with torus-topology
// 
// all bugs are mine: Hanno 2021


#include <vector>
#include "torus.hpp"


namespace torus {

  template <typename T>
  class grid_t
  {
  public:
    grid_t(grid_t&&) = default;
    grid_t& operator=(grid_t&&) = default;
    grid_t(const grid_t&) = default;
    grid_t& operator=(const grid_t&) = default;

    grid_t(size_t S) : S_(S), data_(S * S) {} 
    grid_t(size_t S, T val) : S_(S), data_(S * S, val) {}

    size_t size() const noexcept { return S_ * S_; }
    size_t wh() const noexcept { return S_; }

    const T& operator()(const vec_t& coor) const noexcept
    {
      const auto wp = wrap(coor);
      const auto s = static_cast<float>(S_);
      const auto ix = static_cast<size_t>(s * wp[0]);
      const auto iy = static_cast<size_t>(s * wp[1]);
      assert(ix < S_);
      assert(iy < S_);
      return *(data_.data() + iy * S_ + ix);
    }

    T& operator()(const vec_t& coor) noexcept
    {
      const auto wp = wrap(coor);
      const auto s = static_cast<float>(S_);
      const auto ix = static_cast<size_t>(s * wp[0]);
      const auto iy = static_cast<size_t>(s * wp[1]);
      assert(ix < S_);
      assert(iy < S_);
      return *(data_.data() + iy * S_ + ix); // still do not understand the intelligent way to retrive the data contained in a certain grid cell
    }

    // my scary try to destroy everthing (seeing if two individuals are located on the same cell)
    size_t retriveCell(const vec_t& coor) noexcept
    {
        const auto wp = wrap(coor);
        const auto s = static_cast<float>(S_);
        const auto ix = static_cast<size_t>(s * wp[0]);
        const auto iy = static_cast<size_t>(s * wp[1]);
        assert(ix < S_);
        assert(iy < S_);
        return iy * S_ + ix;
    }

    float pixel_radius() const noexcept
    {
      return 0.5f / S_;
    }

    vec_t pixel_radii() const noexcept
    {
      return { 0.5f / S_, 0.5f / S_ };
    }

    aabb_t pixel(const vec_t& coor) const noexcept
    {
      const auto pr = pixel_radii();
      const auto pc = wrap(coor + pr);   // pixel center
      return { pc, pr };
    }
     
    auto cbegin() const noexcept { return data_.cbegin(); }
    auto cend() const noexcept { return data_.end(); }
    auto begin() const noexcept { return data_.cbegin(); }
    auto end() const noexcept { return data_.cend(); }
    auto begin() noexcept { return data_.begin(); }
    auto end() noexcept { return data_.end(); }

    const T* data() const noexcept { return data_.data(); }
    T* data() noexcept { return data_.data(); } 

  private:
    size_t S_;
    std::vector<T> data_;
  };

}

#endif
