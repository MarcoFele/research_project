#ifndef TORUS_HPP_INCLUDED
#define TORUS_HPP_INCLUDED

// some math-fun on a normalized torus, [0,1), [0,1)
// 
// all bugs are mine: Hanno 2021


#include <cassert>
#include <ostream>
#include <array>
#include <cmath>
#include <limits>


namespace torus {

  using vec_t = std::array<float, 2>;


  // axis-aligned bounding box
  struct aabb_t
  {
    vec_t center;
    vec_t radii;    // half extent in x/y
  };


  inline vec_t operator+(const vec_t& a, const vec_t& b) noexcept
  {
    return { a[0] + b[0], a[1] + b[1] };
  }

  inline vec_t operator+(const vec_t& a, float s) noexcept
  {
    return { a[0] + s, a[1] + s };
  }

  inline vec_t operator+(float s, const vec_t& a) noexcept
  {
    return { a[0] + s, a[1] + s };
  }

  inline vec_t operator-(const vec_t& a, const vec_t& b) noexcept
  {
    return { a[0] - b[0], a[1] - b[1] };
  }

  inline vec_t operator-(const vec_t& a, float s) noexcept
  {
    return { a[0] - s, a[1] - s };
  }

  inline vec_t operator-(float s, const vec_t& a) noexcept
  {
    return { a[0] - s, a[1] - s };
  }

  inline vec_t operator*(const vec_t& a, const vec_t& b) noexcept 
  {
    return { a[0] * b[0], a[1] * b[1] };
  }

  inline vec_t operator*(const vec_t& a, float s) noexcept 
  {
    return { s * a[0], s * a[1] };
  }

  inline vec_t operator*(float s, const vec_t& a) noexcept 
  {
    return a * s;
  }

  inline vec_t operator/(const vec_t& a, const vec_t& b)
  {
    return { a[0] / b[0], a[1] / b[1] };
  }

  inline vec_t operator/(const vec_t& a, float s)
  {
    return { a[0] / s, s * a[1] / s };
  }

  // element wise minimum
  inline vec_t min(const vec_t& a, const vec_t& b) noexcept
  {
    return { std::min(a[0], b[0]), std::min(a[1], b[1]) };
  }

  // element wise maximum
  inline vec_t max(const vec_t& a, const vec_t& b) noexcept
  {
    return { std::max(a[0], b[0]), std::max(a[1], b[1]) };
  }

  inline vec_t abs(const vec_t& a) noexcept
  {
    return { std::abs(a[0]), std::abs(a[1]) };
  }


  inline std::ostream& operator<<(std::ostream& os, const vec_t& pt)
  {
    os << pt[0] << ',' << pt[1];
    return os;
  }

  inline std::ostream& operator<<(std::ostream& os, const aabb_t& box)
  {
    os << box.center << ',' << box.radii;
    return os;
  }


  namespace detail {

    inline float wrap_point_coor(float x) noexcept
    {
      return x - std::floor(x);   // fractional part
    }

    inline float wrap_ofs_coor(float x) noexcept
    {
      if (x < -0.5f) x += 1;
      else if (x >= 0.5f) x -= 1;
      return x;
    }

  }


  inline bool is_wrapped(const vec_t& pt)
  {
    return pt[0] >= 0 && pt[0] <= 1 && pt[1] >= 0 && pt[1] <= 1;
  }


  // wrap the point pt into [0,1), [0,1)
  inline vec_t wrap(const vec_t& pt) noexcept
  {
    return { detail::wrap_point_coor(pt[0]), detail::wrap_point_coor(pt[1]) };
  }


  // minimal offset
  // points shall be wrapped 
  inline vec_t offset(const vec_t& a, const vec_t& b) noexcept
  {
    assert(is_wrapped(a) && is_wrapped(b));
    const auto dir = a - b;
    return { detail::wrap_ofs_coor(dir[0]), detail::wrap_ofs_coor(dir[1]) };
  }


  // minimal distance squared
  // points shall be wrapped
  inline float distance2(const vec_t& a, const vec_t& b) noexcept
  {
    const auto ofs = offset(a, b);
    return ofs[0] * ofs[0] + ofs[1] * ofs[1];
  }


  // minimal distance
  // points shall be wrapped
  inline float distance(const vec_t& a, const vec_t& b)
  {
    return std::sqrt(distance2(a, b));
  }


  // All torus-operations incur rounding errors.
  // For the intersection-tests below, we favor 'false positives'
  // over 'false negatives'. Thus, we bump the radii by a small amount
  // which defaults to reps:  
  constexpr float reps = 4.f * std::numeric_limits<float>::epsilon();


  // returns true if bbox intersects pt
  // bbox.center and pt shall be wrapped.
  inline bool intersects(const aabb_t& bbox, const vec_t& pt, float eps = reps)
  {
    const auto aofs = abs(offset(pt, bbox.center));
    return (aofs[0] <= (bbox.radii[0] + reps)) && (aofs[1] <= (bbox.radii[1] + reps));
  }


  // returns true if a intersects b
  // centers shall be wrapped.
  inline bool intersects(const aabb_t& a, const aabb_t& b, float = reps)
  {
    const auto aofs = abs(offset(b.center, a.center));
    const auto rr = a.radii + b.radii;
    return (aofs[0] <= (rr[0] + reps)) && (aofs[1] <= (rr[1] + reps));
  }


  // returns the minimal bounding box that contains bbox and pt
  // bbox.center and pt shall be wrapped.
  inline aabb_t include(const aabb_t& bbox, const vec_t& pt) noexcept
  {
    assert(bbox.radii[0] >= 0 && bbox.radii[1] >= 0);
    const auto ofs = offset(pt, bbox.center);
    const auto p = bbox.center + ofs;
    const auto lo = min(bbox.center - bbox.radii, p);
    const auto hi = max(bbox.center + bbox.radii, p);
    const aabb_t ret{ wrap(0.5f * (hi + lo)), 0.5f * (hi - lo) };
    return ret;
  }


  // returns the minimal bounding box that contains the boxes a and b
  // centers shall be wrapped.
  inline aabb_t include(const aabb_t& a, const aabb_t& b) noexcept
  {
    assert(a.radii[0] >= 0 && a.radii[1] >= 0);
    assert(b.radii[0] >= 0 && b.radii[1] >= 0);
    const auto ofs = offset(b.center, a.center);
    const auto cb = a.center + ofs;
    const auto lo = min(a.center - a.radii, cb - b.radii);
    const auto hi = max(a.center + a.radii, cb + b.radii);
    const aabb_t ret{ wrap(0.5f * (hi + lo)), 0.5f * (hi - lo) };
    return ret;
  }

}

#endif
