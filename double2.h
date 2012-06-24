#ifndef DOUBLE2_H
#define DOUBLE2_H

#include <cmath>

#include <ostream>
// #include <initializer_list>

#ifdef HAVE_SSE
typedef double v2df __attribute__ ((vector_size (2 * sizeof(double))));
// #else
// # warning "SIMD (SSE) operations disabled."
#endif

// TODO Take advantage of C++11's move semantics to reduce the need
// for temporaries.

struct double2 {
  union {
#ifdef HAVE_SSE
    v2df data;
#endif
    double array[2];
    struct { double x, y; };
  };

  double2() : x(0.0), y(0.0) {}

  double2(double x_, double y_) : x(x_), y(y_) {}

  // double2(std::initializer_list<double> list) {
  //   std::copy(list.begin(), list.end(), array);
  // }

  inline double& operator[](size_t i) {
    return array[i];
  }

  inline const double& operator[](size_t i) const {
    return array[i];
  }

  inline double2& operator+=(const double2& d) {
#ifdef HAVE_SSE
    this->data += d.data;
#else
    x += d.x;
    y += d.y;
#endif
    return *this;
  }

  inline double2& operator-=(const double2& d) {
#ifdef HAVE_SSE
    this->data -= d.data;
#else
    x -= d.x;
    y -= d.y;
#endif
    return *this;
  }

  inline double2& operator*=(const double2& d) {
#ifdef HAVE_SSE
    this->data *= d.data;
#else
    x *= d.x;
    y *= d.y;
#endif
    return *this;
  }

  inline double2& operator/=(const double2& d) {
#ifdef HAVE_SSE
    this->data /= d.data;
#else
    x /= d.x;
    y /= d.y;
#endif
    return *this;
  }

  inline double2& operator*=(const double s) {
#ifdef HAVE_SSE
    this->x *= s;
    this->y *= s;
#else
    x *= s;
    y *= s;
#endif
    return *this;
  }

  inline double2& operator/=(const double s) {
#ifdef HAVE_SSE
    this->x /= s;
    this->y /= s;
#else
    x /= s;
    y /= s;
#endif
    return *this;
  }

  inline double norm2() const {
    // XXX Replace this by calls to __builtin_ia32_dppd
    return x*x + y*y;
  }

  inline double norm() const {
    // XXX Provide an optional variant using rsqrt on platforms that
    // allow it.
    return sqrt(norm2());
  }

  inline void rotate(double theta) {
    double s, c;
#ifdef _GNU_SOURCE
    sincos(theta, &s, &c);
#else
    s = sin(theta);
    c = cos(theta);
#endif
    const double rx = c * x - s * y;
    const double ry = s * x + c * y;
    x = rx;
    y = ry;
  }

  friend std::ostream& operator<<(std::ostream& stream, const double2& d) {
    return stream << d.x << " " << d.y;
  }
};

inline double2 operator+(const double2& a, const double2& b) {
  double2 c = a;
  return c += b;
}

inline double2 operator-(const double2& a, const double2& b) {
  double2 c = a;
  return c -= b;
}

inline double2 operator-(const double2& a) {
  double2 b = a;
  return b *= -1.0;
}

inline double2 operator*(const double2& a, const double2& b) {
  double2 c = a;
  return c *= b;
}

inline double2 operator*(const double s, const double2& v) {
  double2 w = v;
  return w *= s;
}

inline double2 operator*(const double2& v, const double s) {
  double2 w = v;
  return w *= s;
}

inline double2 operator/(const double2& a, const double2& b) {
  double2 c = a;
  return c /= b;
}

inline double2 operator/(const double s, const double2& v) {
  double2 w = v;
  return w /= s;
}

inline double2 operator/(const double2& v, const double s) {
  double2 w = v;
  return w /= s;
}

inline double dot(const double2& u, const double2& v) {
  const double2 d = u * v;
  // XXX Replace this by calls to __builtin_ia32_dppd
  return d.x + d.y;
}

inline double dist(const double2& u, const double2& v) {
  const double2 r = u - v;
  return r.norm();
}

#endif
