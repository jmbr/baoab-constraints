#ifndef LINSPACE_H
#define LINSPACE_H

#include <vector>

template <typename T = double>
std::vector<T> linspace(T a, T b, size_t num_elems) {
  T h = (b - a) / static_cast<T>(num_elems - 1);
  std::vector<T> xs(num_elems);
  typename std::vector<T>::iterator x;
  T val;
  for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
    *x = val;
  return xs;
}

#endif
