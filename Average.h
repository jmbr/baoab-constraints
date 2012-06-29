#ifndef AVERAGE_H
#define AVERAGE_H

class Average {
 public:
  Average();

  void update(double input);

  long double operator()() const;

 private:
  unsigned long long count;
  long double sum;
  long double c;
};

#endif
