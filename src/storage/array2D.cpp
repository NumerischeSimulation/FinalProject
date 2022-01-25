#include "storage/array2D.h"

#include <cassert>

Array2D::Array2D(std::array<int,2> size) :
  size_(size)
{
  // allocate data, initialize to 0
  data_.resize(size_[0]*size_[1], 0.0);
}

//! get the size
std::array<int,2> Array2D::size() const
{
  return size_;
}

double &Array2D::operator()(int x, int y)
{
  // std::cout << "x " << x << "y " << y << std::endl;
  int i = x;
  int j = y;

  const int index = j*size_[0] + i;

  // assert that indices are in range
  assert(0 <= i && i < size_[0]);
  assert(0 <= j && j < size_[1]);
  assert(j*size_[0] + i < (int)data_.size());

  return data_[index];
}

double Array2D::operator()(int x, int y) const
{
  // std::cout << "x " << x << "y " << y << std::endl;
  int i = x;
  int j = y;

  const int index = j*size_[0] + i;

  // assert that indices are in range
  assert(0 <= i && i < size_[0]);
  assert(0 <= j && j < size_[1]);
  assert(j*size_[0] + i < (int)data_.size());

  return data_[index];
}
