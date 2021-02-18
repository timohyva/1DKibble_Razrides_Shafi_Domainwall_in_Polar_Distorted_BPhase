/*
 * 
 * File: CalculateA20InOneElement.cpp
 *
 * 
 * 
 */

#include <iostream>
#include <cstddef>
#include <cstdlib>
#include "ClassCalculateA20InOneElement.hpp"


void CalculateA20::CalculateA20InOneElement(double xi,double xi1)
{
  double a20_idx_0;
  double a20_idx_3;
  double dNi1dx;
  double dNidx;
  double dNidx_tmp;
  
  dNidx_tmp = xi1 - xi;
  dNidx = -1.0 / dNidx_tmp;
  dNi1dx = 1.0 / dNidx_tmp;

  a20_idx_0 = 0.5 * dNidx * dNidx * dNidx_tmp;
  a20_idx_3 = 0.5 * dNi1dx * dNi1dx * dNidx_tmp;
 
  dNidx_tmp = 0.5 * (0.5 * dNi1dx * dNidx * dNidx_tmp + 0.5 * dNidx * dNi1dx * dNidx_tmp);
  
  //  matrix symmetrized 
  AInOneElement[0] = 0.5 * (a20_idx_0 + a20_idx_0);
  AInOneElement[1] = dNidx_tmp;
  AInOneElement[2] = dNidx_tmp;
  AInOneElement[3] = 0.5 * (a20_idx_3 + a20_idx_3);
  
  arr=AInOneElement; //transfer the pointer of array A20 to object data member arr
}

const double * const CalculateA20::returnResult(size_t length) const
{
   
    return this->arr;
}