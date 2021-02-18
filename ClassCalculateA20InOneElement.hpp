/* 
 * header file for Class of A matrix Calculation
 *
 *
 */


#ifndef CALCULATEA20INONEELEMENT_HPP
#define CALCULATEA20INONEELEMENT_HPP

#include <iostream>
#include <cstddef>
#include <cstdlib>

class CalculateA20 {
     
      private:
         double *arr;
         double AInOneElement[4];
         size_t length;    
        
      public:
          // double xi;
          // double xi1;
          
          
          // calculate the A20 for one element
          void CalculateA20InOneElement(double xi,double xi1);
          
          // return the pointer of the result array
          const double * const returnResult(size_t length) const;
       
};
#endif