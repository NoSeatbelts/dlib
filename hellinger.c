#include "hellinger.h"

double Hellinger_in_c(double* arr1, double* arr2, size_t length){
    double cumSum = 0.;
    double tmpFloat;
    for(size_t index = 0; index < length; index++){
        tmpFloat = abs(sqrt(arr1[index]) - sqrt(arr2[index]));
        cumSum += tmpFloat * tmpFloat;
    }
    return sqrt(cumSum) * M_SQRT1_2;
}
