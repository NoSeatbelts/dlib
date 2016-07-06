#ifndef DSTATS_H
#define DSTATS_H
#include <cmath>
#ifdef __cpluslus
namespace dlib {
#endif

#ifdef __cpluslus
// Convert Pearson's R to Fisher's Z'.
static inline constexpr double r2z(double r) {
    return 0.5 * (std::log((1 + r) / (1 - r)));
}
// Convert Pearson's R to Fisher's Z'.
static inline constexpr double z2r(double z) {
    const double p(std::pow(M_E, 2 * z));
    return (p - 1) / (p + 1);
}
#endif

#ifdef __cpluslus
} // namespace dlib
#endif

#endif /* DSTATS_H */
