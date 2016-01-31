#ifndef MATH_UTIL_H
#define MATH_UTIL_H

inline int ipow(int base, int exp)
{
	int result = 1;
	while (exp)
	{
		if (exp & 1)
			result *= base;
		exp >>= 1;
		base *= base;
	}

	return result;
}

#endif /* MATH_UTIL_H */
