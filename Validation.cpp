#include "Validation.h"

bool validate_best(double L,double U,double rel_error,double& val_R)
{
	if (fabs(L) < small_epsilon && fabs(U) < small_epsilon)
	{
		val_R = 0;
		return true;
	}
	if(U-L<=rel_error*(U+L))
	{
		val_R=2*L*U/(L+U);
		return true;
	}

	return false;
}