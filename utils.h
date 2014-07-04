#pragma once

#include <vector>

class utils
{
public:


	bool checkTour(std::vector<int> &t , int size)
    {
		int i, sum = 0;
		for (i = 0; i < t.size(); i++) 
		{
			sum += t[i];
		}
		if (sum != (size - 1) * size / 2) 
		{
			return false;
		}
		return true;
    }


	utils(void)
	{
	}

	~utils(void)
	{
	}
};

