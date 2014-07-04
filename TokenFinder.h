#pragma once
// TokenFinder.h: interface for the TokenFinder class.
//
//////////////////////////////////////////////////////////////////////

#include <string>

class TokenFinder  
{
public:

	TokenFinder(){}

	TokenFinder(const char *del):delimiter(del){}

	int operator ()(const char *&start , const char *&end)
	{
		const char *pPos = 0;

		if(!delimiter.empty())
		{
			pPos =  strstr(end , delimiter.c_str());

			while(pPos == end)
			{
				end+=delimiter.size();
				pPos = strstr(end , delimiter.c_str());
			}
		}

		start = end;
		end = (pPos) ? pPos : end + strlen(start);
		return (end-start);
	}

	~TokenFinder(){}
	
private:

	std::string delimiter;
};
