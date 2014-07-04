#pragma once

#include <iterator>
#include <string>

template<class TFinder>

class Tokenizer : public std::iterator<std::input_iterator_tag , std::string>  
{

public:

	static const char EndOfString;
	static const Tokenizer<TFinder> EOS;
	
	Tokenizer() : start(&EndOfString) , end(&EndOfString) , length(0)
	{
	}

	
	Tokenizer(const char *str , TFinder finder= TFinder())
		: FindToken(finder) , end(str)
	{
		length =FindToken(start , end);
	}
	
	
	Tokenizer& operator =(const char *const str)
	{
		end= str;
		length = FindToken(start , end);
		return *this;
	}

	Tokenizer& operator =(const std::string &str)
	{
		end= str.c_str();
		length = FindToken(start , end);
		return *this;
	}
	
	std::string operator*() 
	{
		//remove white spaces in string
		while(isspace(*start))
		{
			start++;
			length--;
		}
		return std::string(start , length);
	}
	//++object
	Tokenizer& operator ++ () 
	{
		length = FindToken(start , end );
		return *this;
	}
	//object++
	const Tokenizer& operator ++ (int) 
	{
		length = FindToken(start , end );
		return *this;
	}

	bool operator!=(Tokenizer rhs) const 
	{
		return !operator==(rhs);
	}
	
	bool operator == (Tokenizer rhs) const 
	{
		bool bRet=false;
		if(*start == * rhs.start )
		  if( *start == EndOfString || start == rhs.start )
				bRet = true;
	
		  return bRet;
	}

	~Tokenizer(){};

private:
		
	TFinder FindToken;
	const char* start;
	const char* end;
	size_t length;
};

template <typename TFinder>
const char Tokenizer<TFinder>::EndOfString =(char)'\0';


template <typename TFinder>
const Tokenizer<TFinder> Tokenizer<TFinder>::EOS;




