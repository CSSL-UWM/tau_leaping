#include "functions.h"
#include <fstream>
#include <cassert>

using std::vector;
using std::ifstream;

#ifndef _WIN32

#include <fcntl.h>

bool FillRandoms(char* seeds, uint size)
{
	int ff = open("/dev/urandom", O_RDONLY);

	if(!ff)
	{
	   fprintf(stderr, "Can't open /dev/urandom\n");
	   return false;
	}
	
	if(size != read(ff, &seeds[0], size))
	{
		fprintf(stderr, "Can't read random seeds\n");
		return false;
	}

	close(ff);

	return true;
}

#else

#include <windows.h>

bool FillRandoms(char* seeds, unsigned int size)
{
	HMODULE hLib=LoadLibrary("ADVAPI32.DLL");
	if (!hLib) 
	{
		fprintf(stderr, "Can't load ADVAPI32.DLL library\n");
	    return false;
	}

	BOOLEAN (APIENTRY *pfn)(void*, ULONG) = 
		(BOOLEAN (APIENTRY *)(void*,ULONG))GetProcAddress(hLib,"SystemFunction036");
	if (!pfn) 
	{
		FreeLibrary(hLib);
		fprintf(stderr, "Can't find SystemFunction036 function addres in the ADVAPI32.DLL library\n");
	    return false;
	}

	if(!pfn(seeds,size)) 
	{
		FreeLibrary(hLib);
		fprintf(stderr, "Can't find SystemFunction036 function address in the ADVAPI32.DLL library\n");
	    return false;
	}
	FreeLibrary(hLib);

	return true;
}

#endif

void FileToMem(const char *Name, vector<char> &buffer)
{
	vector<char>().swap(buffer);
	ifstream f(Name, std::ios::in | std::ios::binary);

    if(!f)
    {
		char msg[256];
		sprintf(msg, "Can't open file %s", Name);
		throw std::runtime_error(msg);
    }
    ifstream::pos_type size = 0;

    if( f.seekg(0, std::ios::end))
    {
       size = f.tellg();
    }
    if( size && f.seekg(0, std::ios::beg) )
    {
       buffer.resize((int)size+1);
       f.read(&buffer[0], size);
	   *buffer.rbegin()=0;
    }    
    if((int) size==0)
	{

		char msg[256];
		sprintf(msg, "File %s is empty", Name);
		throw std::runtime_error(msg);
	}
}

bool FileExists(const char *fname)
{
	return ifstream(fname) != NULL;
}


