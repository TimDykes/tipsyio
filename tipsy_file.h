#ifndef TIPSY_FILE_H
#define TIPSY_FILE_H

/*
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * 		Tim Dykes 
 *

	File: tipsy_file.h
	Purpose: Encapsulates a Tipsy file, with option to swap endian, and have 1 int
			 padding on the header. Both default to yes.

			Structures gas/dark/star_particle and header all sourced from tipsydefs.h 
			from tipsy tools foudn here http://www-hpcc.astro.washington.edu/tools/tipsy/tipsy.html

 */

// #include "common/debug"
#include <stdlib.h>
#include <math.h>
#include <fstream>

#ifdef USE_MPI
#include "mpi.h"
#endif
			
#define MAXDIM 3

/* Hack in a definition for ErrorMessage to work aroudn missing debug include */
#define ErrorMessage(...) {printf(__VA_ARGS__); exit(-1);}
/* 
	Structures for tipsy particle types & header
*/

struct gas_particle {
    float mass;
    float pos[MAXDIM];
    float vel[MAXDIM];
    float rho;
    float temp;
    float hsmooth;
    float metals ;
    float phi ;
} ;

struct dark_particle {
    float mass;
    float pos[MAXDIM];
    float vel[MAXDIM];
    float eps;
    float phi ;
} ;

struct star_particle {
    float mass;
    float pos[MAXDIM];
    float vel[MAXDIM];
    float metals ;
    float tform ;
    float eps;
    float phi;
} ;

struct header {
    double time ;
    int nbodies ;
    int ndim ;
    int nsph ;
    int ndark ;
    int nstar ;
    int pad ;
} ;

/*
	Abstraction for whole Tipsy file
*/

class TipsyFile{
public:
	//  Data
	header h;
	gas_particle* sph;
	dark_particle* dark;
	star_particle* star;

	// Extra
	const char* name;
#ifdef USE_MPI
	MPI_File src;
#else 
	std::ifstream src;
#endif
	bool swap_endian;
	bool header_read;
	bool open;
	long local_nsph;
	long local_ndark;
	long local_nstar;
#ifdef USE_MPI
	int  rank;
	long sph_start;
	long dark_start;
	long star_start;
#endif

	TipsyFile() { header_read = false; sph = NULL; dark = NULL; star = NULL; open = false;}

	TipsyFile(const char* filename, bool swap = true)
	{
		header_read = false; 
		sph = NULL; 
		dark = NULL; 
		star = NULL;
		open = false;
		open(filename, swap);
	}

	// Create a new tipsy file
	void create()
	{
		// Set up header, particles, etc
		// Also set 'swap endian'
	}

	void open(const char* filename, bool swap = true)
	{
		if(open)
			ErrorMessage("TipsyFile: File %s is already open!\n", filename);
#ifdef USE_MPI
		
#else
		src.open(filename, std::ios::binary);
		if(!src.is_open())
			ErrorMessage("TipsyFile: Cannot open file: %s\n", filename);
#endif

		swap_endian = swap;		
		name = filename;
	}

	void read_header(bool hasPad = true)
	{
		if(header_read)
			return;

		int pad = 0;
		if(!hasPad)
		{
			pad = sizeof(int);
		}
#ifdef USE_MPI
		if(rank==0)
		{
			// Read header on master
		}
#else		
		src.read((char*)&h, sizeof(header)-pad);
#endif
		byteswap(&h.time);
        	byteswap(&h.nbodies);
       		byteswap(&h.ndim);
        	byteswap(&h.nsph);
        	byteswap(&h.ndark);
        	byteswap(&h.nstar);

#ifdef USE_MPI
		// Broadcast header to everyone else
		// MPI_Bcast
#endif

        	header_read = true;
	}

	void report_header()
	{
#ifdef USE_MPI
	if(rank==0)
#endif
		if(header_read)
			printf("TipsyFile Name: %s\ntime: %f nbodies %d ndim %d\nngas: %d, ndark %d, nstar %d\n", \
					name, h.time, h.nbodies, h.ndim, h.nsph, h.ndark, h.nstar);
		else
			printf("TipsyFile: report_header(): havent read header yet\n");
	}

	void read_sph()
	{
#ifdef USE_MPI
		// Set local nsph to nsph/nranks
		local_nsph = h.nsph/nranks;
		start_sph = local_nsph * rank;
		local_nsph = std::min(h.nsph-start_nsph, local_nsph);
#else
		// Set local nsph to header nsph
		local_nsph = h.nsph
#endif
		if(local_nsph)
		{
			sph = (gas_particle*)malloc(local_nsph*sizeof(gas_particle));
			if(!sph) ErrorMessage("Rank %d could not allocate memory", rank);
		}
		// Read sph
		if(local_nsph)
		{
#ifdef USE_MPI
			// Seek to sph plus local sph start
			// MPI Read
#else
			// Seek to file sph start
			// ...
			src.read((char*)sph, local_nsph * sizeof(gas_particle));
#endif
			if(swap_endian)
			{
				gas_particle* pp = sph;
				for(int i = 0; i < local_nsph; i++, pp++)
				{
					for(unsigned j = 0; j < sizeof(gas_particle)/sizeof(float); j++)
						byteswap(&((float*)pp)[j]);
				}
			}
		}
	}

	void read_all(bool hasPad = true)
	{

		if(open)
			ErrorMessage("TipsyFile: read_all(): file is not open\n");

		if(!header_read)
			read_header(hasPad);

/*
		if( (h.nsph > 0 && sph == NULL) || (h.ndark > 0 && dark == NULL) || (h.nstar > 0 && star == NULL))
		{
			// Couldnt allocate, cleanup and quit
			close();
			ErrorMessage("Could not allocate memory for particles...\n");
		}


		// Read dark
		if(h.ndark)
		{
			src.read((char*)dark, h.ndark * sizeof(dark_particle));
			if(swap_endian)
			{
				dark_particle* pp = dark;
				for(int i = 0; i < h.ndark; i++, pp++)
				{
					for(unsigned j = 0; j < sizeof(dark_particle)/sizeof(float); j++)
						byteswap(&((float*)pp)[j]);
				}
			}
		}

		// Read star
		if(h.nstar)
		{
			src.read((char*)star, h.nstar * sizeof(star_particle));
			if(swap_endian)
			{
				star_particle* pp = star;
				for(int i = 0; i < h.nstar; i++, pp++)
				{
					for(unsigned j = 0; j < sizeof(star_particle)/sizeof(float); j++)
						byteswap(&((float*)pp)[j]);
				}
			}
		}
*/

#ifdef USE_MPI
		printf("TipsyFile: read file %s\nnbodies: %i\nnsph:    %i\nndark:   %i\nnstar:   %i\nswapped endian: %s\n", \
			    name, h.nbodies, h.nsph, h.ndark, h.nstar, (swap_endian) ? "yes" : "no");

		// Close file
#else
		printf("TipsyFile: read file %s\nnbodies: %i\nnsph:    %i\nndark:   %i\nnstar:   %i\nswapped endian: %s\n", \
			    name, h.nbodies, h.nsph, h.ndark, h.nstar, (swap_endian) ? "yes" : "no");

		if(src.is_open())
			src.close();
#endif
		
	}



	void close()
	{
		if(sph) free(sph);
		if(dark) free(dark);
		if(star) free(star);
		if(open)
			src.close();
	}

	void write(std::string name, bool hasPad = true)
	{
		// Output file
		std::ofstream out(name.c_str(), std::ios::binary);

		if(!out.is_open())
		{
			ErrorMessage("TipsyFile: write() could not open file %s for output.\n", name.c_str());
		}

		// Header struct includes padding, if we dont want padding dont write the final int.
		// We always swap endianness of header - but if we didnt also swap the endianness 
		// of the file, then we should swap back header

		int nsph = h.nsph;
		int ndark = h.ndark;
		int nstar = h.nstar;

		if(!swap_endian)
		{
			byteswap(&h.time);
	        byteswap(&h.nbodies);
	        byteswap(&h.ndim);
	        byteswap(&h.nsph);
	        byteswap(&h.ndark);
	        byteswap(&h.nstar);
		}
		out.write((char*)&h, sizeof(header) - (!hasPad ? sizeof(int) : 0) );

		// Write sph
		if(nsph > 0)
			out.write((char*)sph, sizeof(gas_particle) * nsph);

		// Write dark
		if(ndark > 0)
			out.write((char*)dark, sizeof(dark_particle) * ndark);

		// Write star
		if(nstar > 0)
			out.write((char*)star, sizeof(star_particle) * nstar);

		out.close();
	};

private:

template <typename T>
void byteswap(T* in)
{
    unsigned size = sizeof(T);
    for(unsigned i = 0; i < size/2; i++)
    {
        // Swap bytes
        ((char*)in)[i]      ^= ((char*)in)[size-i-1];
        ((char*)in)[size-i-1] ^= ((char*)in)[i];
        ((char*)in)[i]      ^= ((char*)in)[size-i-1];
    }
}

};

#endif
