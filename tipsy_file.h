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
	MPI_Offset offset;
	MPI_Status status;
#else 
	std::ifstream src;
#endif
	bool swap_endian;
	bool header_read;
	bool is_open;
	long local_nsph;
	long local_ndark;
	long local_nstar;
	int  rank;

#ifdef USE_MPI
	MPI_Comm comm;
	int comm_size;
	long sph_start;
	long dark_start;
	long star_start;
#endif

	TipsyFile() { header_read = false; sph = NULL; dark = NULL; star = NULL; is_open = false;}

#ifdef USE_MPI
	TipsyFile(const char* filename, MPI_Comm comm, bool swap = true)
#else
	TipsyFile(const char* filename, bool swap = true)
#endif
	{
		header_read = false; 
		sph = NULL; 
		dark = NULL; 
		star = NULL;
		is_open = false;
#ifdef USE_MPI
		open(filename, comm, swap);
#else
		open(filename, swap);
#endif
	}

	// Create a new tipsy file
	void create()
	{
		// Set up header, particles, etc
		// Also set 'swap endian'
	}
#ifdef USE_MPI
	void open(const char* filename, MPI_Comm c, bool swap = true)
#else
	void open(const char* filename, bool swap = true)
#endif
	{
		if(is_open)
			ErrorMessage("TipsyFile: File %s is already open!\n", filename);
#ifdef USE_MPI
		comm = c;
    	MPI_Comm_rank(comm, &rank);
    	MPI_Comm_size(comm, &comm_size);
		MPI_File_open(comm, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &src);
#else
		rank = 0;
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
		int header_size = sizeof(header)-pad;
		if(rank==0)
		{
			int read=0;
			MPI_File_read_at(src, 0, (unsigned char*)&h, header_size, MPI_BYTE, &status); 
			MPI_Get_count(&status, MPI_BYTE, &read);
			printf("TipsyFile: Rank 0 read %d bytes for header size %d\n", read, header_size);
		}
#else		
		src.read((char*)&h, sizeof(header)-pad);
#endif

#ifdef USE_MPI
		// Broadcast header to everyone else
		// MPI_Bcast
        MPI_Bcast( &h, header_size, MPI_Datatype MPI_BYTE, 0, comm);
#endif
		byteswap(&h.time);
        byteswap(&h.nbodies);
       	byteswap(&h.ndim);
        byteswap(&h.nsph);
        byteswap(&h.ndark);
        byteswap(&h.nstar);

		header_read = true;
	}

	void report_header()
	{
#ifdef USE_MPI
	if(rank==0)
#endif
		if(header_read)
			printf("TipsyFile Name: %s\ntime: %f nbodies %d ndim %d\nnsph: %d, ndark %d, nstar %d\n", \
					name, h.time, h.nbodies, h.ndim, h.nsph, h.ndark, h.nstar);
		else
			printf("TipsyFile: report_header(): havent read header yet\n");
	}

	void read_sph()
	{
		if(!is_open)
			ErrorMessage("TipsyFile: read_all(): file is not open\n");

		if(!header_read)
			read_header();

#ifdef USE_MPI
		// Set local nsph to nsph/nranks
		local_nsph = h.nsph/comm_size;
		sph_start = local_nsph * rank;
		local_nsph = std::min(h.nsph-sph_start, local_nsph);
#else
		// Set local nsph to header nsph
		local_nsph = h.nsph;
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

	void read_dark()
	{
		if(!is_open)
			ErrorMessage("TipsyFile: read_all(): file is not open\n");

		if(!header_read)
			read_header();

#ifdef USE_MPI
		// Set local ndark to ndark/nranks
		local_ndark = h.ndark/comm_size;
		dark_start = local_ndark * rank;
		local_ndark = std::min(h.ndark-dark_start, local_ndark);
#else
		// Set local ndark to header ndark
		local_ndark = h.ndark;
#endif
		if(local_ndark)
		{
			dark = (dark_particle*)malloc(local_ndark*sizeof(dark_particle));
			if(!dark) ErrorMessage("Rank %d could not allocate memory", rank);
		}
		// Read dark
		if(local_ndark)
		{
#ifdef USE_MPI
			// Seek to dark plus local dark start
			// MPI Read
#else
			// Seek to file dark start
			// ...
			src.read((char*)dark, local_ndark * sizeof(dark_particle));
#endif
			if(swap_endian)
			{
				dark_particle* pp = dark;
				for(int i = 0; i < local_ndark; i++, pp++)
				{
					for(unsigned j = 0; j < sizeof(dark_particle)/sizeof(float); j++)
						byteswap(&((float*)pp)[j]);
				}
			}
		}
	}

	void read_star()
	{
		if(!is_open)
			ErrorMessage("TipsyFile: read_all(): file is not open\n");

		if(!header_read)
			read_header();

#ifdef USE_MPI
		// Set local nstar to nstar/nranks
		local_nstar = h.nstar/comm_size;
		star_start = local_nstar * rank;
		local_nstar = std::min(h.nstar-star_start, local_nstar);
#else
		// Set local nstar to header nstar
		local_nstar = h.nstar;
#endif
		if(local_nstar)
		{
			star = (star_particle*)malloc(local_nstar*sizeof(star_particle));
			if(!star) ErrorMessage("Rank %d could not allocate memory", rank);
		}
		// Read star
		if(local_nstar)
		{
#ifdef USE_MPI
			// Seek to star plus local star start
			// MPI Read
#else
			// Seek to file star start
			// ...
			src.read((char*)star, local_nstar * sizeof(star_particle));
#endif
			if(swap_endian)
			{
				star_particle* pp = star;
				for(int i = 0; i < local_nstar; i++, pp++)
				{
					for(unsigned j = 0; j < sizeof(star_particle)/sizeof(float); j++)
						byteswap(&((float*)pp)[j]);
				}
			}
		}
	}

	void read_all(bool hasPad = true)
	{

		if(!is_open)
			ErrorMessage("TipsyFile: read_all(): file is not open\n");

		if(!header_read)
			read_header(hasPad);

		read_sph();
		read_dark();
		read_star();

#ifdef USE_MPI
		printf("TipsyFile: read file %s\nnbodies: %i\nnsph:    %i\nndark:   %i\nnstar:   %i\nlocal_nsph:    %i\nlocal_ndark:   %i\nlocal_nstar:   %i\nswapped endian: %s\n", \
			    name, h.nbodies, h.nsph, h.ndark, h.nstar, local_nsph, local_ndark, local_nstar, (swap_endian) ? "yes" : "no");

		// Close with MPI
		MPI_File_close(&src); 
#else
		printf("TipsyFile: read file %s\nnbodies: %i\nnsph:    %i\nndark:   %i\nnstar:   %i\nswapped endian: %s\n", \
			    name, h.nbodies, h.nsph, h.ndark, h.nstar, (swap_endian) ? "yes" : "no");

		if(src.is_open())
			src.close();
#endif
		is_open = false;
		
	}



	void close()
	{
		if(sph) free(sph);
		if(dark) free(dark);
		if(star) free(star);
		if(is_open)
		{
#ifdef USE_MPI
			// Close with MPI
			MPI_File_close(&src); 
#else
			src.close();
#endif
			is_open = false;
		}
	}

	void write(std::string name, bool hasPad = true)
	{
#ifdef USE_MPI
		ErrorMessage("write() not supported in parallel reader yet");
#else
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
#endif
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
