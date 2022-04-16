/**
    bd-anchors: Bidirectional String Anchors
    Copyright (C) 2022 Grigorios Loukides and Solon P. Pissis. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <chrono>
#include <ctime>
#include <random>
#include "grid.h"

typedef grid_point point;
typedef grid_query query;
 
#ifdef _USE_64
typedef int64_t INT;
#endif

#ifdef _USE_32
typedef int32_t INT;
#endif

using namespace std;


#include <boost/functional/hash.hpp>
using namespace boost;


#ifdef _USE_64
#include <divsufsort64.h>                                         // include header for suffix sort
#endif
#ifdef _USE_32
#include <divsufsort.h>                                       	  // include header for suffix sort
#endif

#include <sdsl/bit_vectors.hpp>                                   // include header for bit vectors
#include <sdsl/rmq_support.hpp>

using namespace sdsl;

/* Sorting comparison */
bool sort_sa(const pair<INT,INT> &a,const pair<INT,INT> &b)
{
       return a.first<b.first;
}

/* Kasai et al algorithm for LCP construction */
INT LCParray ( unsigned char * text, INT n, INT * SA, INT * ISA, INT * LCP )
{
        INT i=0, j=0;

        LCP[0] = 0;
        for ( i = 0; i < n; i++ ) // compute LCP[ISA[i]]
                if ( ISA[i] != 0 )
                {
                        if ( i == 0) j = 0;
                        else j = (LCP[ISA[i-1]] >= 2) ? LCP[ISA[i-1]]-1 : 0;
                        while ( text[i+j] == text[SA[ISA[i]-1]+j] )
                                j++;
                        LCP[ISA[i]] = j;
                }

        return ( 1 );
}


/* Booth's algorithm for computing the (leftmost) minimal string rotation */
INT minlexrot( string& X ) 
{    
	//string S = X;
    	//S.append(X.begin(), X.end());
    	//INT n = S.size();
    
    	//unsigned char * seq = (unsigned char *) S.c_str();
    	
    	X.append(X.begin(), X.end());
		INT n=X.size();
		unsigned char * seq = (unsigned char *)X.c_str();
	
	
    	
    	vector <INT> f (n, (INT)-1);
    
    	INT k = 0;
    	for (INT j = 1; j < n; ++j) 
	{
        	unsigned char sj = seq[j];
        	INT i = f[j - k - 1];
        	while (i != (INT)-1 && sj != seq[k + i + 1]) 
		{
            		if (sj < seq[k + i + 1])	k = j - i - 1;
            		i = f[i];
        	}

        	if (i==(INT)-1 && sj != seq[k + i + 1]) 
		{
            		if (sj < seq[k+i+1])	k = j;
            		f[j - k] = -1;
        	} 
		else 
            		f[j - k] = i + 1;
    	}

    	return k;
}

/* Computes the bd-anchors of a string of length n in O(wn) time */
void fast_anchors(string &whole_string, unordered_set<INT> &my_map, INT w)
{
	INT whole_string_len=whole_string.size();
   	INT start_pos = 0;
	while( (start_pos + w) <= whole_string_len )
  	{
		string string_in_win = whole_string.substr(start_pos,w);
		INT anchor_pos = start_pos + minlexrot(string_in_win);
	 	my_map.insert(anchor_pos);
     		start_pos = start_pos + 1; 
  	}
}

/* Constructs the right compacted trie given the anchors and the SA of the whole string */
void right_compacted_trie ( unordered_set<INT> &anchors, INT * SA, INT * LCP, INT n, INT * RSA, INT * RLCP, INT g )
{
	INT ii = 0; //this is the index over RSA[o..g-1] and the RLCP[0..g-1], where g is the number of anchors
	INT minLCP = n; //this is to maintain the minimum LCP in a range over LCP[i..j]: stores the LCP of SA[i] and SA[j]

	for( INT i = 0; i < n; i++ ) // in lex order
	{
		/* If the ith lex suffix is an anchor then add it to the compacted trie (encoded in arrays RSA and RLCP) */
		auto it = anchors.find( SA[i] );
		if( it != anchors.end() )
		{
			RSA[ii] = SA[i];		// store this suffix
			if ( ii == 0 )	RLCP[ii] = 0; 	// if it is the first time the RLCP = LCP = 0
			else
			{
				if ( SA[i-1] == RSA[ii-1] )     // if the immediately prior suffix was added
					RLCP[ii] = LCP[i];	// then the LCP value is the correct one for RLCP
				else
					RLCP[ii] = std::min(minLCP, LCP[i]);	//otherwise, we should take the minimum in the range
			}
  			//cout<<"RSA[i]: "<< RSA[ii] <<" RLCP[i]: "<<RLCP[ii]<<"\n"; getchar();
			minLCP = n; // set this to something high to get the _FIRST_ next minimum value, because a new range STARTS
			ii++;
		}
		else /* Do not add it but remember the minLCP seen so far in a range*/
		{
			if ( LCP[i] < minLCP )
				minLCP = LCP[i];
		}
	}
}

/* Constructs the left compacted trie given the anchors and the SA of the whole string */
void left_compacted_trie ( unordered_set<INT> &anchors, INT * SA, INT * LCP, INT n, INT * RSA, INT * RLCP, INT g )
{
	INT ii = 0;
	INT minLCP = n;
	for( INT i = 0; i < n; i++ ) // in lex order 
	{
		/* If the ith lex suffix is an anchor then add it to the compacted trie (encoded in arrays RSA and RLCP) */
		auto it = anchors.find( ( n - 1 ) - SA[i] );
		if( it != anchors.end() )
		{
			RSA[ii] = SA[i];		// store this suffix
			if ( ii == 0 )	RLCP[ii] = 0; 	// if it is the first time the RLCP = LCP = 0
			else
			{
				if ( SA[i-1] == RSA[ii-1] ) // if the immediately prior suffix was added
					RLCP[ii] = LCP[i];	//then the LCP value is the correct one for RLCP
				else
					RLCP[ii] = std::min(minLCP, LCP[i]);	//otherwise, we should take the minimum in the range
			}
  			
			//cout<<"LSA[i]: "<< RSA[ii] <<" LLCP[i]: "<< RLCP[ii]<<"\n"; getchar();
			minLCP = n; //set this to something high to get the FIRST next minimum value
			ii++;
		}
		else /* Do not add it but remember the minLCP seen so far in a range*/
		{
			if ( LCP[i] < minLCP )
				minLCP = LCP[i];
		}
	}
}


/* A fast RMQ data structure for Manber & Myers algorithm: construct the recursion tree over array arr; at each node (i,j) of the tree we maintain the minimum value in arr[i+1..j] */
INT fast_RMQ ( INT * arr, INT L, INT R, INT n, std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >> &rmq )
{
	if ( R - L > 1 ) // Internal nodes
	{
		INT M = (L+R)/2;
		INT a = fast_RMQ ( arr, L, M, n, rmq ); //Recurse on the left
		INT b = fast_RMQ ( arr, M, R, n, rmq ); //Recurse on the right
		INT value = std::min(a, b); 	     //This is the minimum of arr[L+1..M]

		pair<INT,INT> p(L+1, R);
		pair<pair<INT,INT>,INT> np(p, value);
		rmq.insert(np);
		return value;
	}
	else  // Base case: leaf nodes
	{
		INT value;
		if ( R >= n ) value = 0;
		else value = std::min (arr[L+1], arr[R]);
		pair<INT,INT> p(L+1, R);
		pair<pair<INT,INT>,INT> np(p, value);
		rmq.insert(np);
		return value;
	}
}

/* Computes the length of lcp of two suffixes of two strings */
INT lcp ( string & x, INT M, string & y, INT l )
{
	INT xx = x.size();
	if ( M >= xx ) return 0;
	INT yy = y.size();
	if ( l >= yy ) return 0;

	INT i = 0;
	while ( ( M + i < xx ) && ( l + i < yy ) )
	{
		if ( x[M+i] != y[l+i] )	break;
		i++;
	}
	return i;
}

/* Searching a list of strings using LCP from "Algorithms on Strings" by Crochemore et al. Algorithm takes O(m + log n), where n is the list size and m the length of pattern */
pair<INT,INT> pattern_matching ( string & w, string & a, INT * SA, INT * LCP, std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >> &rmq, INT n )
{
	INT m = w.size(); //length of pattern
	INT N = a.size(); //length of string 
	INT d = -1;
	INT ld = 0;
	INT f = n;
	INT lf = 0;
	
	pair<INT,INT> interval;
	
	while ( d + 1 < f )
	{
		INT i = (d + f)/2;
		std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >>::iterator it;

		/* lcp(i,f) */
		INT lcpif;
		it = rmq.find(make_pair(i+1, f));
		lcpif = it -> second;
		
		/* lcp(d,i) */
		INT lcpdi;
		it = rmq.find(make_pair(d+1, i));
		lcpdi = it -> second;

		if ( ( ld <= lcpif ) && ( lcpif < lf ) )
		{
			d = i;
			ld = lcpif;
		}
		else if ( ( ld <= lf ) && ( lf < lcpif ) ) 	f = i;
		else if ( ( lf <= lcpdi ) && ( lcpdi < ld ) )
		{
			f = i;
			lf = lcpdi;
		}	
		else if ( ( lf <= ld ) && ( ld < lcpdi ) )	d = i;
		else
		{
			INT l = std::max (ld, lf);
			l = l + lcp ( a, SA[i] + l, w, l );
			if ( l == m ) //lower bound is found, let's find the upper bound
		        {
				INT e = i;
				while ( d + 1 < e )
				{
					INT j = (d + e)/2;
					
					/* lcp(j,e) */
					INT lcpje;
					it = rmq.find(make_pair(j+1, e));
					lcpje = it -> second;

					if ( lcpje < m ) 	d = j;
					else 			e = j;
				}

				/* lcp(d,e) */
				INT lcpde;
				it = rmq.find(make_pair(d+1, e));
				lcpde = it -> second;

				if ( lcpde >= m )	d = std::max (d-1,( INT ) -1 );

				e = i;
				while ( e + 1 < f )
				{
					INT j = (e + f)/2;
					
					/* lcp(e,j) */
					INT lcpej;
					it = rmq.find(make_pair(e+1, j));
					lcpej = it -> second;

					if ( lcpej < m ) 	f = j;
					else 			e = j;
				}

				/* lcp(e,f) */
				INT lcpef;
				it = rmq.find(make_pair(e+1, f));
				lcpef = it -> second;

				if ( lcpef >= m )	f = std::min (f+1,n);

				interval.first = d + 1;
				interval.second = f - 1;
				return interval;
				
				
			}
			else if ( ( l == N - SA[i] ) || ( ( SA[i] + l < N ) && ( l != m ) && ( a[SA[i]+l] < w[l] ) ) )
			{
				d = i;
				ld = l;
			}
			else
			{
				f = i;
				lf = l;
			}	

		}
	}

	interval.first = d + 1;
	interval.second = f - 1;
	return interval;
}

/* Computes the length of lcs of two suffixes of two strings */
INT lcs ( string & x, INT M, string & y, INT l )
{
	if ( M < 0 ) return 0;
	INT yy = y.size();
	if ( l >= yy ) return 0;

	INT i = 0;
	while ( ( M - i >= 0 ) && ( l + i < yy ) )
	{
		if ( x[M-i] != y[l+i] )	break;
		i++;
	}
	return i;
}


/* Searching a list of strings using LCP from "Algorithms on Strings" by Crochemore et al. Algorithm takes O(m + log n), where n is the list size and m the length of pattern */
pair<INT,INT> rev_pattern_matching ( string & w, string & a, INT * SA, INT * LCP, std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >> &rmq, INT n )
{
	INT m = w.size(); //length of pattern
	INT N = a.size(); //length of string 
	INT d = -1;
	INT ld = 0;
	INT f = n;
	INT lf = 0;
	
	pair<INT,INT> interval;
	
	while ( d + 1 < f )
	{
		INT i = (d + f)/2;
		INT revSA = N - 1 - SA[i];
		std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >>::iterator it;

		/* lcp(i,f) */
		INT lcpif;
		it = rmq.find(make_pair(i+1, f));
		lcpif = it -> second;
		
		/* lcp(d,i) */
		INT lcpdi;
		it = rmq.find(make_pair(d+1, i));
		lcpdi = it -> second;

		if ( ( ld <= lcpif ) && ( lcpif < lf ) )
		{
			d = i;
			ld = lcpif;
		}
		else if ( ( ld <= lf ) && ( lf < lcpif ) ) 	f = i;
		else if ( ( lf <= lcpdi ) && ( lcpdi < ld ) )
		{
			f = i;
			lf = lcpdi;
		}	
		else if ( ( lf <= ld ) && ( ld < lcpdi ) )	d = i;
		else
		{
			INT l = std::max (ld, lf);
			l = l + lcs ( a, revSA - l, w, l );
			if ( l == m ) //lower bound is found, let's find the upper bound
		        {
				INT e = i;
				while ( d + 1 < e )
				{
					INT j = (d + e)/2;
					
					/* lcp(j,e) */
					INT lcpje;
					it = rmq.find(make_pair(j+1, e));
					lcpje = it -> second;

					if ( lcpje < m ) 	d = j;
					else 			e = j;
				}

				/* lcp(d,e) */
				INT lcpde;
				it = rmq.find(make_pair(d+1, e));
				lcpde = it -> second;

				if ( lcpde >= m )	d = std::max (d-1,( INT ) -1 );

				e = i;
				while ( e + 1 < f )
				{
					INT j = (e + f)/2;
					
					/* lcp(e,j) */
					INT lcpej;
					it = rmq.find(make_pair(e+1, j));
					lcpej = it -> second;

					if ( lcpej < m ) 	f = j;
					else 			e = j;
				}

				/* lcp(e,f) */
				INT lcpef;
				it = rmq.find(make_pair(e+1, f));
				lcpef = it -> second;

				if ( lcpef >= m )	f = std::min (f+1,n);

				interval.first = d + 1;
				interval.second = f - 1;
				return interval;
				
				
			}
			else if ( ( l == N - SA[i] ) || ( ( revSA - l >= 0 ) && ( l != m ) && ( a[revSA - l] < w[l] ) ) )
			{
				d = i;
				ld = l;
			}
			else
			{
				f = i;
				lf = l;
			}	

		}
	}

	interval.first = d + 1;
	interval.second = f - 1;
	return interval;
}

int main(int argc, char** argv)
{
	if( argc != 4 )
 	{
        	cout<<"Wrong arguments!\n";
		cout<<"./index <text_file> <ell> <pattern_file>\n";
 		exit(-1);
 	}

 	ifstream is;
 	is.open (argv[1], ios::in | ios::binary);
 
 	std::string str2(argv[2]);
 
 	ifstream is2; 
 	is2.open (argv[3], ios::in | ios::binary);
 
 	INT ell;
 	std::stringstream(str2)>>ell;
 
 	cout<<"The parameter ell is set to "<<ell<<endl;
 
	ifstream in_file(argv[1], ios::binary);
   	in_file.seekg(0, ios::end);
   	INT file_size = in_file.tellg();

  	vector<unsigned char> text;
  	char c = 0; 
	for (INT i = 1; i < file_size; i++)
	{
		is.read(reinterpret_cast<char*>(&c), 1); 
		text.push_back( (unsigned char) c );
	}
  	is.close();      
  
 	unordered_set<INT> text_anchors;
	string text_string(text.begin(), text.end());
  	
	if( ell < 1 || ell > text_string.size())
      	{
        	cout<<"w must be in [1, "<<text_string.size()<<"]"<<endl;
         	return -1;
      	}
	  
      	fast_anchors(text_string, text_anchors, ell);
	INT g = text_anchors.size();
	cout<<"The text is of length "<< text_string.size() << " and has "<<g<<" bd-anchors of order "<<ell<<endl;
	cout<<"The density is "<<(double) g / text_string.size()<<endl;

  	INT * SA;
  	INT * LCP;
  	INT * invSA;

  	INT n = text_string.size();
  	unsigned char * seq = ( unsigned char * ) text_string.c_str();

  	/* Compute the suffix array */
  	SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
  	if( ( SA == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
        	return ( 0 );
  	}
  	if( divsufsort64( seq, SA,  n ) != 0 )
  	{
  		fprintf(stderr, " Error: SA computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}
  
  	/*Compute the inverse SA array */
  	invSA = ( INT * ) calloc( n , sizeof( INT ) );
 	if( ( invSA == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
        	return ( 0 );
  	}
  	for ( INT i = 0; i < n; i ++ )
  	{
  		invSA [SA[i]] = i;
  	}

  	/* Compute the LCP array */
  	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );
  	if( ( LCP == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
        	return ( 0 );
  	}
  	if( LCParray( seq, n, SA, invSA, LCP ) != 1 )
  	{
		fprintf(stderr, " Error: LCP computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}
  
  	cout<<"SA and LCP constructed"<<endl;

  	/* Constructing right and left compacted tries */
  	INT * RSA;
  	INT * RLCP;

  	RSA = ( INT * ) malloc( ( g ) * sizeof( INT ) );
  	if( ( RSA == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for RSA.\n" );
        	return ( 0 );
  	}

  	RLCP = ( INT * ) malloc( ( g ) * sizeof( INT ) );
  	if( ( RLCP == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for RLCP.\n" );
        	return ( 0 );
  	}

  	right_compacted_trie ( text_anchors, SA, LCP, n, RSA, RLCP, g );  
  	cout<<"Right Compacted trie constructed "<<endl;
  
  	/* We reverse the string for the left direction and also overwrite all other DSs */
  	reverse(text_string.begin(), text_string.end());
  
  	if( divsufsort64( seq, SA,  n ) != 0 )
  	{
  		fprintf(stderr, " Error: SA computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}
  
  	for ( INT i = 0; i < n; i ++ )
  	{
  		invSA [SA[i]] = i;
  	}
  
  	if( LCParray( seq, n, SA, invSA, LCP ) != 1 )
  	{
		fprintf(stderr, " Error: LCP computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}

  	INT * LSA;
  	INT * LLCP;

  	LSA = ( INT * ) malloc( ( g ) * sizeof( INT ) );
  	if( ( LSA == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for LSA.\n" );
        	return ( 0 );
  	}
  	LLCP = ( INT * ) malloc( ( g ) * sizeof( INT ) );
  	if( ( LLCP == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for LLCP.\n" );
        	return ( 0 );
  	}

  	left_compacted_trie ( text_anchors, SA, LCP, n, LSA, LLCP, g );
  	cout<<"Left Compacted trie constructed"<<endl;
	

  	/* After constructing the tries these DSs over the whole string are not needed anymore, our data structure must be of size O(g) */
  	text_anchors.clear();
  	free ( SA );
  	free ( LCP );
  	free ( invSA );
  	cout<<"SA and LCP of the whole string freed"<<endl;
    
  	/* The following RMQ data structures are used for spelling pattern over the LSA and RSA */
  	std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >> lrmq;
  	std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >> rrmq;
  
  	cout<<"Left RMQ DS constructed "<<endl;
  	fast_RMQ ( LLCP, -1, g, g, lrmq ); //construction
  	cout<<"Right RMQ DS constructed "<<endl;
  	fast_RMQ ( RLCP, -1, g, g, rrmq ); //construction

	
  	//Forming coordinate points using LSA and RSA
  	vector<pair<INT,INT>> l;
	vector<pair<INT,INT>> r;
	for ( INT i = 0; i < g; i++ )
  	{
  		l.push_back( make_pair( n-1-LSA[i], i) );
  		r.push_back( make_pair( RSA[i], i ) );
 	}
 	
 	sort(l.begin(),l.end(),sort_sa);
	sort(r.begin(),r.end(),sort_sa);

  	std::vector<point> points;
  	
  	for ( INT i = 0; i < g; i++ )
  	{
 
		point to_insert;
		to_insert.row = l.at(i).second+1;
		to_insert.col = r.at(i).second+1;
		
		to_insert.level = 0;
		to_insert.label = l.at(i).first;
		points.push_back(to_insert); 

  	}
  	
  	
  	//constructing grid with bd-anchors
  	grid construct;
  	construct.build( points, 0 );
  	
	l.clear();
	r.clear();
	points.clear();
	
	cout<<"The range tree is constructed"<<endl;
  
  	cout<<"The whole index is constructed"<<endl;
 
  	reverse(text_string.begin(), text_string.end()); 				//I re-reverse to take the original string
	//string rev_text_string = string(text_string.rbegin(),text_string.rend());	//I store the reverse
  	
	vector<vector<unsigned char> > all_patterns;
  	vector<unsigned char> pattern;
  	c = 0;
  	while (is2.read(reinterpret_cast<char*>(&c), 1))     
  	{
      		if(c=='\n')
      		{
			if(pattern.empty())	break; 				
			all_patterns.push_back(pattern);
			pattern.clear();		
      		}
      		else	pattern.push_back((unsigned char)c);
  	}
  	is2.close();      
  	pattern.clear();

  	
	std::chrono::steady_clock::time_point  begin = std::chrono::steady_clock::now();
 	for(auto &it_pat : all_patterns)   
 	{
   		string pattern(it_pat.begin(), it_pat.end());
   		
		if ( pattern.size() < ell )
		{
			cout<<"Pattern skipped: its length is less than ell!\n";
			continue;
		}
		string first_window = pattern.substr(0, ell);
		INT j = minlexrot( first_window ); 
		string left_pattern = pattern.substr(0, j+1);
  		reverse(left_pattern.begin(), left_pattern.end());
  		//cout<<rev_text_string<<endl;
  		//cout<<left_pattern<<endl;
  		//pair<INT,INT> left_interval = pattern_matching ( left_pattern, rev_text_string, LSA, LLCP, lrmq, g );  
		pair<INT,INT> left_interval = rev_pattern_matching ( left_pattern, text_string, LSA, LLCP, lrmq, g );  
  		//cout<<"Left interval:"<< left_interval.first<<","<<left_interval.second<<endl;
  		//cout<<text_string<<endl;
		string right_pattern = pattern.substr(j, pattern.size()-j);
		//cout<<right_pattern<<endl;
		pair<INT,INT> right_interval = pattern_matching ( right_pattern, text_string, RSA, RLCP, rrmq, g );
		//cout<<"Right interval: "<<right_interval.first<<","<<right_interval.second<<endl;
   
   		
  		if ( left_interval.first <= left_interval.second  && right_interval.first <= right_interval.second )
  		{
  		
  			//Finding rectangle containing bd-anchors in grid
  			grid_query rectangle;
  			rectangle.row1 = left_interval.first+1;
  			rectangle.row2 = left_interval.second+1;
  			rectangle.col1 = right_interval.first+1;
  			rectangle.col2 = right_interval.second+1;
  		
			
			vector<long unsigned int> result;
			construct.search_2d(rectangle, result);
			for(INT i = 0; i<result.size(); i++)
				cout<<pattern<<" found at position "<<RSA[result.at(i)-1]-j<<" of the text"<<endl;
			//cout<<endl;
	
  		}
  		//else	cout<<"No occurrences found!\n";
 	}
	std::chrono::steady_clock::time_point  end = std::chrono::steady_clock::now();
	std::cout <<"Pattern matching took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
  	free ( RSA );
  	free ( RLCP );
  	free ( LSA );
  	free ( LLCP );

	return 0;
}
