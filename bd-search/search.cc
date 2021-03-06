/**
    bd-anchors: Bidirectional String Anchors
    Copyright (C) 2021 Grigorios Loukides and Solon P. Pissis. 

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

#define DEBUG false

#include <unordered_set>
#include <boost/functional/hash.hpp>
#include "utils.h"
#include "edlib.h"

#ifdef _USE_64
#include <divsufsort64.h>                                         // include header for suffix sort
#endif
#ifdef _USE_32
#include <divsufsort.h>                                         // include header for suffix sort
#endif

#include <sdsl/bit_vectors.hpp>                                   // include header for bit vectors
#include <sdsl/rmq_support.hpp>

using namespace sdsl;
using namespace std;
using namespace boost;

class mycomparison
{
 
  public:
 
  bool operator() (const pair<INT,INT>& lhs, const pair<INT,INT>&rhs) const
  {
     return (lhs.first<rhs.first);
  }
};

/* Kasai et al algorithm for O(n)-time LCP construction */
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


/* Booth's O(n)-time algorithm -- slightly adapted for efficiency */
INT minlexrot( string &X, INT *f, INT n)
{
        INT n_d = n<<1;
        for(INT i = 0; i < n_d; ++i)    f[i] = (INT) -1;

        INT k = 0;
        for (INT j = 1; j < n_d; ++j)
        {
                unsigned char sj = X[j%n];
                INT i = f[j - k - 1];
                while (i != (INT)-1 && sj != X[(k + i + 1)%n])
                {
                        if (sj < X[(k + i + 1)%n])        k = j - i - 1;
                        i = f[i];
                }

                if (i == (INT) - 1 && sj != X[(k + i + 1)%n])
                {
                        if (sj < X[(k+i+1)%n])    k = j;
                        f[j - k] = -1;
                }
                else
                        f[j - k] = i + 1;
        }
        return k;
}

/* Computes the bd-anchors of a string of length n in O(ell n) time */
void fast_anchors(string &whole_string, unordered_set<INT> &my_map, INT ell)
{
	INT whole_string_len=whole_string.size();
   	INT start_pos = 0;

    	INT *f = new INT[ell<<1];
	while( (start_pos + ell) <= whole_string_len )
  	{
		string string_in_win = whole_string.substr(start_pos, ell);
		INT anchor_pos = start_pos + minlexrot(string_in_win,f, ell);
	 	my_map.insert(anchor_pos);
     		start_pos = start_pos + 1;
  	}
	delete []f;
}

/* Constructs the right compacted trie given the anchors and the SA of the whole string in O(n) time */
void right_compacted_trie ( unordered_map<INT,INT> &anchors, INT * SA, INT * LCP, INT n, INT * RSA, INT * RLCP, unordered_map<INT,INT> &invRSA, INT g )
{
	INT ii = 0; 
	INT minLCP = n;

	for( INT i = 0; i < n; i++ ) 
	{
		
		auto it = anchors.find( SA[i] );
		if( it != anchors.end() )
		{
			RSA[ii] = SA[i];		
			invRSA[RSA[ ii ]] = ii;
			if ( ii == 0 )	RLCP[ii] = 0; 
			else
			{
				if ( SA[i-1] == RSA[ii-1] ) 
					RLCP[ii] = LCP[i];	
				else
					RLCP[ii] = std::min(minLCP, LCP[i]);
			}
  			
			minLCP = n; 
			ii++;
		}
		else 
		{
			if ( LCP[i] < minLCP )
				minLCP = LCP[i];
		}
	}
}

/* Constructs the left compacted trie given the anchors and the SA of the whole string in O(n) time */
void left_compacted_trie ( unordered_map<INT,INT> &anchors, INT * SA, INT * LCP, INT n, INT * RSA, INT * RLCP, unordered_map<INT,INT> &invRSA, INT g )
{
	INT ii = 0;
	INT minLCP = n;
	for( INT i = 0; i < n; i++ ) 
	{
		
		auto it = anchors.find( ( n - 1 ) - SA[i] );
		if( it != anchors.end() )
		{
			RSA[ii] = SA[i];		
			invRSA[(n-1) - RSA[ii]] = ii;	
			if ( ii == 0 )	RLCP[ii] = 0; 	
			else
			{
				if ( SA[i-1] == RSA[ii-1] ) 
					RLCP[ii] = LCP[i];	
				else
					RLCP[ii] = std::min(minLCP, LCP[i]);	
			}
  			
			
			minLCP = n; 
			ii++;
		}
		else 
		{
			if ( LCP[i] < minLCP )
				minLCP = LCP[i];
		}
	}
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
void rev_pattern_matching ( string & w, string & a, INT * SA, INT * LCP, std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >> &rmq, INT n, INT &pair1, INT &pair2 )
{
	INT m = w.size(); 
	INT N = a.size(); 
	INT d = -1;
	INT ld = 0;
	INT f = n;
	INT lf = 0;
		
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

				pair1=d+1;
				pair2=f-1;
				return;
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
	
	pair1=d+1;
	pair2=f-1;
	return;
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
	else  // Base case: lead nodes
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
void pattern_matching ( string & w, string & a, INT * SA, INT * LCP, std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >> &rmq, INT n, INT &pair1, INT &pair2 )
{
	INT m = w.size(); //length of pattern
	INT N = a.size(); //length of string 
	INT d = -1;
	INT ld = 0;
	INT f = n;
	INT lf = 0;
		
	
	while ( d + 1 < f )
	{
		INT i = (d + f)/2;
		std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >>::iterator it;

		/* lcp(i,f) */
		INT lcpif;
		it = rmq.find(make_pair(i+1, f));
		if ( it == rmq.end() )
		{
			cout<<"1. RMQ failure! Not able to find"<<i+1<<","<<f<<endl;
			break;
		}
		else lcpif = it -> second;
		
		/* lcp(d,i) */
		INT lcpdi;
		it = rmq.find(make_pair(d+1, i));
		if ( it == rmq.end() )
		{
			cout<<"2. RMQ failure! Not able to find"<<d+1<<","<<i<<endl;
			break;
		}
		else lcpdi = it -> second;

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
					if ( it == rmq.end() )
					{
						cout<<"3. RMQ failure! Not able to find"<<j+1<<","<<e<<endl;
						break;
					}
					else lcpje = it -> second;

					if ( lcpje < m ) 	d = j;
					else 			e = j;
				}

				/* lcp(d,e) */
				INT lcpde;
				it = rmq.find(make_pair(d+1, e));
				if ( it == rmq.end() )
				{
					cout<<"4. RMQ failure! Not able to find"<<d+1<<","<<e<<endl;
					break;
				}
				else lcpde = it -> second;

				if ( lcpde >= m )	d = std::max (d-1,( INT ) -1 );

				e = i;
				while ( e + 1 < f )
				{
					INT j = (e + f)/2;
					
					/* lcp(e,j) */
					INT lcpej;
					it = rmq.find(make_pair(e+1, j));
					if ( it == rmq.end() )
					{
						cout<<"5. RMQ failure! Not able to find"<<e+1<<","<<j<<endl;
						break;
					}
					else lcpej = it -> second;

					if ( lcpej < m ) 	f = j;
					else 			e = j;
				}

				/* lcp(e,f) */
				INT lcpef;
				it = rmq.find(make_pair(e+1, f));
				if ( it == rmq.end() )
				{
					cout<<"6. RMQ failure! Not able to find"<<e+1<<","<<f<<endl;
					break;
				}
				else lcpef = it -> second;

				if ( lcpef >= m )	f = std::min (f+1,n);

				pair1 = d + 1;
				pair2 = f - 1;
				return;
				
				
			}
			else if ( ( l == N - SA[i] ) || ( ( SA[i]+l < N ) && ( l != m ) && ( a[SA[i]+l] < w[l] ) ) )
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

	pair1 = d + 1;
	pair2 = f - 1;
	return;
}

/* Computing LIS in O(n log n) time -- function taken from http://comscigate.com/Books/contests/icpc.pdf */
struct IndexValueCompare
{
	inline bool operator() (const pair<INT, INT> &one, const pair<INT, INT> &another)
	{
        	return one.second < another.second;
    	}
};

vector<pair<INT,INT>> fast_lis(const vector<pair<INT, INT> > &sequence)
{
	vector<INT> parent(sequence.size());
    	set<pair<INT, INT>, IndexValueCompare> s;
    	for(INT i = 0; i < sequence.size(); ++i)
	{
        	pair<INT, INT> iv(i, sequence[i].second);
        	if(i == 0)
		{
            		s.insert(iv);
            		continue;
        	}
        	auto index = s.lower_bound(iv);
        	if(index != s.end())
		{
            		if(sequence[i].second < sequence[index->first].second)
			{
                		if(index != s.begin()) 
				{
                    			parent[i] = (--index)->first;
                    			index++;
                		}
                		s.erase(index);
            		}
        	} 
		else
            		parent[i] = s.rbegin()->first;
        	s.insert(iv);
    }
    vector<pair<INT,INT> > result(s.size());
    int index = s.rbegin()->first;
    for(auto iter = s.rbegin(); iter != s.rend(); index = parent[index], ++iter)	result[distance(iter, s.rend()) - 1] = sequence[index];
    return result;
}

/* Pattern matching algorithm based on the bd-anchors; that is Step 1 of the main algorithm */
//Note:THE ARGUMENT whole_string is reversed when this function is called
void fast_anchors_pattern_end(string& whole_string, string &pattern, INT ell, INT * LSA, INT * LLCP, std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >> &lrmq, INT * RSA, INT * RLCP, 
std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >> &rrmq, INT g,  vector<pair<INT,INT> > &pattern_answers, unordered_map<INT,INT> &anchors_and_window_ends_in_pattern,  unordered_map<INT,INT> &invRSA)
{    
	INT pattern_len = pattern.size();
   	INT whole_string_len = whole_string.size();
   	INT start_pos = 0;     
   
   	INT *f = new INT[ell<<1];

   	//start_pos is in pattern
   	while(( start_pos + ell) <= pattern_len )
  	{    
     		string string_in_win = pattern.substr(start_pos,ell).c_str();  
		//anchor_pos is in pattern
	  	//leftmost 
     		INT anchor_pos = start_pos + minlexrot(string_in_win, f, ell);
	
	 	auto it2 = anchors_and_window_ends_in_pattern.find(anchor_pos);
	 	if(it2 == anchors_and_window_ends_in_pattern.end())
	 	{
			anchors_and_window_ends_in_pattern.insert(make_pair(anchor_pos,start_pos+ell-1));
	 	}
	 	else
	 	{
			start_pos++;
		 	continue;		 
	 	}
	 
		if ( ell + start_pos - anchor_pos >= anchor_pos-start_pos ) //if the right part is bigger than the left part then search the right part to get a smaller interval on RSA
		{
	
  			string right_pattern = pattern.substr(anchor_pos,start_pos+ell-1-anchor_pos+1);
			INT pair1=0,pair2=0;
			pattern_matching ( right_pattern, whole_string, RSA, RLCP, rrmq, g, pair1, pair2 );
  			
			if(pair1>pair2)
			{
				start_pos = start_pos + 1; 
				continue;
			}
			INT anchor_pos_to_return=anchor_pos;
			for(INT i = pair1; i <= pair2; i++ ) //this can be an interval of size 10^10 and only one occurrence is valid.
			{
				
				INT index = RSA[i];				
				anchor_pos=anchor_pos_to_return;
				index--;	anchor_pos--;
				while ( ( anchor_pos >= start_pos ) && ( index >= 0 ) && ( whole_string[index] == pattern[anchor_pos] ) )
				{
					index--; anchor_pos--;
				}
				if(anchor_pos==start_pos-1)
				{
					
					pair<INT,INT> my_pair(anchor_pos_to_return,i);
					pattern_answers.push_back(my_pair);
				}
				
			}
		}
		else //othewrise search the left part to get a smaller interval on LSA
		{			
			string left_pattern = pattern.substr(start_pos, anchor_pos-start_pos+1);
			reverse(left_pattern.begin(), left_pattern.end());				
			
			INT rev_pair1=0,rev_pair2=0;
			rev_pattern_matching ( left_pattern, whole_string, LSA, LLCP, lrmq, g, rev_pair1,rev_pair2 );  
	
			if(rev_pair1 > rev_pair2)
			{	
				start_pos = start_pos + 1; 
				continue;
			}
			
			INT anchor_pos_to_return=anchor_pos;
			for(INT i = rev_pair1; i <= rev_pair2; i++ ) //this can be an interval of size 10^10 and only one occurrence is valid.
			{
				INT index = whole_string_len-1-LSA[i];	
				INT index_to_return=index;
				anchor_pos=anchor_pos_to_return;
				
					index++; anchor_pos++;
				
					while ( ( anchor_pos < start_pos+ell ) && ( index < whole_string_len ) && ( whole_string[index] == pattern[anchor_pos] ) )
					{
						anchor_pos++;						
						index++;
					
					}
					
					if(anchor_pos==start_pos+ell)
					{
						pair<INT,INT> my_pair(anchor_pos_to_return,invRSA[whole_string_len-1-LSA[i]]);			
						pattern_answers.push_back(my_pair);
					}
				
				
			}
		}				
	
     		start_pos = start_pos + 1; 
  	}
  
  	delete []f;
}

int main(int argc, char** argv)
{
	if(argc!=7)
 	{
        	cout<<"Wrong arguments!\n";
 		cout<<"./search <dictionary> <ell> <pattern> <K> <tau> <delta>\n";
 		exit(-1);
 	}
 
	 ifstream is;
 	is.open (argv[1], ios::in | ios::binary);
 
 	INT ell;
 	std::stringstream(argv[2])>>ell;
 
 	ifstream is2; 
 
 	is2.open (argv[3], ios::binary);
 
 	INT K;
 	stringstream(argv[4])>>K;
 
 	INT tau;
 	stringstream(argv[5])>>tau;

 	INT delta;
 	stringstream(argv[6])>>delta; 
 
 	cout<<"Parameters: ell = "<<ell<<", K = "<<K<<", tau = "<<tau<<", delta = "<<delta<<endl;
 	
 
  	vector<pair<INT, vector<unsigned char> > > all_strings;
  	vector<unsigned char> input_seq_vec;
  	INT id=0;
  	char c=0; 
	while (is.read(reinterpret_cast<char*>(&c), 1))     
	{
		if(c=='\n')
	    {
			if(input_seq_vec.empty())
				break; 
			vector<unsigned char> v(input_seq_vec);
			pair<INT, vector<unsigned char> > p(id, v);
			all_strings.push_back(p);
			input_seq_vec.clear();
			id++;
	    }
	    else
	    {
			input_seq_vec.push_back((unsigned char)c);
	    }
	}
	is.close();      
  
   
	INT *num_anchors_per_string=(INT *)malloc(id*sizeof(INT));
  	//key is the anchor plus offset and value string_id 
  	// to each anchor 
 	unordered_map<INT, INT> global_anchors;
  	INT total_length_of_previous_strings=0;

	//stores the total length of all previous strings with ids 0,1,...,m-1 in position m
	vector<INT> vec_total_length_of_previous_strings;
  	
	//for each string in the database
  	for(INT j=0; j<id ;++j)
  	{
		unordered_set<INT> anchors;
	  	string current_string(all_strings[j].second.begin(),all_strings[j].second.end()); 
	  	fast_anchors(current_string, anchors, ell);
	
	  	for(auto & it : anchors)
	  	{
			INT anchor = it + total_length_of_previous_strings; 
		  	//we add anchor and string_id
		  	global_anchors.insert(make_pair(anchor,j));
		}
	
	  	//how many anchors each string has
	  	num_anchors_per_string[j]=anchors.size();
	  
		vec_total_length_of_previous_strings.push_back(total_length_of_previous_strings);

	  	total_length_of_previous_strings +=current_string.size();		
  	}
  
  	unsigned char *my_whole_string=(unsigned char*) malloc(total_length_of_previous_strings*sizeof(INT));
  
  	for(INT i=0;i<total_length_of_previous_strings;)
  	{
		for(INT j=0;j<id;++j)  
		{
			string current_string(all_strings[j].second.begin(),all_strings[j].second.end()); 
			for(INT k=0;k<current_string.size();++k)
				my_whole_string[i++]=current_string[k];  
		}
  	}
  	string whole_string( reinterpret_cast<char const*>(my_whole_string), total_length_of_previous_strings ) ;

  	free(my_whole_string);
  
  	INT g=global_anchors.size();
  
  	INT * SA;
  	INT * LCP;
  	INT * invSA;

  	INT n = whole_string.size();
  	unsigned char * seq = ( unsigned char * ) whole_string.c_str();

  	/* Compute the suffix array */
  	SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
  	if( ( SA == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
        	return ( 0 );
  	}
	
	#ifdef _USE_64
  	if( divsufsort64( seq, SA,  n ) != 0 )
  	{
  		fprintf(stderr, " Error: SA computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}
	#endif

	#ifdef _USE_32
  	if( divsufsort( seq, SA,  n ) != 0 )
  	{
  		fprintf(stderr, " Error: SA computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}
	#endif
	
  
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
  
  	/* Constructing right and left compacted tries */
  	INT * RSA;
  	INT * RLCP;
  	unordered_map<INT,INT> invRSA;

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

  	right_compacted_trie ( global_anchors, SA, LCP, n, RSA, RLCP, invRSA, g );  
  
  	//To get rid of global_anchors
  	//I keep ids. ids[p] = string_id
  	INT *ids=(INT *)malloc((g)*sizeof(INT));
  
  	for(INT i=0;i<g;++i)	ids[i]=global_anchors[RSA[i]];
  
  	/* We reverse the string for the left direction and also overwrite all other DSs */
  	reverse(whole_string.begin(), whole_string.end());
  
  	#ifdef _USE_64
  	if( divsufsort64( seq, SA,  n ) != 0 )
  	{
  		fprintf(stderr, " Error: SA computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}
	#endif
	#ifdef _USE_32
  	if( divsufsort( seq, SA,  n ) != 0 )
  	{
  		fprintf(stderr, " Error: SA computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}
	#endif
	
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
  	unordered_map<INT,INT> invLSA;

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
	left_compacted_trie ( global_anchors, SA, LCP, n, LSA, LLCP, invLSA, g );
  	
	global_anchors.clear();
  
  	/* After constructing the tries these DSs over the whole string are not needed anymore, our data structure must be of size O(g) */
  
  	free ( SA );
  	free ( LCP );
  	free ( invSA );
  
  	//After this whole_string will be ok (left to right)
  	reverse(whole_string.begin(), whole_string.end());
 
    
  	/* The following RMQ data structures are used for spelling pattern over the LSA and RSA */
  	std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >> lrmq;
  	std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >> rrmq;
  
	fast_RMQ ( LLCP, -1, g, g, lrmq ); 
  	fast_RMQ ( RLCP, -1, g, g, rrmq ); 
	
  	invLSA.clear();
   
  	vector<vector<unsigned char> > all_patterns;
  
  	vector<unsigned char> new_vec;
  	unsigned char d=0;
 	unsigned char d2=0;
  
	while (is2.read(reinterpret_cast<char*>(&d), sizeof(d)))     
  	{
		if (!is2)
			cerr << "Read failed." << endl;
	
      		if(d=='\n')
      		{
			if(new_vec.empty())
				break;
			vector<unsigned char> v;
			for(auto itv=new_vec.begin();itv!=new_vec.end();++itv)
			{
				v.push_back(*itv);
			}		
			all_patterns.push_back(v);
			new_vec.clear();		
      		}
      		else
      		{
        		new_vec.push_back(d);		
      		}
  	}
  	is2.close();      
  
    
 	auto total_duration=0,total_duration2=0;
 
 	INT pruned_with_tau=0;  
 
 	INT cluster_id=0; 

 	for(auto &it_pat : all_patterns)   
 	{
		INT sum_of_edit_distance_of_answers=0; //for every query, we will sum up the edit distance between its answers and itself
		 /*                       QUERYING                               */
   
   		string pattern(it_pat.begin(),it_pat.end());
         		
   
   		cluster_id++;   
   
		/*  Estimate identity matching heuristic */
	
		unordered_map<INT, INT> anchors_and_window_ends_in_pattern;
		vector<pair<INT,INT> > pattern_answers2;

		//Consruct the input to LIS for each string - key is string_ID 
 		unordered_map<INT, vector<pair<INT,INT> > > all_LIS2;

    		priority_queue<pair<double, INT> > pq_est;	
	  
		auto start_time = chrono::steady_clock::now();
	
		/* Step 1: Find all hits and arrange them per database string */
		fast_anchors_pattern_end(whole_string, pattern, ell, LSA, LLCP, lrmq, RSA, RLCP, rrmq, g, pattern_answers2, anchors_and_window_ends_in_pattern, invRSA);     
 
    		total_duration2 += chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now()-start_time).count();	
	
		
				
  		for(auto &it : pattern_answers2)
  		{
			
			
			//pattern window end obtained from pattern's hit 
			INT pat_window_end=anchors_and_window_ends_in_pattern[it.first];
			INT string_id;
	  		string_id = ids[it.second]; 
		
			//stores last position of previous string (with coords of whole string)
			INT tlps=vec_total_length_of_previous_strings[string_id];
		
	  		INT window_end_of_string_hit=RSA[it.second]+(pat_window_end-it.first); //hit of string + offset (which is the same as the offset in the pattern window) - length of all previous strings
		
			//if string window end > length of current string, then hit is bogus, ignore. 
			if(window_end_of_string_hit-tlps<ell)
			{		
				continue;
			}
			
			auto it3 = all_LIS2.find(string_id);
	  
	    		if(it3 != all_LIS2.end())
	    		{
					
					all_LIS2[string_id].push_back(make_pair(it.first, RSA[it.second]));
				}		  
	    		else
	    		{
					vector<pair<INT,INT> > v;
					v.push_back(make_pair(it.first,RSA[it.second]));
					all_LIS2.insert(make_pair(string_id,v));		  		  
	    		}	  
  		}  
	    

		//Store the output of LIS for each string - key is string_ID 
		unordered_map<INT, vector<pair<INT,INT> > > string_id_LIS;
		
		/* Step 2: Construct the LISs per database string and estimate the identiy score */
   		for(auto &it : all_LIS2)
   		{				 	  

			/*cout<<"HITS:\n";
			for(auto &itx : it.second)
				cout<<"<"<<itx.first<<","<<itx.second<<"> ";
			cout<<endl;
			*/
			//Now we prune <tau
			if(it.second.size()<tau)
			{
				continue;
			}
			
			
			vector<pair<INT,INT> > ret = fast_lis(it.second); 		
			
			
			//remove bogus HITS
			vector<pair<INT, INT> > pruned_LIS;
			for(vector<pair<INT,INT> >::iterator itr=ret.begin();itr!=ret.end();++itr)
			{
				pruned_LIS.push_back(*itr);
				
				auto itr_next=itr;
				itr_next++;
				
				//if next is not the end
				while(itr_next!=ret.end())
				{
					
					if((*itr).first==(*itr_next).first)
					{
						itr++;
						itr_next++;
					}	
					else
					{					
						break;
					}
				}
			}
			string_id_LIS.insert(make_pair(it.first,pruned_LIS));
			
			
			vector<INT> e_is;
			
			
			for(auto &it2 : pruned_LIS)
			{
				e_is.push_back(anchors_and_window_ends_in_pattern[it2.first]);
			
			}
			
						
			INT score=ell;
					
	   		for(vector<INT>::iterator it2=e_is.begin();it2!=e_is.end();++it2)
			{
				
				auto it_next=it2;
				it_next++;
				if(it_next!=e_is.end())
				{
					if(*it2 + ell <= *it_next)	
					{
						score += ell;
			
					}
					else
					{
						score += *it_next - *it2;
			
					}

				}							
			}
			string cur_string(all_strings[it.first].second.begin(),all_strings[it.first].second.end());
			
			
			//double score_prod=-1.0*score;
			double score_prod=score;
			//Add everything into priority queue
			pq_est.push(make_pair(score_prod,it.first));
					
   		}
		
		INT cnt=0;
		double Kth_score=delta; //score of the K-th Initialize with delta so that the first K are added into good_scores
		vector<INT> good_scores;//stores string-ids of the strings with scores >= Kth_score - delta
		
		vector<pair<INT,INT> > score_id;
		
		while(!pq_est.empty())
		{
			
			pair<double,INT> p=pq_est.top();
			
			
			cnt++;
			if(cnt==K)
			{
				Kth_score=p.first;
			}
			
			if(p.first>=Kth_score-delta)
			{
				good_scores.push_back(p.second);		
			}
			else
				break;
			
			pq_est.pop();			
		}	
		
		
		/* Step 3: Mind the gaps! Estimate an upper bound on the edit distance */
		
		
		EdlibAlignResult fed;
		
		priority_queue<pair<double, INT>,vector<pair<INT,INT> >, mycomparison > pq_ub;
			
		for(auto &it : good_scores)
   		{				 	  

		
			vector<INT> e_is;

			//stores the window_ends for string (based on the current string)
			vector<INT> window_ends_string; 

			//stores the total length of previous strings 
			INT tlps=vec_total_length_of_previous_strings[it];

			
			for(auto &it2 : string_id_LIS[it])
			{
				e_is.push_back(anchors_and_window_ends_in_pattern[it2.first]);
				window_ends_string.push_back(it2.second+anchors_and_window_ends_in_pattern[it2.first]-it2.first-tlps);				
			}
			
			string cur_string(all_strings[it].second.begin(),all_strings[it].second.end());	  		
						
			INT ub=0;
			
						
	   		for(vector<INT>::iterator it2=e_is.begin(),it_string_end=window_ends_string.begin();it2!=e_is.end() && it_string_end!=window_ends_string.end();++it2,++it_string_end)
			{
							
				auto it_next=it2;
				it_next++;
				
					
				if(it2==e_is.begin())
				{
								
					INT y1=*it2;				//pat window end 	
					INT x1=y1-ell+1; //pat window start
					
					
					INT win_string_end=(*it_string_end);  //string window end 
					INT win_string_start=win_string_end-ell+1; //string window start
				
							
			
					if(y1==win_string_end && x1>0) //synced  and no overlap
					{
						
						 
						 
								//cout<<pattern.substr(0,x1)<<" "<<cur_string.substr(0,win_string_start)<<" x1="<<x1<<" y1="<<y1<<" win_string_start="<<win_string_start<<endl;														
								fed=edlibAlign(pattern.c_str(),x1,cur_string.c_str(),win_string_start, edlibDefaultAlignConfig());
								ub+=fed.editDistance;								
								//cout<<"1a."<<fed.editDistance;
						 
						
					}
					else
					{						
						if(x1==0)
						{
							//cout<<cur_string.substr(0,win_string_start)<<endl;
														
							ub+=win_string_start;
						}
						else if(win_string_start==0)
						{
							//cout<<pattern.substr(0,x1)<<endl;
							ub+=x1;														
						}
						else
						{
							//cout<<pattern.substr(0,x1)<<" "<<cur_string.substr(0,win_string_start)<<" x1="<<x1<<" y1="<<y1<<" win_string_start="<<win_string_start<<endl;							
							
							fed=edlibAlign(pattern.c_str(),x1,cur_string.c_str(),win_string_start, edlibDefaultAlignConfig());
							ub+=fed.editDistance;															
							//cout<<"2c. "<<fed.editDistance<<endl;
						}
					}
				}
				
				if(it_next!=e_is.end())  
				{
					INT y1=*it2;       //pat window end
					INT x1=y1-ell+1; //pat window start
                    INT y2=*it_next; //next pat window end 
					INT x2=y2-ell+1; //next pat window start
					
                                        
                    INT win_string_end=(*it_string_end);
                    INT win_string_start=(*it_string_end)-ell+1;
					
					auto it_next_string=it_string_end;
					it_next_string++; 
					INT next_win_string_end=(*it_next_string);
					INT next_win_string_start=(*it_next_string)-ell+1;					
					
					if(y2-y1==next_win_string_end-win_string_end) //sync
					{	
						
						if((ell<y2-y1) && (ell<next_win_string_end-win_string_end))       //no overlap 
						{
							//cout<<"3a "<<pattern.substr(y1+1,(x2-1-(y1+1)+1))<<" "<<cur_string.substr(win_string_end+1,(next_win_string_start-1-(win_string_end+1)+1))<<endl;
							//if no overlap then the windows of string won't have overlap as they are synced
							
							fed=edlibAlign(&(pattern.c_str())[y1+1],x2-(y1+1),&(cur_string.c_str())[win_string_end+1],next_win_string_start-(win_string_end+1), edlibDefaultAlignConfig());
							ub+=fed.editDistance;															
							//cout<<"3a. "<<fed.editDistance;
						}
						//else they have overlap do nothing
						
						
							
					}
					else //no sync
					{
						
						if((ell<y2-y1) && (ell<next_win_string_end-win_string_end) ) //none of them overlap
						{
							//cout<<"3b "<<pattern.substr(y1+1,(x2-1-(y1+1)+1))<<" "<<cur_string.substr(win_string_end+1,next_win_string_start-1-(win_string_end+1)+1)<<endl;
														
							fed=edlibAlign(&(pattern.c_str())[y1+1],x2-1-y1,&(cur_string.c_str())[win_string_end+1],next_win_string_start-1-(win_string_end), edlibDefaultAlignConfig());
							ub+=fed.editDistance;
							//cout<<"3b. "<<fed.editDistance;
						}
						else //one of the two has overlap or both have overlap
						{
							if(y2-x1 >= next_win_string_end-win_string_start)  
							{
								
								//cout<<"3c. "<<(y2-x1+1)-(next_win_string_end-win_string_start+1)<<endl;
								ub+=(y2-x1)-(next_win_string_end-win_string_start);
							}
							else
							{
								ub+=next_win_string_end-win_string_start-y2+x1;
							}
						}
					}
				}
				else //not in the beginning and teleutaio it_next == e_is.end() 
				{
					
                    INT y1=*it2;       //pat window end                                       
					INT x1=y1-ell+1; //pat window start

                    INT win_string_end=(*it_string_end);
                    INT win_string_start=win_string_end-ell+1;
                                                           
					if(y1!=pattern.length()-1 && win_string_end!=cur_string.length()-1) //pattern window end is not the last character and string window end is not the last character
					{
						//cout<<"4a. "<<pattern.substr(y1+1,pattern.length()-1-(y1+1)+1)<<" "<<cur_string.substr(win_string_end+1,cur_string.length()-1-(win_string_end+1)+1)<<endl;
						
						fed=edlibAlign(&(pattern.c_str())[y1+1],pattern.length()-(y1+1),&(cur_string.c_str())[win_string_end+1],cur_string.length()-(win_string_end+1), edlibDefaultAlignConfig());
						ub+=fed.editDistance;
						//cout<<"4a. "<<fed.editDistance;
					}
					else if(y1==pattern.length()-1)
					{
						
						//cout<<"4b. "<<cur_string.substr(win_string_end+1,cur_string.length()-1-(win_string_end+1)+1)<<endl;
						ub+=cur_string.length()-(win_string_end+1);
					}
					else if(win_string_end==cur_string.length()-1)
					{
						ub+=pattern.length()-(y1+1);
						//cout<<"4c. "<<pattern.substr(y1+1,((pattern.length()-1-(y1+1))+1))<<endl;
						
					}
				}
			}
			//here ub contains the upper-bound of edit distance for the string
	   		//cout<<"ub="<<ub<<endl;
			
			double score_prod=ub;
		
			if(pq_ub.size()<K) //if I have seen <K elements
			{
				pq_ub.push(make_pair(score_prod,it));				
			}
			else
			{				
				
				if(pq_ub.top().first>(score_prod))
				{		
					pq_ub.pop();
					pq_ub.push(make_pair(score_prod,it));
					
				}			
			}		
			
		}		
		//end step 3
		
		total_duration += chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now()-start_time).count();	
		
		
		cout<<"Pattern: "<<pattern<<endl;
		list<INT> est_ids;
		vector<INT> edit_distances_of_answers;   
		while(!pq_ub.empty())
		{
			pair<double,INT> p=pq_ub.top();
			
			pair<INT, vector<unsigned char> >  v=all_strings[p.second];
			string cur_string(v.second.begin(),v.second.end());
			
			const char* pat_c=pattern.c_str();
			const char* cur_c=cur_string.c_str();
			
			EdlibAlignResult fed=edlibAlign(pat_c,strlen(pat_c),cur_c,strlen(cur_c), edlibDefaultAlignConfig());
			
			INT ed = fed.editDistance;
			edit_distances_of_answers.push_back(ed);
			
			est_ids.push_front(p.second);
			
			cout<<"Ub: "<<p.first<<" Edit distance: "<<ed<<"\n";  //we multiply with -1 to get the real ub														
			cout<<"Id: "<< v.first<<"     String: "<<cur_string<<"\n";
			
			sum_of_edit_distance_of_answers+=ed;
			pq_ub.pop();
		}
		if(est_ids.size()<K)
		{
			cout<<"Try a smaller K.\n";
			exit(-1);
		}	
 	} 
 
	
 	cout << "Elapsed total time (ms) : "<<total_duration<<endl;
 	cout<< "Avg time per query (ms) : "<<(double)total_duration/all_patterns.size()<<endl;
   

 	free ( RSA );
  	free ( RLCP );
  	free ( LSA );
  	free ( LLCP );
  	free(ids);
	free(num_anchors_per_string);		
	invRSA.clear();
	    
	return 0;
}

